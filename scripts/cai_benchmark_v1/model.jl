using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Oceananigans.Coriolis: ConstantCartesianCoriolis
using Oceananigans.TurbulenceClosures: Smagorinsky
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ

# Breeze passes momentum fields to Coriolis, while the generic Cartesian
# operator expects velocity-style field names. This wrapper only renames fields;
# interpolation, signs, coefficients and the energy-conserving stencil are native.
struct MomentumCoriolis{C} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    rotation::C
end
@inline momentum_aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))

# Conservative temperature-conduction flux divergence. Boundary fluxes are part
# of this operator, not duplicated in prognostic theta boundary conditions.
@inline function conduction(i,j,k,T,p)
    im=ifelse(i==1,p.nx,i-1);ip=ifelse(i==p.nx,1,i+1)
    jm=ifelse(j==1,p.ny,j-1);jp=ifelse(j==p.ny,1,j+1)
    km=max(1,k-1);kp=min(p.nz,k+1)
    @inbounds begin
        tc=T[i,j,k]
        horizontal=p.conductivity*((T[ip,j,k]-2tc+T[im,j,k])/p.dx^2+(T[i,jp,k]-2tc+T[i,jm,k])/p.dy^2)
        qbottom=ifelse(k==1,p.flux,-p.conductivity*(tc-T[i,j,km])/p.dz)
        qtop=ifelse(k==p.nz,-2p.conductivity*(p.top_temperature-tc)/p.dz,-p.conductivity*(T[i,j,kp]-tc)/p.dz)
    end
    return horizontal+(qbottom-qtop)/p.dz
end
@inline function energy_source(i,j,k,grid,clock,fields,p)
    sigma2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,grid,fields.u,fields.v,fields.w)
    nu=p.smag_length_squared*sqrt(2sigma2)
    @inbounds heating=2fields.ρᵈ[i,j,k]*nu*sigma2
    return conduction(i,j,k,fields.T,p)+heating
end
function write_toml(path,d)
    open(path,"w") do io;TOML.print(io,d);end
end
function build_model(cfg;arch=GPU(), perturb=true)
    nx,ny,nz=cfg["nx"],cfg["ny"],cfg["nz"]
    dx=4/nx;dy=4/ny;dz=1/nz
    grid=RectilinearGrid(arch,Float64;size=(nx,ny,nz),x=(0,4),y=(0,4),z=(0,1),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=8.,energy_reference_temperature=0.)
    theta_ref(z)=(5-4z)^0.2
    dynamics=CompressibleDynamics(SplitExplicitTimeDiscretization();standard_pressure=1.,base_pressure=25.,reference_potential_temperature=theta_ref)
    # Cai Eq A7: cd² sqrt(dx dy) dz sqrt(2 S:S). Match its filter length
    # with native Smagorinsky's (dx dy dz)^(2/3) convention.
    ell2=0.28^2*sqrt(dx*dy)*dz
    coefficient=sqrt(ell2)/(dx*dy*dz)^(1/3)
    closure=Smagorinsky(;coefficient,Pr=Inf)
    p=(;nx,ny,nz,dx,dy,dz,conductivity=cfg["flux"]/4,flux=cfg["flux"],top_temperature=1.,smag_length_squared=ell2)
    coriolis=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=cfg["omega"],latitude=cfg["latitude"]))
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis,closure,tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(;ρE=Forcing(energy_source;discrete_form=true,parameters=p)))
    z=((1:nz) .- 0.5)/nz
    rng=MersenneTwister(cfg["seed"]);noise=rand(rng,nx,ny,nz) .- 0.5;noise .-=mean(noise,dims=(1,2))
    theta=[theta_ref(z[k])*(1+(perturb ? cfg["perturbation"] : 0)*noise[i,j,k]*sinpi(z[k])) for i=1:nx,j=1:ny,k=1:nz]
    ref=model.dynamics.reference_state
    set!(model;θ=theta,ρ=ref.density,qᵗ=0,u=0,v=0,w=0)
    rho=Array(interior(model.dynamics.dry_density))
    dye=[.8exp(-((sinpi(((i-.5)*dx-1)/4)*4/pi)^2+(sinpi(((j-.5)*dy-2)/4)*4/pi)^2)/(.45^2)) for i=1:nx,j=1:ny,k=1:nz]
    set!(model.tracers.dye,dye .* rho)
    return model,p
end
mutable struct SavedStep
    value::Float64
end
(s::SavedStep)(sim)=(s.value=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;value=s.value)
Oceananigans.restore_prognostic_state!(s::SavedStep,state::NamedTuple)=(s.value=state.value;s)
function state_arrays(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .* (uf .+ circshift(uf,(-1,0,0)));v=.5 .* (vf .+ circshift(vf,(0,-1,0)));w=.5 .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
    rho=Array(interior(m.dynamics.dry_density));T=Array(interior(m.temperature));nu=Array(interior(m.closure_fields.νₑ))
    return (;u,v,w,rho,T,nu,uf,vf,wf)
end
function run_case(cfg,out;arch=GPU(),pickup=false)
    isdir(out)&&!isempty(readdir(out))&&error("Refuse overwrite $out")
    mkpath(out);model,p=build_model(cfg;arch)
    sim=Simulation(model;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cfg["max_dt"])
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(model;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    started=time();z=((1:p.nz) .- 0.5)*p.dz;planes=[argmin(abs.(z .- h)) for h in (.2,.5,.8)]
    f=2cfg["omega"]*sind(cfg["latitude"])
    meta=merge(cfg,Dict("status"=>"RUNNING","equations"=>"dry fully compressible ideal gas gamma=5/3","f"=>f,"fy"=>2cfg["omega"]*cosd(cfg["latitude"]),"beta"=>0.,"heights"=>z[planes],"conductivity"=>p.conductivity,"source_sha256"=>bytes2hex(sha256(read(@__FILE__))),"breeze_commit"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b","paper"=>"https://arxiv.org/pdf/2111.11003","replication_status"=>"preliminary lower-resolution adaptation; full local SGS stress rather than paper's simplified fluctuating-viscosity treatment; finite-volume WENO rather than spectral-horizontal scheme; explicit conduction; native acoustic divergence damping","initial_horizontal_rms"=>0.,"initial_large_vortex"=>false,"spatial_dimensions"=>[4.,4.,1.]))
    write_toml(joinpath(out,"metadata.toml"),meta)
    diagnostic=joinpath(out,"diagnostics.csv")
    open(diagnostic,"w") do io
        println(io,"time,iteration,dt,wall_seconds,mass,total_energy,rms_velocity,rms_w,max_mach,mean_nu,top_heat_flux,bottom_heat_flux,barotropic_ke,skew_zeta_lower,skew_zeta_mid,skew_zeta_upper,min_density,min_temperature,wall_w_max")
    end
    function output(sim)
        a=state_arrays(sim.model);t=sim.model.clock.time
        mass=mean(a.rho)*16
        total=mean(a.rho .* (1.5 .* a.T .+ 8reshape(z,1,1,:) .+ .5 .* (a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2)))*16
        rms=sqrt(mean(a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2));wrms=sqrt(mean(a.w .^ 2));mach=maximum(sqrt.(a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2) ./ sqrt.((5/3) .* a.T))
        zetaf=(a.vf .- circshift(a.vf,(1,0,0))) ./ p.dx .- (a.uf .- circshift(a.uf,(0,1,0))) ./ p.dy
        zeta=.25 .* (zetaf .+ circshift(zetaf,(-1,0,0)) .+ circshift(zetaf,(0,-1,0)) .+ circshift(zetaf,(-1,-1,0)))
        skew=[mean(zeta[:,:,k] .^ 3)/max(mean(zeta[:,:,k] .^ 2)^1.5,1e-300) for k in planes]
        weights=a.rho ./ sum(a.rho,dims=3);ub=dropdims(sum(weights .* a.u,dims=3),dims=3);vb=dropdims(sum(weights .* a.v,dims=3),dims=3)
        qtop=mean(2p.conductivity .* (a.T[:,:,end] .- 1) ./ p.dz)
        vals=(t,sim.model.clock.iteration,sim.Δt,time()-started,mass,total,rms,wrms,mach,mean(a.nu),qtop,p.flux,.5mean(ub .^ 2 .+ vb .^ 2),skew...,minimum(a.rho),minimum(a.T),maximum(abs,a.wf[:,:,[1,end]]))
        all(isfinite,vals)||error("nonfinite diagnostic");minimum(a.rho)>0&&minimum(a.T)>0||error("Nonpositive thermodynamics")
        open(diagnostic,"a") do io;println(io,join(vals,","));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            dye=Array(interior(sim.model.tracers.dye)) ./ a.rho
            filename=joinpath(out,@sprintf("frame_%05d.jld2",round(Int,t/cfg["movie_interval"])))
            isfile(filename)&&error("Duplicate frame")
            jldsave(filename;t,z,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),density_profile=vec(mean(a.rho,dims=(1,2))),temperature_profile=vec(mean(a.T,dims=(1,2))),w_rms_profile=vec(sqrt.(mean(a.w .^ 2,dims=(1,2)))),f)
        end
        if isapprox(t/cfg["checkpoint_interval"],round(t/cfg["checkpoint_interval"]);atol=1e-7) || t>=cfg["stop_time"]
            mkpath(joinpath(out,"fields3d"))
            jldsave(joinpath(out,"fields3d",@sprintf("state_%07d.jld2",round(Int,100t)));t,z,
                u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),
                density=Float32.(a.rho),temperature=Float32.(a.T),viscosity=Float32.(a.nu))
        end
        @printf("t=%.3f iteration=%d dt=%.4g wall=%.1fs rms=%.4g Mach=%.3g qtop=%.5g\n",t,sim.model.clock.iteration,sim.Δt,time()-started,rms,mach,qtop);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    if pickup!==false;set!(sim;checkpoint=pickup);sim.Δt=saved.value;end
    run!(sim;checkpoint_at_end=true)
    meta["status"]=model.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    meta["final_time"]=model.clock.time;meta["wall_elapsed"]=time()-started;meta["iterations"]=model.clock.iteration
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
