using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Breeze.AtmosphereModels: thermodynamic_density
using Oceananigans.BoundaryConditions: fill_halo_regions!
const GRAV=.8; const CP=1.; const RGAS=.4; const THETA0=1.
function parameters(c)
    (;nx=c["nx"],ny=c["ny"],nz=c["nz"],dx=c["length"]/c["nx"],dy=c["length"]/c["ny"],dz=c["height"]/c["nz"],height=c["height"],dc=c["convective_depth"],f=c["f"],nu=c["nu"],kappa=c["kappa"],flux=c["bottom_theta_flux"],cap_delta=c["cap_delta_theta"],tau=c["cap_relaxation_time"])
end
@inline zcell(k,p)=(k-.5)*p.dz
@inline function cap_shape(k,p)
    x=clamp((zcell(k,p)-p.dc)/(p.height-p.dc),0.,1.)
    x*x*(3-2x)
end
@inline cap_target(k,p)=THETA0+p.cap_delta*cap_shape(k,p)
@inline function cooling(i,j,k,g,clock,F,p)
    @inbounds theta=F.T[i,j,k]/p.Tref[i,j,k]
    @inbounds return -p.rho[i,j,k]*cap_shape(k,p)*(theta-cap_target(k,p))/p.tau
end
function build_model(c;arch=GPU(),rest=false)
    p=parameters(c)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,c["length"]),y=(0,c["length"]),z=(0,p.height),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=RGAS,dry_air_molar_mass=1.,dry_air_heat_capacity=CP,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    ref=ReferenceState(grid,constants;base_pressure=.4,standard_pressure=.4,potential_temperature=THETA0)
    p=merge(p,(rho=ref.density,Tref=ref.temperature))
    dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64)
    bc=FieldBoundaryConditions(bottom=FluxBoundaryCondition(rest ? 0. : p.flux),top=FluxBoundaryCondition(0.))
    cf=rest ? NamedTuple() : (ρθ=Forcing(cooling;discrete_form=true,parameters=p),)
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=FPlane(f=p.f),closure=ScalarDiffusivity(ν=p.nu,κ=(ρθ=p.kappa,ρqᵛ=0.,dye=0.)),tracers=(:dye,),boundary_conditions=(ρθ=bc,),forcing=cf,momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))))
    rng=MersenneTwister(c["seed"]);noise=rand(rng,p.nx,p.ny,p.nz) .- .5;noise .-= mean(noise;dims=(1,2))
    theta=[rest ? THETA0 : cap_target(k,p)+c["initial_delta_theta"]*max(0.,1-zcell(k,p)/p.dc)+c["noise_amplitude"]*noise[i,j,k]*sinpi(clamp(zcell(k,p)/p.dc,0.,1.)) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model;θ=theta,qᵗ=0.,u=0.,v=0.,w=0.)
    rho=Array(interior(ref.density))
    dye=[.1+.05cospi(2*(i-.5)/p.nx)*cospi(2*(j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model.tracers.dye,dye .* rho)
    return model,p
end
function write_toml(path,d)
    open(path,"w") do io;TOML.print(io,d);end
end
mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)
function arrays(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .* (uf .+ circshift(uf,(-1,0,0)));v=.5 .* (vf .+ circshift(vf,(0,-1,0)));w=.5 .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
    rho=Array(interior(m.dynamics.reference_state.density));theta=Array(interior(thermodynamic_density(m.formulation))) ./ rho
    T=Array(interior(m.temperature));mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    (;uf,vf,wf,u,v,w,rho,theta,T,mx,my,mz)
end
function run_case(c,out;arch=GPU(),rest=false,pickup=false)
    isdir(out)&&!isempty(readdir(out))&&error("Refuse overwrite $out")
    mkpath(out);m,p=build_model(c;arch,rest)
    z=[zcell(k,p) for k=1:p.nz];planes=[argmin(abs.(z .- zz)) for zz in (.2,.6,.9)];cz=findall(z .< p.dc)
    cap=min(c["max_dt"],.08/(max(p.nu,p.kappa)*(p.dx^-2+p.dy^-2+p.dz^-2)),.1/max(abs(p.f),1e-20),.1p.tau)
    sim=Simulation(m;Δt=c["initial_dt"],stop_time=c["stop_time"],wall_time_limit=c["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.3,max_Δt=cap)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(c["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    meta=merge(c,Dict("status"=>"RUNNING","architecture"=>string(typeof(arch)),"model_sha256"=>bytes2hex(sha256(read(@__FILE__))),"video_heights"=>z[planes],"time_units"=>"nondimensional; not Earth days","dt_cap"=>cap,"physics"=>"dry anelastic, fixed adiabatic density, constant f, native density-weighted viscosity and potential-temperature diffusion; no viscous heating; cap theta relaxation; no top thermal flux","budget"=>"Conservative potential-temperature content budget monitored; not a closed total-energy/vorticity budget certificate","continuation"=>"Manual review required; no branch launch or auto-extension"))
    write_toml(joinpath(out,"metadata.toml"),meta)
    file=joinpath(out,"diagnostics.csv")
    open(file,"w") do io;println(io,"time,iteration,dt,wall_seconds,theta_content,theta_source,budget_residual_over_flux,urms_conv,wrms_conv,Ro_depth,Re_depth,turnovers_depth,lambda_w,lambda_cells,Ro_gradient,max_mach,max_relative_theta,max_mass_divergence,wall_w_max,barotropic_ke,convective_theta_flux_mid,cap_min_dtheta_dz");end
    started=time();old=Ref((t=NaN,mtheta=NaN,source=NaN,urms=0.,turnovers=0.))
    function output(sim)
        a=arrays(sim.model);t=Float64(sim.model.clock.time);vol=c["length"]^2*p.height
        mtheta=sum(a.rho .* a.theta)*p.dx*p.dy*p.dz
        rates=reshape([cap_shape(k,p)/p.tau for k=1:p.nz],1,1,:);targets=reshape([cap_target(k,p) for k=1:p.nz],1,1,:)
        source=rest ? 0. : c["length"]^2*p.flux-sum(a.rho .* (a.theta .- targets) .* rates)*p.dx*p.dy*p.dz
        prev=old[];residual=isfinite(prev.t)&&t>prev.t ? ((mtheta-prev.mtheta)/(t-prev.t)-.5*(source+prev.source))/max(c["length"]^2*p.flux,1e-30) : 0.
        urms=sqrt(mean(a.u[:,:,cz].^2 .+ a.v[:,:,cz].^2 .+ a.w[:,:,cz].^2));wrms=sqrt(mean(a.w[:,:,cz].^2))
        turnovers=prev.turnovers+(isfinite(prev.t) ? .5*(urms+prev.urms)*(t-prev.t)/p.dc : 0.)
        old[]=(t=t,mtheta=mtheta,source=source,urms=urms,turnovers=turnovers)
        dwx=(circshift(a.w,(-1,0,0)).-a.w)./p.dx;dwy=(circshift(a.w,(0,-1,0)).-a.w)./p.dy
        grad=sqrt(mean(dwx[:,:,cz].^2 .+ dwy[:,:,cz].^2));ell=wrms/max(grad,1e-30);lambda=2pi*ell
        mach=maximum(sqrt.(a.u.^2 .+ a.v.^2 .+ a.w.^2)./sqrt.((5/3)*RGAS .* a.T))
        div=(circshift(a.mx,(-1,0,0)).-a.mx)./p.dx .+(circshift(a.my,(0,-1,0)).-a.my)./p.dy .+(a.mz[:,:,2:end].-a.mz[:,:,1:end-1])./p.dz
        md=maximum(abs,div);wallw=maximum(abs,a.wf[:,:,[1,end]])
        ub=dropdims(sum(a.rho .* a.u;dims=3)./sum(a.rho;dims=3);dims=3);vb=dropdims(sum(a.rho .* a.v;dims=3)./sum(a.rho;dims=3);dims=3)
        thmean=mean(a.theta;dims=(1,2));convflux=vec(mean(a.rho .* a.w .* (a.theta .- thmean);dims=(1,2)))
        profile=vec(thmean);capgrad=minimum(diff(profile)[findall(z[1:end-1] .> .85)]./p.dz)
        vals=(t,sim.model.clock.iteration,sim.Δt,time()-started,mtheta,source,residual,urms,wrms,urms/max(abs(p.f)*p.dc,1e-30),urms*p.dc/p.nu,turnovers,lambda,lambda/min(p.dx,p.dy),grad/max(abs(p.f),1e-30),mach,maximum(abs,a.theta .- 1),md,wallw,.5mean(ub.^2 .+ vb.^2),convflux[argmin(abs.(z .- p.dc/2))],capgrad)
        open(file,"a") do io;println(io,join(vals,','));end
        all(isfinite,vals)||error("Nonfinite diagnostic")
        minimum(a.T)>0||error("Nonpositive temperature")
        maximum(abs,a.theta.-1)<.1||error("Reference-state departure exceeds 10%; anelastic pilot invalid")
        md<1e-8||error("Mass-flux projection residual exceeds 1e-8: $md")
        mach<.3||error("Mach exceeds 0.3; review approximation")
        wallw<1e-10||error("Wall impermeability failed")
        (rest || t<5 || abs(residual)<.05)||error("Potential-temperature budget discrepancy >5% input; do not auto-continue: $residual")
        if isapprox(t/c["movie_interval"],round(t/c["movie_interval"]);atol=1e-7)||t>=c["stop_time"]
            zf=(a.vf.-circshift(a.vf,(1,0,0)))./p.dx .-(a.uf.-circshift(a.uf,(0,1,0)))./p.dy
            zeta=.25 .* (zf.+circshift(zf,(-1,0,0)).+circshift(zf,(0,-1,0)).+circshift(zf,(-1,-1,0)))
            dye=Array(interior(m.tracers.dye))./a.rho
            jldsave(joinpath(out,@sprintf("frame_%08d.jld2",round(Int,100t)));t,z,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),theta_slices=Float32.(a.theta[:,:,planes]),theta_profile=profile,w_rms_profile=vec(sqrt.(mean(a.w.^2;dims=(1,2)))),convective_theta_flux_profile=convflux,ubar=Float32.(ub),vbar=Float32.(vb),f=p.f)
        end
        @printf("t=%.3f step=%d dt=%.4g wall=%.1f turns=%.3f urms=%.4g flux/F=%.3g div=%.2g budget=%.3g\n",t,m.clock.iteration,sim.Δt,time()-started,turnovers,urms,convflux[argmin(abs.(z.-p.dc/2))]/p.flux,md,residual);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(c["diagnostic_interval"]);name=:output)
    if pickup!==false;set!(sim;checkpoint=pickup);sim.Δt=saved.dt;end
    run!(sim;checkpoint_at_end=true)
    meta["status"]=m.clock.time>=c["stop_time"]-1e-8 ? "TARGET_REACHED_REVIEW_REQUIRED" : "WALL_LIMIT_REVIEW_REQUIRED"
    meta["final_time"]=Float64(m.clock.time);meta["wall_elapsed"]=time()-started;meta["turnovers_depth"]=old[].turnovers
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
