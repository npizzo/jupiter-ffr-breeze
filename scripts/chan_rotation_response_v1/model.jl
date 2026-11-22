using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
using Breeze.CompressibleEquations: NoDivergenceDamping

# Chan--Mayr dimensionless units: R=Ttop=rhotop=ptop=Lz=1, g=12.
using Oceananigans.Coriolis: ConstantCartesianCoriolis
using Oceananigans.TurbulenceClosures: Smagorinsky
struct MomentumCoriolis{C} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    rotation::C
end
@inline momentum_aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
const CV=1.5; const RGAS=1.; const GRAV=12.; const DC=.78
function reference(z)
    ti=1+GRAV*(1-DC)/25
    if z>=DC
        t=1+GRAV*(1-z)/25; p=t^25
    else
        t=ti+GRAV*(DC-z)/2.5; p=ti^25*(t/ti)^2.5
    end
    return (;T=t,rho=p/t,p,theta=t/p^.4)
end
function parameters(cfg)
    return (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],dx=4/cfg["nx"],dy=4/cfg["ny"],dz=1/cfg["nz"],
        omega=cfg["omega"],f=2cfg["omega"]*sind(cfg["latitude"]),pr=cfg["pr"],kc=.0004*cfg["bottom_flux"]/(GRAV/2.5),flux=cfg["bottom_flux"],
        tau=cfg["cooling_time"],heat=cfg["heat"],tt=1.,
        target_omega=cfg["target_omega"],ramp_start=cfg["ramp_start"],ramp_duration=cfg["ramp_duration"],
        rotation=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=cfg["omega"],latitude=cfg["latitude"])))
end
@inline zcell(k,p)=(k-.5)*p.dz
@inline function cool_rate(k,p)
    q=clamp((zcell(k,p)-DC)/(1-DC),0.,1.)
    return q*q*(3-2q)/p.tau
end
@inline function conduction(i,j,k,T,p)
    im=ifelse(i==1,p.nx,i-1);ip=ifelse(i==p.nx,1,i+1)
    jm=ifelse(j==1,p.ny,j-1);jp=ifelse(j==p.ny,1,j+1)
    km=max(1,k-1);kp=min(p.nz,k+1)
    @inbounds begin
        tc=T[i,j,k]
        horizontal=p.kc*((T[ip,j,k]-2tc+T[im,j,k])/p.dx^2+(T[i,jp,k]-2tc+T[i,jm,k])/p.dy^2)
        qb=ifelse(k==1,p.flux,-p.kc*(tc-T[i,j,km])/p.dz)
        qt=ifelse(k==p.nz,-2p.kc*(p.tt-tc)/p.dz,-p.kc*(T[i,j,kp]-tc)/p.dz)
    end
    return horizontal+(qb-qt)/p.dz
end
@inline entropy(i,j,k,F)=CV*log(F.T[i,j,k])-log(F.ρᵈ[i,j,k])
@inline conductance(i,j,k,F,p)=F.ρᵈ[i,j,k]*F.T[i,j,k]*p.nu[i,j,k]/p.pr
@inline function entropy_flux(i,j,k,ii,jj,kk,F,p,delta)
    @inbounds return -.5*(conductance(i,j,k,F,p)+conductance(ii,jj,kk,F,p))*(entropy(ii,jj,kk,F)-entropy(i,j,k,F))/delta
end
@inline function sgs_heating(i,j,k,F,p)
    im=ifelse(i==1,p.nx,i-1);ip=ifelse(i==p.nx,1,i+1)
    jm=ifelse(j==1,p.ny,j-1);jp=ifelse(j==p.ny,1,j+1)
    qx=(entropy_flux(im,j,k,i,j,k,F,p,p.dx)-entropy_flux(i,j,k,ip,j,k,F,p,p.dx))/p.dx
    qy=(entropy_flux(i,jm,k,i,j,k,F,p,p.dy)-entropy_flux(i,j,k,i,jp,k,F,p,p.dy))/p.dy
    qb=k==1 ? 0. : entropy_flux(i,j,k-1,i,j,k,F,p,p.dz)
    qt=k==p.nz ? 0. : entropy_flux(i,j,k,i,j,k+1,F,p,p.dz)
    return qx+qy+(qb-qt)/p.dz
end
@inline function divergence(i,j,k,F,p)
    @inbounds return (F.u[i+1,j,k]-F.u[i,j,k])/p.dx+(F.v[i,j+1,k]-F.v[i,j,k])/p.dy+(F.w[i,j,k+1]-F.w[i,j,k])/p.dz
end
@inline trace_stress(i,j,k,g,F,p)=-2p.nu[i,j,k]/3*F.ρᵈ[i,j,k]*divergence(i,j,k,F,p)
# Parameter continuation of the Coriolis coefficient, not a physically accelerated frame.
# Native base rotation plus its time-dependent linear increment preserves the native
# momentum interpolation and includes BOTH horizontal and vertical rotation components.
@inline function omega_at(t,p)
    q=clamp((t-p.ramp_start)/p.ramp_duration,0.,1.)
    return p.omega+(p.target_omega-p.omega)*q*q*(3-2q)
end
@inline rotation_fraction(t,p)=omega_at(t,p)/p.omega-1
@inline rotation_x(i,j,k,g,c,F,p)=-rotation_fraction(c.time,p)*Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,p.rotation,F)
@inline rotation_y(i,j,k,g,c,F,p)=-rotation_fraction(c.time,p)*Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,p.rotation,F)
@inline rotation_z(i,j,k,g,c,F,p)=-rotation_fraction(c.time,p)*Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,p.rotation,F)
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)+rotation_x(i,j,k,g,c,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)+rotation_y(i,j,k,g,c,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)+rotation_z(i,j,k,g,c,F,p)
@inline function heat_source(i,j,k,g,c,F,p)
    div=divergence(i,j,k,F,p)
    strain2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3
    # Paper Eq. 3 is -rho*(T-Ttop)/tau, WITHOUT a cv factor.
    @inbounds return 2p.nu[i,j,k]*F.ρᵈ[i,j,k]*strain2+sgs_heating(i,j,k,F,p)+p.heat*(conduction(i,j,k,F.T,p)-F.ρᵈ[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
end
function build_model(cfg;arch=GPU(),perturb=true)
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,4),y=(0,4),z=(0,1),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    # Bind the forcing to the actual cached viscosity used by the native closure.
    nu=CenterField(grid);p=merge(p,(;nu))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    th(z)=reference(z).theta
    dynamics=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping=NoDivergenceDamping());standard_pressure=1.,base_pressure=reference(0.).p,reference_potential_temperature=th)
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    rot=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=p.omega,latitude=cfg["latitude"]))
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=rot,closure=Smagorinsky(coefficient=cfg["smagorinsky_coefficient"],Pr=Inf),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    model.closure_fields=(νₑ=nu,)
    rng=MersenneTwister(cfg["seed"])
    noise=rand(rng,p.nx,p.ny,p.nz) .- .5;noise .-= mean(noise,dims=(1,2))
    theta=[th(zcell(k,p))*(1+(perturb ? cfg["perturbation"] : 0.)*noise[i,j,k]*sinpi(clamp(zcell(k,p)/DC,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model;θ=theta,ρ=model.dynamics.reference_state.density,qᵗ=0,u=0,v=0,w=0)
    seed_info=apply_seed!(model,p,cfg;perturb)
    p=merge(p,(;seed_info))
    rho=Array(interior(model.dynamics.dry_density))
    dye=[.1+.05cospi(2*(i-.5)*p.dx/4)*cospi(2*(j-.5)*p.dy/4) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
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
    u=0.5  .* (uf .+ circshift(uf,(-1,0,0)));v=0.5  .* (vf .+ circshift(vf,(0,-1,0)));w=0.5  .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
    rho=Array(interior(m.dynamics.dry_density));T=Array(interior(m.temperature))
    (;uf,vf,wf,u,v,w,rho,T)
end
function energy_budget_notice(residual,t)
    if t>=1 && abs(residual)>=.1
        @warn "Energy residual exceeds 10% of heating input; continuing exploratory integration and retaining the budget diagnostic" time=t residual_over_heating=residual maxlog=1
    end
    nothing
end
function run_case(cfg,out;arch=GPU(),pickup=false)
    isdir(out)&&!isempty(readdir(out))&&error("Refuse overwrite $out")
    mkpath(out);m,p=build_model(cfg;arch)
    rho0=Array(interior(m.dynamics.dry_density));z=[zcell(k,p) for k=1:p.nz]
    # Explicit diffusion, cooling and rotation caps supplement the advective CFL.
    thermal_dt=.08minimum(rho0)*CV/(p.kc*(p.dx^-2+p.dy^-2+p.dz^-2))
    viscous_dt=Inf
    cap=min(cfg["max_dt"],thermal_dt,viscous_dt,.15p.tau,.1/max(2max(abs(p.omega),abs(p.target_omega)),1e-20))
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    # SGS coefficient is time dependent; enforce thermal diffusion stability each step.
    function diffusion_limit(sim)
        numax=maximum(interior(sim.model.closure_fields.νₑ))
        maxdiff=max(4numax/3, (2.5/CV)*numax/p.pr)
        sim.Δt=min(sim.Δt,.08/(max(maxdiff,1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diffusion_limit,IterationInterval(1);name=:diffusion_limit)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    planes=[argmin(abs.(z .- zz)) for zz in (.2,.5,.9)];cz=findall(0 .< z .< DC)
    meta=merge(cfg,Dict("status"=>"RUNNING","f"=>p.f,"conductivity"=>p.kc,"dt_cap"=>cap,"video_heights"=>z[planes],"initial_density_contrast"=>maximum(rho0)/minimum(rho0),"model_sha256"=>bytes2hex(sha256(read(@__FILE__))),"paper"=>"https://doi.org/10.1016/j.epsl.2013.03.046","fidelity"=>"Chan--Mayr physical reconstruction, reduced grid; Cs=.16 cubic-volume filter, entropy diffusion, smoothstep cooling explicitly chosen; potential-temperature WENO scheme is not paper conservative total-energy discretization","budget_status"=>"Exploratory finite-amplitude evolution: 10%-of-heating residual is diagnostic only; no closed energy budget or equilibrated sign-transition claim","time_unit"=>"box-height / top isothermal sound speed, not Earth days"))
    meta["base_physics_from"]="chan_eddies_v4"
    meta["only_physical_change"]="Time-dependent full Coriolis coefficient in ramp branch; no Euler force, beta, shear, or heating change"
    meta["vorticity_display_reference_f"]=p.f
    meta["energy_stop_policy"]="Warn at 10% of heating input; do not stop on this residual alone. Retain finite-state, positive rho/T, mass, Mach and timestep safeguards."
    meta["initialization"]=p.seed_info
    meta["experiment_scope"]="Matched response of an existing mixed-sign vortex state; no spontaneous formation or equilibrium claim."
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,iteration,dt,wall_seconds,mass,total_energy,urms_conv,Ro,Re,max_mach,barotropic_ke,skew_lower,skew_mid,skew_upper,top_heat_flux,cooling_power,min_density,min_temperature,wall_w_max,nu_mean,external_power,energy_residual_over_input,centered_energy,centered_residual_over_input,kinetic_energy,cumulative_external_heat_since_restart,cumulative_energy_defect_since_restart,energy_defect_over_restart_ke,omega_actual,f_actual");end
    started=time();initial_mass=mean(rho0)*16
    if !isempty(p.seed_info);write_toml(joinpath(out,"INITIALIZATION.toml"),p.seed_info);end
    tracker=Ref((initial_energy=NaN,initial_ke=NaN,heat=0.))
    previous=Ref((time=NaN,energy=NaN,power=NaN,centered=NaN))
    function output(sim)
        a=arrays(sim.model);t=Float64(sim.model.clock.time)
        mass=mean(a.rho)*16;centered_energy=mean(a.rho .* (CV .* a.T .+ GRAV .* reshape(z,1,1,:) .+ 0.5  .* (a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2)))*16
        rhox=.5 .* (a.rho .+ circshift(a.rho,(1,0,0)))
        rhoy=.5 .* (a.rho .+ circshift(a.rho,(0,1,0)))
        rhoz=.5 .* (a.rho[:,:,1:end-1] .+ a.rho[:,:,2:end])
        ke=(sum(rhox .* a.uf.^2)+sum(rhoy .* a.vf.^2)+sum(rhoz .* a.wf[:,:,2:end-1].^2))*.5*p.dx*p.dy*p.dz
        energy=mean(a.rho .* (CV .* a.T .+ GRAV .* reshape(z,1,1,:)))*16+ke
        urms=sqrt(mean(a.u[:,:,cz] .^ 2 .+ a.v[:,:,cz] .^ 2 .+ a.w[:,:,cz] .^ 2));omega_now=omega_at(t,p);co=urms/max(2abs(omega_now)*DC,1e-30)
        numean=mean(Array(interior(sim.model.closure_fields.νₑ))[:,:,cz]);re=urms*DC/max(numean,1e-30)
        mach=maximum(sqrt.(a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2) ./ sqrt.((5/3)*RGAS .* a.T))
        zf=(a.vf .- circshift(a.vf,(1,0,0))) ./ p.dx .- (a.uf .- circshift(a.uf,(0,1,0))) ./ p.dy
        zeta=0.25  .* (zf .+ circshift(zf,(-1,0,0)) .+ circshift(zf,(0,-1,0)) .+ circshift(zf,(-1,-1,0)))
        skew=[mean(zeta[:,:,k] .^ 3)/max(mean(zeta[:,:,k] .^ 2)^1.5,1e-30) for k in planes]
        weight=a.rho ./ sum(a.rho;dims=3);ub=dropdims(sum(weight .* a.u;dims=3);dims=3);vb=dropdims(sum(weight .* a.v;dims=3);dims=3)
        qtop=mean(2p.kc .* (a.T[:,:,end] .- p.tt) ./ p.dz)
        cooling=mean(a.rho .* (a.T .- p.tt) .* reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))*16
        power=p.heat*(16*(p.flux-qtop)-cooling)
        old=previous[]
        residual=isfinite(old.time)&&t>old.time ? ((energy-old.energy)/(t-old.time)-.5*(power+old.power))/max(16p.flux,1e-30) : 0.
        centered_residual=isfinite(old.time)&&t>old.time ? ((centered_energy-old.centered)/(t-old.time)-.5*(power+old.power))/max(16p.flux,1e-30) : 0.
        budget=tracker[]
        if !isfinite(budget.initial_energy)
            budget=(initial_energy=energy,initial_ke=ke,heat=0.)
        elseif isfinite(old.time) && t>old.time
            budget=(initial_energy=budget.initial_energy,initial_ke=budget.initial_ke,heat=budget.heat+.5*(power+old.power)*(t-old.time))
        end
        tracker[]=budget
        cumulative_defect=energy-budget.initial_energy-budget.heat
        previous[]=(time=t,energy=energy,power=power,centered=centered_energy)
        vals=(t,sim.model.clock.iteration,sim.Δt,time()-started,mass,energy,urms,co,re,mach,.5mean(ub .^ 2 .+ vb .^ 2),skew...,qtop,cooling,minimum(a.rho),minimum(a.T),maximum(abs,a.wf[:,:,[1,end]]),numean,power,residual,centered_energy,centered_residual,ke,budget.heat,cumulative_defect,cumulative_defect/max(budget.initial_ke,1e-30),omega_now,p.f*omega_now/p.omega)
        all(isfinite,vals)||error("Nonfinite diagnostic")
        minimum(a.rho)>0&&minimum(a.T)>0||error("Nonpositive state")
        abs(mass/initial_mass-1)<1e-7||error("Mass conservation failed")
        mach<3||error("Mach exceeds 3: stop and inspect")
        open(path,"a") do io;println(io,join(vals,','));end
        energy_budget_notice(residual,t)
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]||isapprox(t,cfg["restart_time"];atol=1e-9)
            dye=Array(interior(sim.model.tracers.dye)) ./ a.rho
            jldsave(joinpath(out,@sprintf("frame_%07d.jld2",round(Int,100t)));t,z,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),density_profile=vec(mean(a.rho;dims=(1,2))),temperature_profile=vec(mean(a.T;dims=(1,2))),w_rms_profile=vec(sqrt.(mean(a.w .^ 2;dims=(1,2)))),enthalpy_flux_profile=vec(mean(2.5 .* a.rho .* a.w .* a.T;dims=(1,2))),f=p.f,omega_actual=omega_now,f_actual=p.f*omega_now/p.omega)
        end
        if cfg["save_3d"] && (isapprox(t/cfg["field_interval"],round(t/cfg["field_interval"]);atol=1e-7)||t>=cfg["stop_time"])
            mkpath(joinpath(out,"fields3d"))
            jldsave(joinpath(out,"fields3d",@sprintf("state_%07d.jld2",round(Int,t)));t,z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),density=Float32.(a.rho),temperature=Float32.(a.T))
        end
        if sim.model.clock.iteration==0 || t>=cfg["stop_time"] || isapprox(t,round(t);atol=1e-7)
        @printf("t=%.3f step=%d dt=%.4g wall=%.1fs Re=%.3g Ro=%.3g Mach=%.3g masserr=%.3g\n",t,sim.model.clock.iteration,sim.Δt,time()-started,re,co,mach,mass/initial_mass-1);flush(stdout)
        end
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=min(saved.dt,cap)
        @assert isapprox(m.clock.time,cfg["restart_time"];atol=1e-9) "Unexpected checkpoint time"
        meta["restart_time_actual"]=Float64(m.clock.time)
        meta["restart_source"]=pickup
        meta["restart_density_contrast"]=maximum(interior(m.dynamics.dry_density))/minimum(interior(m.dynamics.dry_density))
        write_toml(joinpath(out,"metadata.toml"),meta)
        output(sim) # Persist actual restored state and initialize the cumulative budget at restart.
    end
    run!(sim;checkpoint_at_end=true)
    meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    meta["final_time"]=Float64(m.clock.time);meta["wall_elapsed"]=time()-started;meta["iterations"]=m.clock.iteration
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
