using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
using Breeze.CompressibleEquations: NoDivergenceDamping

# Chan--Mayr dimensionless units: R=Ttop=rhotop=ptop=Lz=1, g=12.
using Oceananigans.Coriolis: ConstantCartesianCoriolis
using Oceananigans.TurbulenceClosures: Smagorinsky
# Preserve the existing horizontal rotation component; vary only vertical f(y).
struct ChannelCoriolis{V,H} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    vertical::V
    horizontal::H
end
@inline momentum_aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::ChannelCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.vertical,momentum_aliases(M))+Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.horizontal,momentum_aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::ChannelCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.vertical,momentum_aliases(M))+Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.horizontal,momentum_aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::ChannelCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.vertical,momentum_aliases(M))+Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.horizontal,momentum_aliases(M))
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
    return (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],lx=cfg["lx"],ly=cfg["ly"],dx=cfg["lx"]/cfg["nx"],dy=cfg["ly"]/cfg["ny"],dz=1/cfg["nz"],
        omega=cfg["omega"],f=2cfg["omega"]*sind(cfg["latitude"]),pr=cfg["pr"],kc=.0004*cfg["bottom_flux"]/(GRAV/2.5),flux=cfg["bottom_flux"],
        entropy_diffusion_factor=get(cfg,"entropy_diffusion_factor",1.),jet_speed=cfg["jet_speed"],jet_width=cfg["jet_width"],jet_tau=cfg["jet_tau"],fh=get(cfg,"horizontal_rotation_factor",1.)*2cfg["omega"]*cosd(cfg["latitude"]),beta=cfg["beta"],tau=cfg["cooling_time"],heat=cfg["heat"],tt=1.)
end
@inline zcell(k,p)=(k-.5)*p.dz
@inline function cool_rate(k,p)
    q=clamp((zcell(k,p)-DC)/(1-DC),0.,1.)
    return q*q*(3-2q)/p.tau
end
@inline function conduction(i,j,k,T,p)
    im=ifelse(i==1,p.nx,i-1);ip=ifelse(i==p.nx,1,i+1)
    jm=max(1,j-1);jp=min(p.ny,j+1)
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
    jm=max(1,j-1);jp=min(p.ny,j+1)
    qx=(entropy_flux(im,j,k,i,j,k,F,p,p.dx)-entropy_flux(i,j,k,ip,j,k,F,p,p.dx))/p.dx
    qy=((j==1 ? 0. : entropy_flux(i,jm,k,i,j,k,F,p,p.dy))-(j==p.ny ? 0. : entropy_flux(i,j,k,i,jp,k,F,p,p.dy)))/p.dy
    qb=k==1 ? 0. : entropy_flux(i,j,k-1,i,j,k,F,p,p.dz)
    qt=k==p.nz ? 0. : entropy_flux(i,j,k,i,j,k+1,F,p,p.dz)
    return qx+qy+(qb-qt)/p.dz
end
@inline function divergence(i,j,k,F,p)
    @inbounds return (F.u[i+1,j,k]-F.u[i,j,k])/p.dx+(F.v[i,j+1,k]-F.v[i,j,k])/p.dy+(F.w[i,j,k+1]-F.w[i,j,k])/p.dz
end
@inline trace_stress(i,j,k,g,F,p)=-2p.nu[i,j,k]/3*F.ρᵈ[i,j,k]*divergence(i,j,k,F,p)
@inline function jet_profile(y,p)
    q=clamp((abs(y)-2.5)/1.0,0.,1.)
    window=1-(6q^5-15q^4+10q^3)
    p.jet_speed*(y/p.jet_width)*exp(.5*(1-(y/p.jet_width)^2))*window
end
@inline jet_acceleration(j,k,p)=p.jet_speed==0 ? 0. : (jet_profile((j-.5)*p.dy-p.ly/2,p)-p.umean[1,j,k])/p.jet_tau
@inline jet_momentum(i,j,k,F,p)=.5*(F.ρᵈ[i-1,j,k]+F.ρᵈ[i,j,k])*jet_acceleration(j,k,p)
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)+jet_momentum(i,j,k,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
@inline function heat_source(i,j,k,g,c,F,p)
    div=divergence(i,j,k,F,p)
    strain2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3
    # Paper Eq. 3 is -rho*(T-Ttop)/tau, WITHOUT a cv factor.
    @inbounds return 2p.nu[i,j,k]*F.ρᵈ[i,j,k]*strain2+p.entropy_diffusion_factor*sgs_heating(i,j,k,F,p)+p.heat*(conduction(i,j,k,F.T,p)-F.ρᵈ[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
end
function build_model(cfg;arch=GPU(),perturb=true)
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,p.lx),y=(0,p.ly),z=(0,1),halo=(4,4,4),topology=(Periodic,Bounded,Bounded))
    # Bind the forcing to the actual cached viscosity used by the native closure.
    nu=CenterField(grid);umean=Field{Nothing,Center,Center}(grid);p=merge(p,(;nu,umean))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    th(z)=reference(z).theta
    dynamics=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping=NoDivergenceDamping());standard_pressure=1.,base_pressure=reference(0.).p,reference_potential_temperature=th)
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    rot=ChannelCoriolis(BetaPlane(;f₀=p.f-p.ly*p.beta/2,β=p.beta,scheme=Oceananigans.Coriolis.EnergyConserving()),ConstantCartesianCoriolis(;fx=0.,fy=p.fh,fz=0.))
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=rot,closure=Smagorinsky(coefficient=cfg["smagorinsky_coefficient"],Pr=Inf),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    model.closure_fields=(νₑ=nu,)
    rng=MersenneTwister(cfg["seed"])
    noise=rand(rng,p.nx,p.ny,p.nz) .- .5;noise .-= mean(noise,dims=(1,2))
    theta=[th(zcell(k,p))*(1+(perturb ? cfg["perturbation"] : 0.)*noise[i,j,k]*sinpi(clamp(zcell(k,p)/DC,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model;θ=theta,ρ=model.dynamics.reference_state.density,qᵗ=0,u=0,v=0,w=0)
    seed_info=apply_seed!(model,p,cfg;perturb)
    p=merge(p,(;seed_info))
    rho=Array(interior(model.dynamics.dry_density))
    dye=[.1+.05cospi(2*(i-.5)*p.dx/p.lx)*cospi(2*(j-.5)*p.dy/p.ly) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
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
    u=0.5  .* (uf .+ circshift(uf,(-1,0,0)));v=0.5  .* (vf[:,1:end-1,:] .+ vf[:,2:end,:]);w=0.5  .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
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
    cap=min(cfg["max_dt"],thermal_dt,viscous_dt,.15p.tau,.1/max(2abs(p.omega),1e-20))
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    averaged_u=Field(Average(m.velocities.u,dims=1))
    function update_jet_mean(model)
        compute!(averaged_u)
        interior(p.umean).=interior(averaged_u)
        nothing
    end
    update_jet_mean(m)
    add_callback!(sim,update_jet_mean,IterationInterval(1);name=:jet_mean,callsite=Oceananigans.UpdateStateCallsite())
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    # SGS coefficient is time dependent; enforce thermal diffusion stability each step.
    function diffusion_limit(sim)
        numax=maximum(interior(sim.model.closure_fields.νₑ))
        maxdiff=max(4numax/3, (2.5/CV)*numax/p.pr)
        sim.Δt=min(sim.Δt,.08/(max(maxdiff,1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diffusion_limit,IterationInterval(1);name=:diffusion_limit)
    health=open(joinpath(out,"step_health.csv"),"w")
    println(health,"iteration,time,min_rho,min_rhotheta,min_temperature,min_pressure,max_speed,max_nu")
    function guard(model)
        stats=(minimum(interior(model.dynamics.dry_density)),minimum(interior(model.formulation.potential_temperature_density)),minimum(interior(model.temperature)),minimum(interior(model.dynamics.pressure)),max(maximum(abs,interior(model.velocities.u)),maximum(abs,interior(model.velocities.v)),maximum(abs,interior(model.velocities.w))),maximum(interior(model.closure_fields.νₑ)))
        println(health,join((model.clock.iteration,model.clock.time,stats...),','))
        if !all(isfinite,stats) || min(stats[1:4]...)<=0
            flush(health);save_bad(model,out,"Nonfinite or nonpositive state in state-update guard")
            error("First invalid thermodynamic or velocity state captured")
        end
    end
    add_callback!(sim,guard,IterationInterval(1);name=:health_guard,callsite=Oceananigans.UpdateStateCallsite())
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    planes=[argmin(abs.(z .- zz)) for zz in (.2,.5,.9)];cz=findall(0 .< z .< DC)
    meta=merge(cfg,Dict("status"=>"RUNNING","f"=>p.f,"conductivity"=>p.kc,"dt_cap"=>cap,"video_heights"=>z[planes],"initial_density_contrast"=>maximum(rho0)/minimum(rho0),"model_sha256"=>bytes2hex(sha256(read(@__FILE__))),"paper"=>"https://doi.org/10.1016/j.epsl.2013.03.046","fidelity"=>"Chan--Mayr physical reconstruction, reduced grid; Cs=.16 cubic-volume filter, entropy diffusion, smoothstep cooling explicitly chosen; potential-temperature WENO scheme is not paper conservative total-energy discretization","budget_status"=>"Exploratory finite-amplitude evolution: 10%-of-heating residual is diagnostic only; no closed energy budget or equilibrated sign-transition claim","time_unit"=>"box-height / top isothermal sound speed, not Earth days"))
    meta["physics_unchanged_from"]="rotation_channel_v1 except widened y=8, balanced bounded jets, and optional zonal-mean relaxation; original heating/cap/SGS retained"
    meta["beta_definition"]="Vertical f=f_center+beta*(y-Ly/2); horizontal Coriolis component retained constant. Jupiter beta/f geometry only; dimensionless thermodynamics are not a dimensional Jupiter atmosphere."
    meta["energy_stop_policy"]="Warn at 10% of heating input; do not stop on this residual alone. Retain finite-state, positive rho/T, mass, Mach and timestep safeguards."
    meta["initialization"]=p.seed_info
    meta["experiment_scope"]="Finite-amplitude mixed-sign rotation comparison; not spontaneous formation or thermal equilibrium."
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,iteration,dt,wall_seconds,mass,total_energy,urms_conv,Ro,Re,max_mach,barotropic_ke,skew_lower,skew_mid,skew_upper,top_heat_flux,cooling_power,min_density,min_temperature,wall_w_max,nu_mean,external_power,energy_residual_over_input,centered_energy,centered_residual_over_input,kinetic_energy,cumulative_external_heat_since_restart,cumulative_energy_defect_since_restart,energy_defect_over_restart_ke,jet_power,jet_profile_rms_error,nonzonal_velocity_rms");end
    started=time();initial_mass=mean(rho0)*(p.lx*p.ly)
    if !isempty(p.seed_info);write_toml(joinpath(out,"INITIALIZATION.toml"),p.seed_info);end
    tracker=Ref((initial_energy=NaN,initial_ke=NaN,heat=0.))
    previous=Ref((time=NaN,energy=NaN,power=NaN,centered=NaN))
    function output(sim)
        a=arrays(sim.model);t=Float64(sim.model.clock.time)
        mass=mean(a.rho)*(p.lx*p.ly);centered_energy=mean(a.rho .* (CV .* a.T .+ GRAV .* reshape(z,1,1,:) .+ 0.5  .* (a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2)))*(p.lx*p.ly)
        rhox=.5 .* (a.rho .+ circshift(a.rho,(1,0,0)))
        rhoy=yface_density(a.rho)
        rhoz=.5 .* (a.rho[:,:,1:end-1] .+ a.rho[:,:,2:end])
        ke=(sum(rhox .* a.uf.^2)+sum(rhoy .* a.vf.^2)+sum(rhoz .* a.wf[:,:,2:end-1].^2))*.5*p.dx*p.dy*p.dz
        energy=mean(a.rho .* (CV .* a.T .+ GRAV .* reshape(z,1,1,:)))*(p.lx*p.ly)+ke
        urms=sqrt(mean(a.u[:,:,cz] .^ 2 .+ a.v[:,:,cz] .^ 2 .+ a.w[:,:,cz] .^ 2));co=urms/max(2abs(p.omega)*DC,1e-30)
        numean=mean(Array(interior(sim.model.closure_fields.νₑ))[:,:,cz]);re=urms*DC/max(numean,1e-30)
        mach=maximum(sqrt.(a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2) ./ sqrt.((5/3)*RGAS .* a.T))
        zeta=channel_vorticity(a.uf,a.vf,p)
        skew=[mean(zeta[:,:,k] .^ 3)/max(mean(zeta[:,:,k] .^ 2)^1.5,1e-30) for k in planes]
        weight=a.rho ./ sum(a.rho;dims=3);ub=dropdims(sum(weight .* a.u;dims=3);dims=3);vb=dropdims(sum(weight .* a.v;dims=3);dims=3)
        qtop=mean(2p.kc .* (a.T[:,:,end] .- p.tt) ./ p.dz)
        cooling=mean(a.rho .* (a.T .- p.tt) .* reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))*(p.lx*p.ly)
        mean_u=mean(a.uf;dims=1)
        target=[jet_profile((j-.5)*p.dy-p.ly/2,p) for i=1:1,j=1:p.ny,k=1:p.nz]
        acc=p.jet_speed==0 ? 0 .* mean_u : (target.-mean_u)./p.jet_tau
        jetpower=sum(rhox.*a.uf.*acc)*p.dx*p.dy*p.dz
        jeterror=sqrt(mean((mean_u.-target).^2))
        nonzonal=sqrt(mean((a.uf.-mean_u).^2)+mean(a.v.^2)+mean(a.w.^2))
        power=p.heat*((p.lx*p.ly)*(p.flux-qtop)-cooling)+jetpower
        old=previous[]
        residual=isfinite(old.time)&&t>old.time ? ((energy-old.energy)/(t-old.time)-.5*(power+old.power))/max((p.lx*p.ly)*p.flux,1e-30) : 0.
        centered_residual=isfinite(old.time)&&t>old.time ? ((centered_energy-old.centered)/(t-old.time)-.5*(power+old.power))/max((p.lx*p.ly)*p.flux,1e-30) : 0.
        budget=tracker[]
        if !isfinite(budget.initial_energy)
            budget=(initial_energy=energy,initial_ke=ke,heat=0.)
        elseif isfinite(old.time) && t>old.time
            budget=(initial_energy=budget.initial_energy,initial_ke=budget.initial_ke,heat=budget.heat+.5*(power+old.power)*(t-old.time))
        end
        tracker[]=budget
        cumulative_defect=energy-budget.initial_energy-budget.heat
        previous[]=(time=t,energy=energy,power=power,centered=centered_energy)
        vals=(t,sim.model.clock.iteration,sim.Δt,time()-started,mass,energy,urms,co,re,mach,.5mean(ub .^ 2 .+ vb .^ 2),skew...,qtop,cooling,minimum(a.rho),minimum(a.T),max(maximum(abs,a.wf[:,:,[1,end]]),maximum(abs,a.vf[:,[1,end],:])),numean,power,residual,centered_energy,centered_residual,ke,budget.heat,cumulative_defect,cumulative_defect/max(budget.initial_ke,1e-30),jetpower,jeterror,nonzonal)
        all(isfinite,vals)||error("Nonfinite diagnostic")
        minimum(a.rho)>0&&minimum(a.T)>0||error("Nonpositive state")
        abs(mass/initial_mass-1)<1e-7||error("Mass conservation failed")
        mach<get(cfg,"mach_stop",.5)||error("Mach exceeds diagnostic acceptance limit: stop and inspect")
        open(path,"a") do io;println(io,join(vals,','));end
        energy_budget_notice(residual,t)
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            dye=Array(interior(sim.model.tracers.dye)) ./ a.rho
            jldsave(joinpath(out,@sprintf("frame_%07d.jld2",round(Int,100t)));t,z,lx=p.lx,ly=p.ly,zonal_u=Float32.(dropdims(mean_u;dims=1)),jet_target=Float32.(dropdims(target;dims=1)),video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),density_profile=vec(mean(a.rho;dims=(1,2))),temperature_profile=vec(mean(a.T;dims=(1,2))),w_rms_profile=vec(sqrt.(mean(a.w .^ 2;dims=(1,2)))),enthalpy_flux_profile=vec(mean(2.5 .* a.rho .* a.w .* a.T;dims=(1,2))),f=p.f)
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
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt;update_jet_mean(m)
        @assert isapprox(m.clock.time,cfg["restart_time"];atol=1e-9) "Unexpected checkpoint time"
        meta["restart_time_actual"]=Float64(m.clock.time)
        meta["budget_initialization"]="First scheduled refreshed output after run initialization; no manual stale restart sample"
        write_toml(joinpath(out,"metadata.toml"),meta)
    end
    try
        run!(sim;checkpoint_at_end=true)
    catch err
        save_bad(m,out,sprint(showerror,err))
        rethrow()
    finally
        close(health)
    end
    meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    meta["final_time"]=Float64(m.clock.time);meta["wall_elapsed"]=time()-started;meta["iterations"]=m.clock.iteration
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
