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
        tau=cfg["cooling_time"],heat=cfg["heat"],tt=1.)
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
# Dry specific entropy uses thermodynamic density p_ref/T, NOT continuity rho0.
# cp log(T)-log(p_ref) = cp log(theta) at standard pressure 1.
@inline entropy(i,j,k,F,p)=2.5log(F.T[i,j,k])-log(p.pref[i,j,k])
@inline conductance(i,j,k,F,p)=p.rho[i,j,k]*F.T[i,j,k]*p.nu[i,j,k]/p.pr
@inline function entropy_flux(i,j,k,ii,jj,kk,F,p,delta)
    @inbounds return -.5*(conductance(i,j,k,F,p)+conductance(ii,jj,kk,F,p))*(entropy(ii,jj,kk,F,p)-entropy(i,j,k,F,p))/delta
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
@inline trace_stress(i,j,k,g,F,p)=-2p.nu[i,j,k]/3*p.rho[i,j,k]*divergence(i,j,k,F,p)
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
@inline function heat_source(i,j,k,g,c,F,p)
    div=divergence(i,j,k,F,p)
    strain2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3
    # Paper Eq. 3 is -rho*(T-Ttop)/tau, WITHOUT a cv factor.
    @inbounds return 2p.nu[i,j,k]*p.rho[i,j,k]*strain2+sgs_heating(i,j,k,F,p)+p.heat*(conduction(i,j,k,F.T,p)-p.rho[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
end
function write_toml(path,d)
    open(path,"w") do io; TOML.print(io,d); end
end
function arrays(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .* (uf .+ circshift(uf,(-1,0,0)));v=.5 .* (vf .+ circshift(vf,(0,-1,0)));w=.5 .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
    rho=repeat(Array(interior(m.dynamics.reference_state.density)),size(u,1),size(u,2),1)
    T=Array(interior(m.temperature));theta=Array(interior(m.formulation.potential_temperature))
    (;uf,vf,wf,u,v,w,rho,T,theta)
end
function constraint(m,p)
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    d=(circshift(mx,(-1,0,0)).-mx)/p.dx+(circshift(my,(0,-1,0)).-my)/p.dy+(mz[:,:,2:end].-mz[:,:,1:end-1])/p.dz
    scale=max(maximum(abs,mx)/p.dx,maximum(abs,my)/p.dy,maximum(abs,mz)/p.dz,1.)
    (;absolute=maximum(abs,d),relative=maximum(abs,d)/scale)
end
function build_model(cfg;arch=GPU(),rest=false)
    # Construct the original compressible IC using its unchanged source.
    cc=rest ? merge(cfg,Dict("initial_rms"=>0.,"perturbation"=>0.)) : cfg
    cm,cp=CompressibleControl.build_model(cc;arch)
    ca=CompressibleControl.arrays(cm)
    theta=Array(interior(cm.formulation.potential_temperature))
    dye=Array(interior(cm.tracers.dye))./ca.rho
    grid=cm.grid;constants=cm.thermodynamic_constants
    ref=ReferenceState(grid,constants;base_pressure=reference(0.).p,standard_pressure=1.,potential_temperature=z->reference(z).theta)
    # The compressible reference uses a discrete Exner integration. Copy its
    # actual pressure/density including halos; independently reintegrating the
    # analytic profile with ReferenceState would not be a matched control.
    copyto!(parent(ref.pressure),parent(cm.dynamics.reference_state.pressure))
    copyto!(parent(ref.density),parent(cm.dynamics.reference_state.density))
    parent(ref.temperature) .= parent(ref.pressure) ./ parent(ref.density)
    nu=CenterField(grid);p=merge(parameters(cfg),(;nu,pref=ref.pressure,rho=ref.density))
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    rot=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=p.omega,latitude=cfg["latitude"]))
    m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64),thermodynamic_constants=constants,coriolis=rot,closure=Smagorinsky(coefficient=cfg["smagorinsky_coefficient"],Pr=Inf),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    m.closure_fields=(νₑ=nu,)
    set!(m;θ=theta,qᵗ=0.,u=ca.uf,v=ca.vf,w=ca.wf)
    a=arrays(m);set!(m.tracers.dye,dye.*a.rho)
    Oceananigans.TimeSteppers.update_state!(m)
    referr=maximum(maximum(abs,Array(interior(getproperty(ref,k))).-Array(interior(getproperty(cm.dynamics.reference_state,k)))) for k in (:density,:pressure))
    verr=maximum((maximum(abs,a.uf.-ca.uf),maximum(abs,a.vf.-ca.vf),maximum(abs,a.wf.-ca.wf)))
    terr=maximum(abs,a.theta.-theta)
    @assert referr<1e-10 && verr<1e-12 && terr<1e-12
    @assert constraint(m,p).relative<1e-12
    thref=Array(interior(ref.pressure)).^.6 ./ Array(interior(ref.density));Tref=Array(interior(ref.temperature))
    Pref=Array(interior(ref.pressure));P=Array(interior(cm.dynamics.pressure))
    audit=Dict("reference_max_abs_error"=>referr,"velocity_max_abs_error"=>verr,"theta_max_abs_error"=>terr,"tracer_max_abs_error"=>maximum(abs,Array(interior(m.tracers.dye))./a.rho.-dye),"anelastic_initial_constraint_relative"=>constraint(m,p).relative,"temperature_mapping_max_fraction"=>maximum(abs,a.T./ca.T.-1),"temperature_mapping_rms_fraction"=>sqrt(mean((a.T./ca.T.-1).^2)),"compressible_pressure_departure_max"=>maximum(abs,P./Pref.-1),"compressible_density_departure_max"=>maximum(abs,ca.rho./a.rho.-1),"theta_departure_max"=>maximum(abs,theta./thref.-1),"anelastic_thermodynamic_density_departure_max"=>maximum(abs,Tref./a.T.-1),"compressible_initialization"=>cp.seed_info,"mapping"=>"Identical face velocities, specific potential temperature/entropy and dye; identical hydrostatic reference. Continuity density frozen to rho0; thermodynamic pressure fixed to pref. Initial T and thermodynamic density cannot also match.","mass_weighting"=>"All stresses, flux conductances, cooling and prognostic weights use rho0. Entropy is cp*log(T)-log(pref), not cv*log(T)-log(rho0). This is the anelastic approximation of the original constitutive laws.")
    return m,merge(p,(;audit))
end
mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)
function run_case(cfg,out;arch=GPU(),pickup=false)
    ispath(out)&&error("Refuse overwrite $out");mkpath(out)
    m,p=build_model(cfg;arch);write_toml(joinpath(out,"INITIAL_MAPPING.toml"),p.audit)
    z=[zcell(k,p) for k=1:p.nz];planes=[argmin(abs.(z.-zz)) for zz in (.2,.5,.9)]
    cap=min(cfg["max_dt"],.1/(2abs(p.omega)))
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    function diffusion_limit(sim)
        nu_max=maximum(interior(sim.model.closure_fields.νₑ));rho_min=minimum(interior(m.dynamics.reference_state.density))
        diffusivity=max(4nu_max/3,(2.5/CV)*nu_max/p.pr,p.kc/(rho_min*CV))
        sim.Δt=min(sim.Δt,.08/(max(diffusivity,1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diffusion_limit,IterationInterval(1);name=:diffusion_limit)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    meta=merge(cfg,Dict("status"=>"RUNNING","equations"=>"Breeze AnelasticDynamics; dry ideal gas; full constant Cartesian rotation","time_unit"=>"box height / top isothermal sound speed, not Earth days","video_heights"=>z[planes],"f"=>p.f,"initial_mapping"=>p.audit,"long_target"=>1000.,"scope"=>"Bounded validation pilot; no sign-selection conclusion","energy_budget"=>"No compressible total-energy closure claim: anelastic pressure/buoyancy work differs; report constraint, theta and resolved kinetic energy."))
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,iteration,dt,wall_seconds,kinetic_energy,urms,max_mach,skew_lower,skew_mid,skew_upper,constraint_absolute,constraint_relative,theta_departure_max,thermodynamic_density_departure_max,min_temperature,wall_w_max,top_heat_flux,cooling_power,theta_mass,barotropic_ke,nu_mean");end
    started=time();theta_ref=Array(interior(m.dynamics.reference_state.pressure)).^.6 ./ Array(interior(m.dynamics.reference_state.density));Tref=Array(interior(m.dynamics.reference_state.temperature))
    function output(sim)
        a=arrays(sim.model);t=Float64(m.clock.time);cons=constraint(m,p)
        rz=.5 .* (a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
        ke=.5p.dx*p.dy*p.dz*(sum(a.rho.*(a.uf.^2 .+a.vf.^2))+sum(rz.*a.wf[:,:,2:end-1].^2))
        speed2=a.u.^2 .+a.v.^2 .+a.w.^2;urms=sqrt(mean(speed2));mach=maximum(sqrt.(speed2 ./((5/3).*a.T)))
        zf=(a.vf.-circshift(a.vf,(1,0,0)))./p.dx.-(a.uf.-circshift(a.uf,(0,1,0)))./p.dy
        zeta=.25 .* (zf.+circshift(zf,(-1,0,0)).+circshift(zf,(0,-1,0)).+circshift(zf,(-1,-1,0)))
        skew=[mean(zeta[:,:,k].^3)/max(mean(zeta[:,:,k].^2)^1.5,1e-30) for k in planes]
        weight=a.rho./sum(a.rho;dims=3);ub=dropdims(sum(weight.*a.u;dims=3);dims=3);vb=dropdims(sum(weight.*a.v;dims=3);dims=3)
        qtop=mean(2p.kc.*(a.T[:,:,end].-p.tt)./p.dz)
        cooling=16mean(a.rho.*(a.T.-p.tt).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
        vals=(t,m.clock.iteration,sim.Δt,time()-started,ke,urms,mach,skew...,cons.absolute,cons.relative,maximum(abs,a.theta./theta_ref.-1),maximum(abs,Tref./a.T.-1),minimum(a.T),maximum(abs,a.wf[:,:,[1,end]]),qtop,cooling,16mean(a.rho.*a.theta),.5mean(ub.^2 .+vb.^2),mean(interior(m.closure_fields.νₑ)))
        @assert all(isfinite,vals) && minimum(a.T)>0 && mach<1 "Nonfinite, nonpositive T or Mach >= 1: inspect approximation"
        @assert cons.relative<1e-8 && maximum(abs,a.wf[:,:,[1,end]])<1e-10 "Anelastic constraint or wall flux failed"
        open(path,"a") do io;println(io,join(vals,','));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            dye=Array(interior(m.tracers.dye))./a.rho
            jldsave(joinpath(out,@sprintf("frame_%07d.jld2",round(Int,100t)));t,z,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),density_profile=vec(mean(a.rho;dims=(1,2))),temperature_profile=vec(mean(a.T;dims=(1,2))),w_rms_profile=vec(sqrt.(mean(a.w.^2;dims=(1,2)))),enthalpy_flux_profile=vec(mean(2.5 .*a.rho.*a.w.*a.T;dims=(1,2))),f=p.f)
        end
        if cfg["save_3d"]&&(isapprox(t/cfg["field_interval"],round(t/cfg["field_interval"]);atol=1e-7)||t>=cfg["stop_time"])
            mkpath(joinpath(out,"fields3d"));jldsave(joinpath(out,"fields3d",@sprintf("state_%09d.jld2",round(Int,1000t)));t,z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),density=Float32.(a.rho),temperature=Float32.(a.T),theta=Float32.(a.theta))
        end
        @printf("t=%.3f step=%d dt=%.4g wall=%.1fs Mach=%.4g constraint=%.3g theta_depart=%.3g\n",t,m.clock.iteration,sim.Δt,time()-started,mach,cons.relative,maximum(abs,a.theta./theta_ref.-1));flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt
        @assert isapprox(m.clock.time,cfg["restart_time"];atol=1e-10)
    end
    run!(sim;checkpoint_at_end=true)
    meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    meta["final_time"]=Float64(m.clock.time);meta["wall_elapsed"]=time()-started
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
