using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
using Breeze.CompressibleEquations: NoDivergenceDamping

# Units: d=g=cp=rho(z=0)=1, R=0.4, cv=0.6. No Jovian dimensional mapping.
const ZBOT=-0.85; const ZTOP=1.15; const RGAS=0.4; const CV=0.6
function reference(z,xi)
    tt=xi/RGAS
    t0=tt+ZTOP/(2RGAS)
    if z>=0
        t=tt+(ZTOP-z)/(2RGAS);rho=t/t0
    else
        t=t0-z/(4RGAS);rho=(t/t0)^3
    end
    p=rho*RGAS*t
    return (;T=t,rho,p,theta=t/p^RGAS)
end
function parameters(cfg)
    r=reference(.5,cfg["xi"])
    ds=.1/(RGAS*r.T)
    chi=sqrt(ds/(cfg["ra"]*cfg["pr"]))
    nu=cfg["pr"]*chi
    kc=chi*r.rho
    flux=kc/(2RGAS)
    return (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],dx=4/cfg["nx"],dy=4/cfg["ny"],dz=2/cfg["nz"],nu,chi,kc,flux,
            f=nu*sqrt(cfg["ta"]),tt=cfg["xi"]/RGAS,tau=cfg["cooling_time"],width=cfg["cooling_width"],heat=cfg["heat"])
end
@inline zcell(k,p)=ZBOT+(k-.5)*p.dz
@inline kcell(k,p)=ifelse(zcell(k,p)<0,2p.kc,p.kc)
# Effective series conductivity across the interface, including cut cells.
@inline function kface(k,p)
    zb=zcell(k-1,p);zt=zcell(k,p)
    lower=clamp(-zb,0.,p.dz)
    return p.dz/(lower/(2p.kc)+(p.dz-lower)/p.kc)
end
@inline function cool_rate(k,p)
    q=clamp((zcell(k,p)-(1-p.width))/(2p.width),0.,1.)
    return q*q*(3-2q)/p.tau
end
@inline function conduction(i,j,k,T,p)
    im=ifelse(i==1,p.nx,i-1);ip=ifelse(i==p.nx,1,i+1)
    jm=ifelse(j==1,p.ny,j-1);jp=ifelse(j==p.ny,1,j+1)
    km=max(1,k-1);kp=min(p.nz,k+1)
    @inbounds begin
        tc=T[i,j,k]
        horizontal=kcell(k,p)*((T[ip,j,k]-2tc+T[im,j,k])/p.dx^2+(T[i,jp,k]-2tc+T[i,jm,k])/p.dy^2)
        qb=ifelse(k==1,p.flux,-kface(k,p)*(tc-T[i,j,km])/p.dz)
        qt=ifelse(k==p.nz,-2p.kc*(p.tt-tc)/p.dz,-kface(k+1,p)*(T[i,j,kp]-tc)/p.dz)
    end
    return horizontal+(qb-qt)/p.dz
end
@inline function divergence(i,j,k,F,p)
    @inbounds return (F.u[i+1,j,k]-F.u[i,j,k])/p.dx+(F.v[i,j+1,k]-F.v[i,j,k])/p.dy+(F.w[i,j,k+1]-F.w[i,j,k])/p.dz
end
# Native ScalarDiffusivity has stress 2 rho nu sym(grad u). Supply only
# the missing -2/3 rho nu div(u) I, without overwriting package methods.
@inline trace_stress(i,j,k,g,F,p)=-2p.nu/3*F.ρᵈ[i,j,k]*divergence(i,j,k,F,p)
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
@inline function heat_source(i,j,k,g,c,F,p)
    div=divergence(i,j,k,F,p)
    strain2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3
    @inbounds return 2p.nu*F.ρᵈ[i,j,k]*strain2+p.heat*(conduction(i,j,k,F.T,p)-CV*F.ρᵈ[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
end
function build_model(cfg;arch=GPU(),perturb=true)
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,4),y=(0,4),z=(ZBOT,ZTOP),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=RGAS,dry_air_molar_mass=1.,dry_air_heat_capacity=1.,gravitational_acceleration=1.,energy_reference_temperature=0.)
    th(z)=reference(z,cfg["xi"]).theta
    dynamics=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping=NoDivergenceDamping());standard_pressure=1.,base_pressure=reference(0.,cfg["xi"]).p,reference_potential_temperature=th)
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=FPlane(f=p.f),closure=ScalarDiffusivity(ν=p.nu,κ=0.),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    rng=MersenneTwister(cfg["seed"])
    noise=rand(rng,p.nx,p.ny,p.nz) .- 0.5;noise .-= mean(noise,dims=(1,2))
    theta=[th(zcell(k,p))*(1+(perturb ? cfg["perturbation"] : 0.)*noise[i,j,k]*max(0.,sinpi(clamp(zcell(k,p),0.,1.)))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model;θ=theta,ρ=model.dynamics.reference_state.density,qᵗ=0,u=0,v=0,w=0)
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
function run_case(cfg,out;arch=GPU(),pickup=false)
    isdir(out)&&!isempty(readdir(out))&&error("Refuse overwrite $out")
    mkpath(out);m,p=build_model(cfg;arch)
    rho0=Array(interior(m.dynamics.dry_density));z=[zcell(k,p) for k=1:p.nz]
    # Explicit diffusion, cooling and rotation caps supplement the advective CFL.
    thermal_dt=.18minimum(rho0)*CV/(2p.kc*(p.dx^-2+p.dy^-2+p.dz^-2))
    viscous_dt=.18/(p.nu*(p.dx^-2+p.dy^-2+p.dz^-2))
    cap=min(cfg["max_dt"],thermal_dt,viscous_dt,.15p.tau,.15/max(abs(p.f),1e-20))
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    planes=[argmin(abs.(z .- zz)) for zz in (.2,.5,.8)];cz=findall(0 .< z .< 1)
    meta=merge(cfg,Dict("status"=>"RUNNING","f"=>p.f,"nu"=>p.nu,"chi_mid"=>p.chi,"conductivity_cz"=>p.kc,"bottom_flux"=>p.flux,"dt_cap"=>cap,"video_heights"=>z[planes],"initial_density_contrast"=>maximum(rho0)/minimum(rho0),"model_sha256"=>bytes2hex(sha256(read(@__FILE__))),"paper"=>"https://doi.org/10.1088/0004-637X/742/1/34","fidelity"=>"Reconstruction matching reported Ra, Pr, Ta, xi and layer geometry; cooling time/width chosen explicitly; published normalized flux and density contrast not assumed reproduced; WENO finite-volume differs from Pencil discretization","budget_status"=>"Energy balance monitored; no claim of closed discrete vorticity budget in this split-explicit compressible driver"))
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,iteration,dt,wall_seconds,mass,total_energy,urms_conv,Co,Re,max_mach,barotropic_ke,skew_lower,skew_mid,skew_upper,top_heat_flux,cooling_power,min_density,min_temperature,wall_w_max");end
    started=time();initial_mass=mean(rho0)*32
    function output(sim)
        a=arrays(sim.model);t=Float64(sim.model.clock.time)
        mass=mean(a.rho)*32;energy=mean(a.rho .* (CV .* a.T .+ reshape(z,1,1,:) .+ 0.5  .* (a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2)))*32
        urms=sqrt(3mean(a.w[:,:,cz] .^ 2));co=abs(p.f)/max(2pi*urms,1e-30);re=urms/(2pi*p.nu)
        mach=maximum(sqrt.(a.u .^ 2 .+ a.v .^ 2 .+ a.w .^ 2) ./ sqrt.((5/3)*RGAS .* a.T))
        zf=(a.vf .- circshift(a.vf,(1,0,0))) ./ p.dx .- (a.uf .- circshift(a.uf,(0,1,0))) ./ p.dy
        zeta=0.25  .* (zf .+ circshift(zf,(-1,0,0)) .+ circshift(zf,(0,-1,0)) .+ circshift(zf,(-1,-1,0)))
        skew=[mean(zeta[:,:,k] .^ 3)/max(mean(zeta[:,:,k] .^ 2)^1.5,1e-30) for k in planes]
        weight=a.rho ./ sum(a.rho;dims=3);ub=dropdims(sum(weight .* a.u;dims=3);dims=3);vb=dropdims(sum(weight .* a.v;dims=3);dims=3)
        qtop=mean(2p.kc .* (a.T[:,:,end] .- p.tt) ./ p.dz)
        cooling=mean(CV .* a.rho .* (a.T .- p.tt) .* reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))*32
        vals=(t,sim.model.clock.iteration,sim.Δt,time()-started,mass,energy,urms,co,re,mach,.5mean(ub .^ 2 .+ vb .^ 2),skew...,qtop,cooling,minimum(a.rho),minimum(a.T),maximum(abs,a.wf[:,:,[1,end]]))
        all(isfinite,vals)||error("Nonfinite diagnostic")
        minimum(a.rho)>0&&minimum(a.T)>0||error("Nonpositive state")
        abs(mass/initial_mass-1)<1e-7||error("Mass conservation failed")
        mach<.7||error("Mach exceeds 0.7: stop and inspect; do not continue automatically")
        open(path,"a") do io;println(io,join(vals,','));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            dye=Array(interior(sim.model.tracers.dye)) ./ a.rho
            jldsave(joinpath(out,@sprintf("frame_%07d.jld2",round(Int,100t)));t,z,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),density_profile=vec(mean(a.rho;dims=(1,2))),temperature_profile=vec(mean(a.T;dims=(1,2))),w_rms_profile=vec(sqrt.(mean(a.w .^ 2;dims=(1,2)))),f=p.f)
        end
        if cfg["save_3d"] && (isapprox(t/cfg["field_interval"],round(t/cfg["field_interval"]);atol=1e-7)||t>=cfg["stop_time"])
            mkpath(joinpath(out,"fields3d"))
            jldsave(joinpath(out,"fields3d",@sprintf("state_%07d.jld2",round(Int,t)));t,z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),density=Float32.(a.rho),temperature=Float32.(a.T))
        end
        @printf("t=%.3f step=%d dt=%.4g wall=%.1fs Re=%.3g Co=%.3g Mach=%.3g masserr=%.3g\n",t,sim.model.clock.iteration,sim.Δt,time()-started,re,co,mach,mass/initial_mass-1);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    if pickup!==false;set!(sim;checkpoint=pickup);sim.Δt=saved.dt;end
    run!(sim;checkpoint_at_end=true)
    meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    meta["final_time"]=Float64(m.clock.time);meta["wall_elapsed"]=time()-started;meta["iterations"]=m.clock.iteration
    write_toml(joinpath(out,"metadata.toml"),meta)
    return sim,meta
end
