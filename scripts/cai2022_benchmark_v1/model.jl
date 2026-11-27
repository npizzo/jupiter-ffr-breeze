using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Breeze.CompressibleEquations: NoDivergenceDamping
using Oceananigans.TurbulenceClosures: Smagorinsky
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ

# Cai 2022: H=R=Ttop=rhotop=1, eta=4, m=1, g=(m+1)eta=8.
using Oceananigans.Coriolis: ConstantCartesianCoriolis
struct MomentumCoriolis{C} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    rotation::C
end
@inline momentum_aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.rotation,momentum_aliases(M))
const CV=1.5; const RGAS=1.; const GRAV=8.; const DC=1.
function reference(z)
    t=1+4*(1-z); p=t^2
    (;T=t,rho=t,p,theta=t/p^.4)
end
function parameters(cfg)
    (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],lx=4.,ly=4.,lz=1.,dx=4/cfg["nx"],dy=4/cfg["ny"],dz=1/cfg["nz"],
      omega=cfg["omega"],f=2cfg["omega"],fh=0.,kc=cfg["bottom_flux"]/4,flux=cfg["bottom_flux"],
      heat=cfg["heat"],tt=1.,constant_heat=false,constant_heating=0.)
end
@inline zcell(k,p)=(k-.5)*p.dz
@inline cool_rate(k,p)=0.
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
@inline mass_density(i,j,k,F,p)=continuity_at(i,j,k,F,p.rho0)
@inline continuity_at(i,j,k,F,::Nothing)=F.ρᵈ[i,j,k]
@inline continuity_at(i,j,k,F,rho0)=rho0[i,j,k]
@inline heat_source(i,j,k,g,c,F,p)=p.heat*conduction(i,j,k,F.T,p)+stress_work(i,j,k,g,c,F,p)

@inline divu(i,j,k,F,p)=(@inbounds (F.u[i+1,j,k]-F.u[i,j,k])/p.dx+(F.v[i,j+1,k]-F.v[i,j,k])/p.dy+(F.w[i,j,k+1]-F.w[i,j,k])/p.dz)
@inline trace_stress(i,j,k,g,F,p)=(@inbounds -2/3*p.nu[i,j,k]*F.ρᵈ[i,j,k]*divu(i,j,k,F,p))
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
@inline sxx(i,j,k,F,p)=(@inbounds 2*F.ρᵈ[i,j,k]*p.nu[i,j,k]*((F.u[i+1,j,k]-F.u[i,j,k])/p.dx-divu(i,j,k,F,p)/3))
@inline syy(i,j,k,F,p)=(@inbounds 2*F.ρᵈ[i,j,k]*p.nu[i,j,k]*((F.v[i,j+1,k]-F.v[i,j,k])/p.dy-divu(i,j,k,F,p)/3))
@inline szz(i,j,k,F,p)=(@inbounds 2*F.ρᵈ[i,j,k]*p.nu[i,j,k]*((F.w[i,j,k+1]-F.w[i,j,k])/p.dz-divu(i,j,k,F,p)/3))
@inline avgxy(a,i,j,k)=(@inbounds (a[i,j,k]+a[i-1,j,k]+a[i,j-1,k]+a[i-1,j-1,k])/4)
@inline avgxz(a,i,j,k)=(@inbounds (a[i,j,k]+a[i-1,j,k]+a[i,j,k-1]+a[i-1,j,k-1])/4)
@inline avgyz(a,i,j,k)=(@inbounds (a[i,j,k]+a[i,j-1,k]+a[i,j,k-1]+a[i,j-1,k-1])/4)
@inline sxy(i,j,k,F,p)=(@inbounds avgxy(F.ρᵈ,i,j,k)*avgxy(p.nu,i,j,k)*((F.u[i,j,k]-F.u[i,j-1,k])/p.dy+(F.v[i,j,k]-F.v[i-1,j,k])/p.dx))
@inline sxz(i,j,k,F,p)= k==1 || k==p.nz+1 ? 0. : (@inbounds avgxz(F.ρᵈ,i,j,k)*avgxz(p.nu,i,j,k)*((F.u[i,j,k]-F.u[i,j,k-1])/p.dz+(F.w[i,j,k]-F.w[i-1,j,k])/p.dx))
@inline syz(i,j,k,F,p)= k==1 || k==p.nz+1 ? 0. : (@inbounds avgyz(F.ρᵈ,i,j,k)*avgyz(p.nu,i,j,k)*((F.v[i,j,k]-F.v[i,j,k-1])/p.dz+(F.w[i,j,k]-F.w[i,j-1,k])/p.dy))
# Flux of viscous stress work. No separate +rho epsilon in a TOTAL-energy equation.
@inline workx(i,j,k,F,p)=(@inbounds F.u[i,j,k]*(sxx(i-1,j,k,F,p)+sxx(i,j,k,F,p))/2+
    ((F.v[i-1,j,k]+F.v[i,j,k])*sxy(i,j,k,F,p)+(F.v[i-1,j+1,k]+F.v[i,j+1,k])*sxy(i,j+1,k,F,p))/4+
    ((F.w[i-1,j,k]+F.w[i,j,k])*sxz(i,j,k,F,p)+(F.w[i-1,j,k+1]+F.w[i,j,k+1])*sxz(i,j,k+1,F,p))/4)
@inline worky(i,j,k,F,p)=(@inbounds F.v[i,j,k]*(syy(i,j-1,k,F,p)+syy(i,j,k,F,p))/2+
    ((F.u[i,j-1,k]+F.u[i,j,k])*sxy(i,j,k,F,p)+(F.u[i+1,j-1,k]+F.u[i+1,j,k])*sxy(i+1,j,k,F,p))/4+
    ((F.w[i,j-1,k]+F.w[i,j,k])*syz(i,j,k,F,p)+(F.w[i,j-1,k+1]+F.w[i,j,k+1])*syz(i,j,k+1,F,p))/4)
@inline workz(i,j,k,F,p)= k==1 || k==p.nz+1 ? 0. : (@inbounds F.w[i,j,k]*(szz(i,j,k-1,F,p)+szz(i,j,k,F,p))/2+
    ((F.u[i,j,k-1]+F.u[i,j,k])*sxz(i,j,k,F,p)+(F.u[i+1,j,k-1]+F.u[i+1,j,k])*sxz(i+1,j,k,F,p))/4+
    ((F.v[i,j,k-1]+F.v[i,j,k])*syz(i,j,k,F,p)+(F.v[i,j+1,k-1]+F.v[i,j+1,k])*syz(i,j+1,k,F,p))/4)
@inline stress_work(i,j,k,g,c,F,p)=(workx(i+1,j,k,F,p)-workx(i,j,k,F,p))/p.dx+(worky(i,j+1,k,F,p)-worky(i,j,k,F,p))/p.dy+(workz(i,j,k+1,F,p)-workz(i,j,k,F,p))/p.dz

density_field(m)=Breeze.AtmosphereModels.continuity_density(m.dynamics)
function write_toml(path,d)
    open(path,"w") do io;TOML.print(io,d);end
end
function arrays(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .*(uf.+circshift(uf,(-1,0,0)));v=.5 .*(vf.+circshift(vf,(0,-1,0)));w=.5 .*(wf[:,:,1:end-1].+wf[:,:,2:end])
    rho=Array(interior(density_field(m)))
    size(rho,1)==1 && (rho=repeat(rho,size(u,1),size(u,2),1))
    T=Array(interior(m.temperature));theta=T./(rho.*T).^0.4
    (;uf,vf,wf,u,v,w,rho,T,theta)
end
function constraint(m,p)
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    d=(circshift(mx,(-1,0,0)).-mx)/p.dx+(circshift(my,(0,-1,0)).-my)/p.dy+diff(mz;dims=3)/p.dz
    scale=max(maximum(abs,mx)/p.dx,maximum(abs,my)/p.dy,maximum(abs,mz)/p.dz,1.)
    maximum(abs,d)/scale
end
function construct(grid,dyn,constants,p,cfg)
    rot=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=p.omega,latitude=90.))
    formulation=TotalEnergyFormulation(reference_temperature=0.)
    # Native Smagorinsky uses Delta=(dx*dy*dz)^(1/3). Match Cai Eq A7 exactly on this uniform grid.
    cs=cfg["deardorff_coefficient"]*sqrt(sqrt(p.dx*p.dy)*p.dz)/(p.dx*p.dy*p.dz)^(1/3)
    closure=Smagorinsky(coefficient=cs,Pr=Inf)
    nu=CenterField(grid); pp=merge(p,(;nu))
    force(fun)=Forcing(fun;discrete_form=true,parameters=pp)
    schemes=(ρE=WENO(order=cfg["advection_order"]),ρqᵛ=WENO(order=cfg["advection_order"],bounds=(0,1)),dye=WENO(order=cfg["advection_order"],bounds=(0,1)))
    model=AtmosphereModel(grid;dynamics=dyn,formulation,thermodynamic_constants=constants,coriolis=rot,closure,tracers=(:dye,),
        momentum_advection=WENO(order=cfg["advection_order"]),scalar_advection=schemes,
        forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    model.closure_fields=(νₑ=nu,)
    return model,pp
end
function build_model(cfg;arch=GPU())
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,4),y=(0,4),z=(0,1),halo=(6,6,6),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    dyn=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping=NoDivergenceDamping());standard_pressure=1.,base_pressure=reference(0.).p,reference_potential_temperature=z->reference(z).theta)
    cp=merge(p,(;rho0=nothing,anelastic=false))
    cm,cp=construct(grid,dyn,constants,cp,cfg)
    rng=MersenneTwister(cfg["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
    theta=[reference(zcell(k,p)).theta*(1+get(cfg,"perturbation",1e-5)*noise[i,j,k]*sinpi(clamp(zcell(k,p)/DC,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(cm;θ=theta,ρ=cm.dynamics.reference_state.density,qᵗ=0.,u=0.,v=0.,w=0.)
    dye=[.1+.05cospi(2*(i-.5)/p.nx)*cospi(2*(j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(cm.tracers.dye,dye.*Array(interior(density_field(cm))))
    info=Dict{String,Any}("initial_velocity_rms"=>0.,"noise_seed"=>cfg["seed"],"relative_theta_noise_coefficient"=>cfg["perturbation"],"description"=>"Hydrostatic reference plus horizontally zero-mean thermal noise; no seeded eddies")
    if get(cfg,"smooth_velocity",false)
        temp=Array(interior(cm.temperature))
        set!(cm;ρ=cm.dynamics.reference_state.density,T=temp,qᵗ=0.,
             u=(x,y,z)->.02cospi(y/2), v=(x,y,z)->.01sinpi(x/2),
             w=(x,y,z)->.005sinpi(z)*cospi(x/2))
    end
    Oceananigans.TimeSteppers.update_state!(cm)
    cm,merge(cp,(;seed_info=info))
end

mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)
