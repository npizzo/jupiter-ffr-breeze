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
        omega=cfg["omega"],f=2cfg["omega"]*sind(cfg["latitude"]),pr=cfg["pr"],kc=get(cfg,"conductivity",.0004*.03125/(GRAV/2.5)),flux=cfg["bottom_flux"],
        bottom_factor=get(cfg,"bottom_factor",cfg["heat"]),conduction_factor=get(cfg,"conduction_factor",cfg["heat"]),top_factor=get(cfg,"top_factor",cfg["heat"]),cap_factor=get(cfg,"cap_factor",cfg["heat"]),cap_jet=get(cfg,"jet_design","deep")=="cap",cap_low=get(cfg,"cap_low",.82),cap_high=get(cfg,"cap_high",.96),
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
        horizontal=p.conduction_factor*p.kc*((T[ip,j,k]-2tc+T[im,j,k])/p.dx^2+(T[i,jp,k]-2tc+T[i,jm,k])/p.dy^2)
        qb=ifelse(k==1,p.bottom_factor*p.flux,-p.conduction_factor*p.kc*(tc-T[i,j,km])/p.dz)
        qt=ifelse(k==p.nz,-2p.top_factor*p.kc*(p.tt-tc)/p.dz,-p.conduction_factor*p.kc*(T[i,j,kp]-tc)/p.dz)
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
@inline function cap_weight(z,p)
    q=clamp((z-p.cap_low)/(p.cap_high-p.cap_low),0.,1.)
    q^3*(10-15q+6q^2)
end
@inline function cap_weight_z(z,p)
    q=clamp((z-p.cap_low)/(p.cap_high-p.cap_low),0.,1.)
    30q^2*(1-q)^2/(p.cap_high-p.cap_low)
end
@inline jet_target(y,z,p)=jet_profile(y,p)*(p.cap_jet ? cap_weight(z,p) : 1.)
@inline jet_acceleration(j,k,p)=p.jet_speed==0 ? 0. : (jet_target((j-.5)*p.dy-p.ly/2,zcell(k,p),p)-p.umean[1,j,k])/p.jet_tau
@inline jet_momentum(i,j,k,F,p)=.5*(F.ρᵈ[i-1,j,k]+F.ρᵈ[i,j,k])*jet_acceleration(j,k,p)
@inline fx(i,j,k,g,c,F,p)=∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)+jet_momentum(i,j,k,F,p)
@inline fy(i,j,k,g,c,F,p)=∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
@inline fz(i,j,k,g,c,F,p)=∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
@inline external_heating(i,j,k,F,p)=conduction(i,j,k,F.T,p)-p.cap_factor*F.ρᵈ[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p)
@inline function heat_source(i,j,k,g,c,F,p)
    div=divergence(i,j,k,F,p)
    strain2=ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3
    # Paper Eq. 3 is -rho*(T-Ttop)/tau, WITHOUT a cv factor.
    @inbounds return 2p.nu[i,j,k]*F.ρᵈ[i,j,k]*strain2+p.entropy_diffusion_factor*sgs_heating(i,j,k,F,p)+external_heating(i,j,k,F,p)
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
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=rot,closure=Smagorinsky(coefficient=cfg["smagorinsky_coefficient"],Pr=Inf),tracers=(:dye,),momentum_advection=(get(cfg,"momentum_scheme","weno5")=="centered2" ? Centered(order=2) : WENO(order=5)),scalar_advection=(ρθ=(get(cfg,"scalar_scheme","weno5")=="centered2" ? Centered(order=2) : WENO(order=5)),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
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
