using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Breeze.CompressibleEquations: NoDivergenceDamping

# Chan--Mayr dimensionless units: R=Ttop=rhotop=ptop=Lz=1, g=12.
using Oceananigans.Coriolis: ConstantCartesianCoriolis
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
    (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],lx=4.,ly=4.,lz=1.,dx=4/cfg["nx"],dy=4/cfg["ny"],dz=1/cfg["nz"],
      omega=cfg["omega"],f=2cfg["omega"]*sind(cfg["latitude"]),fh=2cfg["omega"]*cosd(cfg["latitude"]),
      kc=.0004*cfg["bottom_flux"]/(GRAV/2.5),flux=cfg["bottom_flux"],tau=cfg["cooling_time"],heat=cfg["heat"],tt=1.)
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
@inline mass_density(i,j,k,F,p)=continuity_at(i,j,k,F,p.rho0)
@inline continuity_at(i,j,k,F,::Nothing)=F.ρᵈ[i,j,k]
@inline continuity_at(i,j,k,F,rho0)=rho0[i,j,k]
@inline heat_source(i,j,k,g,c,F,p)=p.heat*(conduction(i,j,k,F.T,p)-mass_density(i,j,k,F,p)*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
density_field(m)=Breeze.AtmosphereModels.continuity_density(m.dynamics)
function write_toml(path,d)
    open(path,"w") do io;TOML.print(io,d);end
end
function arrays(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .*(uf.+circshift(uf,(-1,0,0)));v=.5 .*(vf.+circshift(vf,(0,-1,0)));w=.5 .*(wf[:,:,1:end-1].+wf[:,:,2:end])
    rho=Array(interior(density_field(m)))
    size(rho,1)==1 && (rho=repeat(rho,size(u,1),size(u,2),1))
    T=Array(interior(m.temperature));theta=Array(interior(m.formulation.potential_temperature))
    (;uf,vf,wf,u,v,w,rho,T,theta)
end
function constraint(m,p)
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    d=(circshift(mx,(-1,0,0)).-mx)/p.dx+(circshift(my,(0,-1,0)).-my)/p.dy+diff(mz;dims=3)/p.dz
    scale=max(maximum(abs,mx)/p.dx,maximum(abs,my)/p.dy,maximum(abs,mz)/p.dz,1.)
    maximum(abs,d)/scale
end
function construct(grid,dyn,constants,p,cfg)
    rot=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=p.omega,latitude=cfg["latitude"]))
    AtmosphereModel(grid;dynamics=dyn,thermodynamic_constants=constants,coriolis=rot,closure=nothing,tracers=(:dye,),
        momentum_advection=WENO(order=9),
        scalar_advection=(ρθ=WENO(order=9),ρqᵛ=WENO(order=9,bounds=(0,1)),dye=WENO(order=9,bounds=(0,1))),
        forcing=(ρE=Forcing(heat_source;discrete_form=true,parameters=p),))
end
function build_model(cfg;arch=GPU())
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,4),y=(0,4),z=(0,1),halo=(6,6,6),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    dyn=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping=NoDivergenceDamping());standard_pressure=1.,base_pressure=reference(0.).p,reference_potential_temperature=z->reference(z).theta)
    cp=merge(p,(;rho0=nothing,anelastic=false))
    cm=construct(grid,dyn,constants,cp,cfg)
    rng=MersenneTwister(cfg["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
    theta=[reference(zcell(k,p)).theta*(1+cfg["perturbation"]*noise[i,j,k]*sinpi(clamp(zcell(k,p)/DC,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(cm;θ=theta,ρ=cm.dynamics.reference_state.density,qᵗ=0.,u=0.,v=0.,w=0.)
    dye=[.1+.05cospi(2*(i-.5)/p.nx)*cospi(2*(j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(cm.tracers.dye,dye.*Array(interior(density_field(cm))))
    info=Dict{String,Any}("initial_velocity_rms"=>0.,"noise_seed"=>cfg["seed"],"relative_theta_noise_coefficient"=>cfg["perturbation"],"description"=>"Hydrostatic reference plus horizontally zero-mean thermal noise; no seeded eddies")
    cfg["equations"]=="compressible" && return cm,merge(cp,(;seed_info=info))
    ref=ReferenceState(grid,constants;base_pressure=reference(0.).p,standard_pressure=1.,potential_temperature=z->reference(z).theta)
    copyto!(parent(ref.pressure),parent(cm.dynamics.reference_state.pressure));copyto!(parent(ref.density),parent(cm.dynamics.reference_state.density))
    parent(ref.temperature).=parent(ref.pressure)./parent(ref.density)
    ap=merge(p,(;rho0=ref.density,anelastic=true))
    am=construct(grid,AnelasticDynamics(ref;pressure_solver_float_type=Float64),constants,ap,cfg)
    shared_theta=Array(interior(cm.formulation.potential_temperature))
    set!(am;θ=shared_theta,qᵗ=0.,u=0.,v=0.,w=0.)
    set!(am.tracers.dye,dye.*Array(interior(density_field(am))))
    Oceananigans.TimeSteppers.update_state!(am)
    a=arrays(am);c=arrays(cm)
    info["theta_mapping_max_error"]=maximum(abs,a.theta.-c.theta)
    info["temperature_mapping_max_relative"]=maximum(abs,a.T./c.T.-1)
    info["temperature_mapping_mass_rms_relative"]=sqrt(sum(a.rho.*(a.T./c.T.-1).^2)/sum(a.rho))
    info["reference_density_max_error"]=maximum(abs,Array(interior(ref.density)).-Array(interior(cm.dynamics.reference_state.density)))
    info["reference_pressure_max_error"]=maximum(abs,Array(interior(ref.pressure)).-Array(interior(cm.dynamics.reference_state.pressure)))
    @assert info["theta_mapping_max_error"]<1e-12 && info["reference_density_max_error"]==0 && info["reference_pressure_max_error"]==0
    @assert info["temperature_mapping_max_relative"]<1e-4 && constraint(am,ap)<1e-10
    am,merge(ap,(;seed_info=info))
end
mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)
