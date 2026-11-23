using CUDA
using Breeze, Oceananigans, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans: interior
using Oceananigans.TurbulenceClosures.Smagorinskys: ΣᵢⱼΣᵢⱼᶜᶜᶜ
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
using Breeze.CompressibleEquations: NoDivergenceDamping, ThermalDivergenceDamping, UpperSponge, CubicRamp

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
    z>1 && return (;T=1.,rho=exp(-GRAV*(z-1)),p=exp(-GRAV*(z-1)),theta=exp(.4GRAV*(z-1)))
    ti=1+GRAV*(1-DC)/25
    if z>=DC
        t=1+GRAV*(1-z)/25; p=t^25
    else
        t=ti+GRAV*(DC-z)/2.5; p=ti^25*(t/ti)^2.5
    end
    return (;T=t,rho=p/t,p,theta=t/p^.4)
end
function parameters(cfg)
    return (;nx=cfg["nx"],ny=cfg["ny"],nz=cfg["nz"],lx=cfg["lx"],ly=cfg["ly"],dx=cfg["lx"]/cfg["nx"],dy=cfg["ly"]/cfg["ny"],lz=get(cfg,"lz",1.),dz=get(cfg,"lz",1.)/cfg["nz"],
        omega=cfg["omega"],f=2cfg["omega"]*sind(cfg["latitude"]),pr=cfg["pr"],kc=get(cfg,"conductivity",.0004*.03125/(GRAV/2.5)),flux=cfg["bottom_flux"],
        bottom_factor=get(cfg,"bottom_factor",cfg["heat"]),conduction_factor=get(cfg,"conduction_factor",cfg["heat"]),top_factor=get(cfg,"top_factor",cfg["heat"]),cap_factor=get(cfg,"cap_factor",cfg["heat"]),cap_jet=get(cfg,"jet_design","deep")=="cap",cap_low=get(cfg,"cap_low",.82),cap_high=get(cfg,"cap_high",.96),
        horizontal_sponge=get(cfg,"horizontal_sponge",false),sponge_rate=get(cfg,"sponge_rate",8.),sponge_depth=get(cfg,"sponge_depth",.18),mean_restoration=get(cfg,"mean_restoration",false),neutral_coriolis=get(cfg,"neutral_coriolis",false),entropy_diffusion_factor=get(cfg,"entropy_diffusion_factor",1.),jet_speed=cfg["jet_speed"],jet_width=cfg["jet_width"],jet_tau=cfg["jet_tau"],fh=get(cfg,"horizontal_rotation_factor",1.)*2cfg["omega"]*cosd(cfg["latitude"]),beta=cfg["beta"],tau=cfg["cooling_time"],heat=cfg["heat"],tt=1.)
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
@inline jet_momentum(i,j,k,F,p)=.5*(mass_density(i-1,j,k,F,p)+mass_density(i,j,k,F,p))*jet_acceleration(j,k,p)
@inline function upper_rate(k,p)
    q=clamp((zcell(k,p)-(p.lz-p.sponge_depth))/p.sponge_depth,0.,1.)
    p.sponge_rate*q*q*(3-2q)
end
@inline sponge_x(i,j,k,F,p)=p.horizontal_sponge ? -upper_rate(k,p)*0.5*(F.ρᵈ[i-1,j,k]+F.ρᵈ[i,j,k])*(F.u[i,j,k]-jet_target((j-.5)*p.dy-p.ly/2,zcell(k,p),p)) : 0.
@inline sponge_y(i,j,k,F,p)=p.horizontal_sponge ? -upper_rate(k,p)*0.5*(F.ρᵈ[i,j-1,k]+F.ρᵈ[i,j,k])*F.v[i,j,k] : 0.
@inline mass_density(i,j,k,F,p)=continuity_at(i,j,k,F,p.rho0)
@inline continuity_at(i,j,k,F,::Nothing)=F.ρᵈ[i,j,k]
@inline continuity_at(i,j,k,F,rho0)=rho0[i,j,k]
@inline fx(i,j,k,g,c,F,p)=jet_momentum(i,j,k,F,p)+rotation_x(i,j,k,F,p)
@inline fy(i,j,k,g,c,F,p)=rotation_y(i,j,k,F,p)
@inline fz(i,j,k,g,c,F,p)=rotation_z(i,j,k,F,p)
@inline external_heating(i,j,k,F,p)=conduction(i,j,k,F.T,p)-p.cap_factor*mass_density(i,j,k,F,p)*(p.mean_restoration ? p.tmean[1,1,k]-p.tt : F.T[i,j,k]-p.tt)*cool_rate(k,p)
@inline heat_source(i,j,k,g,c,F,p)=external_heating(i,j,k,F,p)
function build_compressible(cfg;arch=GPU(),perturb=true)
    p=parameters(cfg)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,p.lx),y=(0,p.ly),z=(0,p.lz),halo=(6,6,6),topology=(Periodic,Bounded,Bounded))
    # Bind the forcing to the actual cached viscosity used by the native closure.
    nu=CenterField(grid);umean=Field{Nothing,Center,Center}(grid);tmean=Field{Nothing,Nothing,Center}(grid);p=merge(p,(;nu,umean,tmean,anelastic=false,rho0=nothing))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    th(z)=reference(z).theta
    damping=get(cfg,"damping",false) ? ThermalDivergenceDamping(coefficient=.1) : NoDivergenceDamping()
    sponge=get(cfg,"sponge",false) ? UpperSponge(damping_rate=get(cfg,"sponge_rate",8.),depth=get(cfg,"sponge_depth",.18),ramp=CubicRamp()) : nothing
    dynamics=CompressibleDynamics(SplitExplicitTimeDiscretization(;damping,sponge,substeps=get(cfg,"substeps",0)==0 ? nothing : cfg["substeps"]);standard_pressure=1.,base_pressure=reference(0.).p,reference_potential_temperature=th)
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    rot=ChannelCoriolis(BetaPlane(;f₀=p.f-p.ly*p.beta/2,β=p.beta,scheme=Oceananigans.Coriolis.EnergyConserving()),ConstantCartesianCoriolis(;fx=0.,fy=p.fh,fz=0.))
    model=AtmosphereModel(grid;dynamics,thermodynamic_constants=constants,coriolis=nothing,closure=nothing,tracers=(:dye,),momentum_advection=WENO(order=9),scalar_advection=(ρθ=WENO(order=9),ρqᵛ=WENO(order=9,bounds=(0,1)),dye=WENO(order=9,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    rng=MersenneTwister(cfg["seed"])
    noise=rand(rng,p.nx,p.ny,p.nz) .- .5;noise .-= mean(noise,dims=(1,2))
    theta=[th(zcell(k,p))*(1+(perturb ? cfg["perturbation"] : 0.)*noise[i,j,k]*sinpi(clamp(zcell(k,p)/DC,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model;θ=theta,ρ=model.dynamics.reference_state.density,qᵗ=0,u=0,v=0,w=0)
    seed_info=(cfg["initial_rms"]==0 && p.jet_speed==0) ? Dict{String,Any}("initial_rms"=>0.) : apply_seed!(model,p,cfg;perturb)
    p=merge(p,(;seed_info))
    rho=Array(interior(model.dynamics.dry_density))
    dye=[.1+.05cospi(2*(i-.5)*p.dx/p.lx)*cospi(2*(j-.5)*p.dy/p.ly) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(model.tracers.dye,dye .* rho)
    interior(p.tmean).=sum(interior(model.temperature);dims=(1,2))./(p.nx*p.ny)
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
    rho=Array(interior(density_field(m)));size(rho,1)==1 && (rho=repeat(rho,size(u,1),size(u,2),1));T=Array(interior(m.temperature))
    (;uf,vf,wf,u,v,w,rho,T)
end

@inline uc(i,j,k,F)=.5*(F.u[i,j,k]+F.u[i+1,j,k])
@inline vc(i,j,k,F)=.5*(F.v[i,j,k]+F.v[i,j+1,k])
@inline wc(i,j,k,F)=.5*(F.w[i,j,k]+F.w[i,j,k+1])
@inline rx_cell(i,j,k,F,p)=mass_density(i,j,k,F,p)*((p.f+p.beta*((j-.5)*p.dy-p.ly/2))*vc(i,j,k,F)-p.fh*wc(i,j,k,F))
@inline ry_cell(i,j,k,F,p)=-mass_density(i,j,k,F,p)*(p.f+p.beta*((j-.5)*p.dy-p.ly/2))*uc(i,j,k,F)
@inline rz_cell(i,j,k,F,p)=mass_density(i,j,k,F,p)*p.fh*uc(i,j,k,F)
@inline rotation_x(i,j,k,F,p)=p.neutral_coriolis ? .5*(rx_cell(i-1,j,k,F,p)+rx_cell(i,j,k,F,p)) : 0.
@inline rotation_y(i,j,k,F,p)=p.neutral_coriolis ? .5*(ry_cell(i,j-1,k,F,p)+ry_cell(i,j,k,F,p)) : 0.
@inline rotation_z(i,j,k,F,p)=p.neutral_coriolis ? .5*(rz_cell(i,j,k-1,F,p)+rz_cell(i,j,k,F,p)) : 0.

density_field(m)=Breeze.AtmosphereModels.continuity_density(m.dynamics)
pressure_field(m)=m.dynamics isa AnelasticDynamics ? m.dynamics.reference_state.pressure : m.dynamics.pressure
function constraint(m,p)
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    d=(circshift(mx,(-1,0,0)).-mx)/p.dx+diff(my;dims=2)/p.dy+diff(mz;dims=3)/p.dz
    scale=max(maximum(abs,mx)/p.dx,maximum(abs,my)/p.dy,maximum(abs,mz)/p.dz,1.)
    maximum(abs,d)/scale
end
function build_model(cfg;arch=GPU(),perturb=true)
    cm,cp=build_compressible(cfg;arch,perturb)
    cfg["equations"]=="compressible" && return cm,cp
    a=arrays(cm);theta=Array(interior(cm.formulation.potential_temperature));dye=Array(interior(cm.tracers.dye))./a.rho
    grid=cm.grid;constants=cm.thermodynamic_constants
    ref=ReferenceState(grid,constants;base_pressure=reference(0.).p,standard_pressure=1.,potential_temperature=z->reference(z).theta)
    copyto!(parent(ref.pressure),parent(cm.dynamics.reference_state.pressure));copyto!(parent(ref.density),parent(cm.dynamics.reference_state.density))
    parent(ref.temperature).=parent(ref.pressure)./parent(ref.density)
    nu=CenterField(grid);umean=Field{Nothing,Center,Center}(grid);tmean=Field{Nothing,Nothing,Center}(grid)
    p=merge(parameters(cfg),(;nu,umean,tmean,anelastic=true,rho0=ref.density))
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64),thermodynamic_constants=constants,coriolis=nothing,closure=nothing,tracers=(:dye,),momentum_advection=WENO(order=9),scalar_advection=(ρθ=WENO(order=9),ρqᵛ=WENO(order=9,bounds=(0,1)),dye=WENO(order=9,bounds=(0,1))),forcing=(ρu=force(fx),ρv=force(fy),ρw=force(fz),ρE=force(heat_source)))
    set!(m;θ=theta,qᵗ=0.,u=a.uf,v=a.vf,w=a.wf)
    b=arrays(m);set!(m.tracers.dye,dye.*b.rho)
    interior(p.tmean).=sum(interior(m.temperature);dims=(1,2))./(p.nx*p.ny)
    Oceananigans.TimeSteppers.update_state!(m)
    err=maximum((maximum(abs,a.uf.-b.uf),maximum(abs,a.vf.-b.vf),maximum(abs,a.wf.-b.wf)))
    @assert err<1e-12 && constraint(m,p)<1e-10
    info=Dict{String,Any}("velocity_mapping_error"=>err,"theta_mapping_error"=>maximum(abs,Array(interior(m.formulation.potential_temperature)).-theta),"initial_constraint"=>constraint(m,p),"temperature_mapping_rms"=>sqrt(mean((b.T./a.T.-1).^2)),"temperature_mapping_max"=>maximum(abs,b.T./a.T.-1),"mapping"=>"Same reference pressure/density, face velocities, specific theta and tracer; rho0 continuity and pref thermodynamics replace evolving compressible density/pressure. T cannot also match.")
    return m,merge(p,(;seed_info=info))
end
