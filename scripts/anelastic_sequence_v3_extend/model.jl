using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, LinearAlgebra, SHA, Printf
using Oceananigans: interior
using Oceananigans.BoundaryConditions: fill_halo_regions!
@inline smooth(x)=clamp(x,0.,1.)^2*(3-2clamp(x,0.,1.))
@inline shape(z,p)=smooth((z-p.dc)/(p.height-p.dc))
@inline theta_reference(z,p)=exp(p.caplog*shape(z,p))
function params(c)
    (; nx=c["nx"],ny=c["ny"],nz=c["nz"],length=c["length"],height=c["height"],dc=c["convective_depth"],g=c["g"],cp=c["cp"],R=c["R"],f=c["f"],nu=c["nu"],kappa=c["kappa"],flux=c["bottom_theta_flux"],caplog=c["cap_log_theta"],tau=c["cap_relaxation_time"],dz=c["height"]/c["nz"],bl=c["bottom_layer_depth"])
end
@inline function cap_cooling(i,j,k,grid,clock,F,p)
    @inbounds theta=F.T[i,j,k]/p.Tref[i,j,k]*p.theta_ref[i,j,k]
    @inbounds return -p.rho[i,j,k]*shape((k-.5)*p.dz,p)*(theta-p.theta_ref[i,j,k])/p.tau
end
function build(c;arch=CPU(),thermal=true,rest=false,inviscid=false,initial=true)
    p=params(c)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,p.length),y=(0,p.length),z=(0,p.height),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=p.R,dry_air_molar_mass=1.,dry_air_heat_capacity=p.cp,gravitational_acceleration=p.g,energy_reference_temperature=0.)
    theta_fn=let basep=p; z->theta_reference(z,basep); end
    ref=ReferenceState(grid,constants;base_pressure=p.R,standard_pressure=p.R,potential_temperature=theta_fn)
    theta_ref=Field{Nothing,Nothing,Center}(grid);set!(theta_ref,theta_fn);fill_halo_regions!(theta_ref)
    p=merge(p,(rho=ref.density,Tref=ref.temperature,theta_ref=theta_ref))
    # Prescribed density-potential-temperature flux; dry potential temperature is prognostic.
    bc=FieldBoundaryConditions(bottom=FluxBoundaryCondition(thermal ? p.flux : 0.),top=FluxBoundaryCondition(0.))
    forcing=thermal ? (ρθ=Forcing(cap_cooling;discrete_form=true,parameters=p),) : NamedTuple()
    m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64),thermodynamic_constants=constants,coriolis=FPlane(f=p.f),closure=ScalarDiffusivity(ν=inviscid ? 0. : p.nu,κ=(ρθ=inviscid ? 0. : p.kappa,ρqᵛ=0.,dye=0.)),tracers=(:dye,),boundary_conditions=(ρθ=bc,),forcing,momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))))
    if initial
        z=[(k-.5)*p.dz for k=1:p.nz];rho=vec(Array(interior(ref.density)))
        offset=thermal ? p.flux*p.tau/(sum(rho .* shape.(z,Ref(p)))*p.dz) : 0.
        rng=MersenneTwister(c["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
        theta=[theta_reference(z[k],p)+(rest ? 0. : offset+p.flux*p.bl/p.kappa*exp(-z[k]/p.bl)+c["noise_amplitude"]*noise[i,j,k]*sinpi(clamp(z[k]/p.dc,0.,1.))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
        set!(m;θ=theta,qᵗ=0.,u=0.,v=0.,w=0.)
        set!(m.tracers.dye,[rho[k]*(.1+.05*cospi(2*(i-.5)/p.nx)*cospi(2*(j-.5)/p.ny)) for i=1:p.nx,j=1:p.ny,k=1:p.nz])
    end
    return m,p
end
function fields(m)
    (;u=Array(interior(m.velocities.u)),v=Array(interior(m.velocities.v)),w=Array(interior(m.velocities.w)),theta=Array(interior(m.formulation.potential_temperature)),T=Array(interior(m.temperature)),rho=Array(interior(m.dynamics.reference_state.density)),mx=Array(interior(m.momentum.ρu)),my=Array(interior(m.momentum.ρv)),mz=Array(interior(m.momentum.ρw)))
end
function invariants(m,p)
    a=fields(m);dx=p.length/p.nx;dy=p.length/p.ny
    div=(circshift(a.mx,(-1,0,0)).-a.mx)/dx+(circshift(a.my,(0,-1,0)).-a.my)/dy+(a.mz[:,:,2:end].-a.mz[:,:,1:end-1])/p.dz
    (; mass=sum(a.rho.*a.theta)*dx*dy*p.dz,div=maximum(abs,div),wall=maximum(abs,a.w[:,:,[1,end]]),speed=max(maximum(abs,a.u),maximum(abs,a.v),maximum(abs,a.w)),Tmin=minimum(a.T))
end
function evolve!(m,stop;dt=.01,callback=nothing,interval=1)
    sim=Simulation(m;Δt=dt,stop_time=stop,wall_time_limit=180.)
    if callback!==nothing;add_callback!(sim,callback,IterationInterval(interval);name=:audit);end
    run!(sim)
    @assert m.clock.time>=stop-1e-9 "CPU validation wall limit"
    sim
end
function record(path,d)
    open(path,"w") do io;TOML.print(io,d);end
end
