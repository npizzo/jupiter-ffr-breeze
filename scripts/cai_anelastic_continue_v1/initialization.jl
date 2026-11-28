function build_anelastic(c;arch)
    # Obtain exactly the compressible discrete vertical reference on a cheap CPU grid.
    small=merge(c,Dict("nx"=>16,"ny"=>16,"perturbation"=>0.))
    donor,_=build_model(small;arch=CPU());r=donor.dynamics.reference_state
    rp=copy(parent(r.pressure));rr=copy(parent(r.density));donor=nothing;GC.gc()
    p=parameters(c)
    grid=RectilinearGrid(arch,Float64;size=(p.nx,p.ny,p.nz),x=(0,4),y=(0,4),z=(0,1),halo=(6,6,6),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=GRAV,energy_reference_temperature=0.)
    ref=ReferenceState(grid,constants;base_pressure=reference(0.).p,standard_pressure=1.,potential_temperature=z->reference(z).theta)
    copyto!(parent(ref.pressure),rp);copyto!(parent(ref.density),rr)
    parent(ref.temperature).=parent(ref.pressure)./parent(ref.density)
    m,p=construct(grid,AnelasticDynamics(ref;pressure_solver_float_type=Float64),constants,merge(p,(;rho0=ref.density,anelastic=true)),c)
    m,p
end

function build_science(c;arch)
    m,p=build_anelastic(c;arch)
    rng=MersenneTwister(c["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
    theta=[reference(zcell(k,p)).theta*(1+c["perturbation"]*noise[i,j,k]*sinpi(zcell(k,p))) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(m;θ=theta,qᵗ=0.,u=0.,v=0.,w=0.)
    dye=[.1+.05cospi(2*(i-.5)/p.nx)*cospi(2*(j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(m.tracers.dye,dye.*Array(interior(p.rho0)))
    if get(c,"smooth_velocity",false)
        set!(m;u=(x,y,z)->.02cospi(y/2),v=(x,y,z)->.01sinpi(x/2),w=0.)
    end
    Oceananigans.TimeSteppers.update_state!(m)
    info=Dict("initial_velocity_rms"=>0.,"noise_seed"=>c["seed"],"relative_theta_noise_coefficient"=>c["perturbation"],"description"=>"Same rest-start, mean-free thermal perturbation as Cai compressible configuration; fixed discrete density/pressure reference copied from compressible builder")
    m,merge(p,(;seed_info=info))
end
