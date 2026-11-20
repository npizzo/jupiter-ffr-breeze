using Breeze: AtmosphereModels
using Oceananigans.BoundaryConditions: fill_halo_regions!
function operator_checks(cfg)
    small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>40))
    m,p=build_model(small;arch=CPU(),perturb=false)
    @assert p.nu === m.closure_fields.νₑ
    @assert isapprox(p.kc*(GRAV/2.5)/p.flux,.0004;rtol=1e-12)
    @assert isapprox(p.f,2p.omega*sind(cfg["latitude"]);rtol=1e-12)
    @assert reference(1.).T==1. && reference(1.).rho==1.
    @assert abs(reference(DC-1e-8).p/reference(DC+1e-8).p-1)<1e-6
    @assert abs(reference(.1).theta/reference(.7).theta-1)<1e-12
    T=[reference(zcell(k,p)).T+.03cospi(2*(i-.5)/p.nx) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    qt=mean(2p.kc .* (T[:,:,end] .- p.tt) ./ p.dz)
    @assert abs(mean(Q)-(p.flux-qt))<1e-12
    # Conservative entropy diffusion must telescope even for nonuniform rho,T,nu.
    rho=[reference(zcell(k,p)).rho for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(p.nu,(x,y,z)->.001*(1+.2cospi(x/2)))
    F=(T=T,ρᵈ=rho)
    qs=[sgs_heating(i,j,k,F,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    @assert abs(mean(qs))<1e-11
    # Constant entropy gives exactly zero SGS heating, despite stratification.
    Tad=rho .^ (1/CV);Fa=(T=Tad,ρᵈ=rho)
    qa=[sgs_heating(i,j,k,Fa,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    @assert maximum(abs,qa)<1e-9
    # Stress trace correction and viscous heating for one-dimensional compression.
    set!(m;u=(x,y,z)->sinpi(x/2),v=0,w=0)
    Oceananigans.TimeSteppers.update_state!(m)
    set!(p.nu,.001);fill_halo_regions!(p.nu)
    F=Oceananigans.fields(m);err=0.
    for k=3:p.nz-2,j=1:p.ny,i=1:p.nx
        native=-AtmosphereModels.∂ⱼ_𝒯₁ⱼ(i,j,k,m.grid,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing)
        got=native+fx(i,j,k,m.grid,m.clock,F,p)
        want=4/3*.001*F.ρᵈ[i,j,k]*(-4sinpi(1/p.nx)^2/p.dx^2)*F.u[i,j,k]
        err=max(err,abs(got-want))
    end
    @assert err<1e-10
    # Adapter must retain the horizontal rotation component too.
    M=(ρu=ones(24,24,48),ρv=2ones(24,24,48),ρw=3ones(24,24,48))
    c=m.coriolis;g=m.grid;fy=2p.omega*cosd(cfg["latitude"]);fzv=p.f
    @assert isapprox(Oceananigans.Coriolis.x_f_cross_U(6,6,6,g,c,M),3fy-2fzv;atol=1e-12)
    @assert isapprox(Oceananigans.Coriolis.y_f_cross_U(6,6,6,g,c,M),fzv;atol=1e-12)
    @assert isapprox(Oceananigans.Coriolis.z_f_cross_U(6,6,6,g,c,M),-fy;atol=1e-12)
    println("PARAMETERS_FLUX_SGS_STRESS_ROTATION_PASS ",err);flush(stdout)
    write_toml(joinpath(ARGS[1],"PROFILE.toml"),Dict("density_contrast"=>reference(0.).rho,"pressure_scale_heights_lower"=>log(reference(0.).p/reference(DC).p),"pressure_scale_heights_upper"=>log(reference(DC).p),"pressure_scale_heights_total"=>log(reference(0.).p),"radiative_conductivity"=>p.kc,"note"=>"Printed polytropes/g/interface give ~6.2 scale heights TOTAL, not 6.2 in lower region as text states. Use explicit printed model parameters."))
    return err
end
