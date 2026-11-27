using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, LinearAlgebra
using Oceananigans: interior
using Oceananigans.Coriolis: ConstantCartesianCoriolis
const FFT = Oceananigans.Solvers.FFTW

struct MomentumCoriolis{C} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    rotation::C
end
@inline aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.rotation,aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.rotation,aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.rotation,aliases(M))
write_toml(path,x)=open(io->TOML.print(io,x),path,"w")
density_field(m)=Breeze.AtmosphereModels.continuity_density(m.dynamics)

@inline capshape(z,c)=begin
    q=clamp((z-c["cap_base"]+.5c["transition"])/c["transition"],0.,1.)
    q*q*(3-2q)
end
function capintegral(z,c)
    a=c["cap_base"]-.5c["transition"];b=c["cap_base"]+.5c["transition"]
    z<=a && return 0.
    z>=b && return z-c["cap_base"]
    q=(z-a)/c["transition"]
    c["transition"]*(q^3-.5q^4)
end
function reference(grid,constants,c,N2)
    θb=c["bottom_temperature"]/(c["bottom_density"]*c["bottom_temperature"])^.4
    ReferenceState(grid,constants;base_pressure=c["bottom_density"]*c["bottom_temperature"],standard_pressure=1.,potential_temperature=z->θb*exp(N2*capintegral(z,c)/c["gravity"]))
end
# Rigid-lid hydrostatic displacement problem -(a_z/rho)_z=N²*a/(rho*c²).
# Dirichlet a=0 at lids; it remains regular when N²=0 in the lower region.
function internal_speed(rho,N2faces,dz)
    n=length(rho)-1
    A=SymTridiagonal((1 ./rho[1:end-1]+1 ./rho[2:end])/dz^2,-1 ./rho[2:end-1]/dz^2)
    rhof=0.5 .*(rho[1:end-1]+rho[2:end]);B=Diagonal(N2faces./rhof)
    L=cholesky(A).L;C=Symmetric(Matrix(L\B/L'))
    sqrt(maximum(eigvals(C)))
end
@inline function thermal_source(i,j,k,grid,clock,fields,p)
    z=(k-.5)*p.dz
    # Uniform horizontally, imposed enthalpy per volume per time.
    p.heating*(ifelse(z<p.thickness,p.flux/p.thickness,0.)-ifelse(z>p.lz-p.thickness,p.flux/p.thickness,0.))
end

# Resolved random phases on several scales. Odd nx+ny makes the entire field
# change sign under a half-box diagonal translation: exactly equal vorticity signs.
function horizontal_seed(c)
    nx=c["nx"];ny=c["ny"];dx=c["lx"]/nx
    modes=[(1,2),(2,1),(-1,2),(-2,1),(2,3),(3,2),(-2,3),(-3,2)]
    rng=MersenneTwister(c["seed"]);phases=2pi.*rand(rng,length(modes))
    psi=zeros(nx,ny)
    for ((kx,ky),ph) in zip(modes,phases),j=1:ny,i=1:nx
        psi[i,j]+=cos(2pi*(kx*(i-1)/nx+ky*(j-1)/ny)+ph)/(kx^2+ky^2)
    end
    u=-(circshift(psi,(0,-1))-psi)/dx;v=(circshift(psi,(-1,0))-psi)/dx
    ζ=(v-circshift(v,(1,0)))/dx-(u-circshift(u,(0,1)))/dx
    # Exact recorded discrete length convention; U is the maximum-plane face RMS.
    L=sqrt(mean(u.^2+v.^2)/mean(ζ.^2));U=c["rossby"]*c["f"]*L
    fac=U/sqrt(mean(u.^2+v.^2));psi.*=fac;u.*=fac;v.*=fac
    pc=0.25 .*(psi+circshift(psi,(-1,0))+circshift(psi,(0,-1))+circshift(psi,(-1,-1)))
    uc=0.5 .*(u+circshift(u,(-1,0)));vc=0.5 .*(v+circshift(v,(0,-1)))
    ddx(a)=(circshift(a,(-1,0))-circshift(a,(1,0)))/(2dx)
    ddy(a)=(circshift(a,(0,-1))-circshift(a,(0,1)))/(2dx)
    advx=uc.*ddx(uc)+vc.*ddy(uc);advy=uc.*ddx(vc)+vc.*ddy(vc)
    rhs=-(ddx(advx)+ddy(advy));rh=FFT.fft(rhs);nh=zeros(ComplexF64,nx,ny)
    for j=1:ny,i=1:nx
        λ=-4*(sinpi((i-1)/nx)^2+sinpi((j-1)/ny)^2)/dx^2
        if λ<0;nh[i,j]=rh[i,j]/λ;end
    end
    nl=real.(FFT.ifft(nh))
    (;psi,pc,u,v,nl,L,U,turnover=L/U,ν=U*L/c["reynolds"])
end

function fields_cpu(m)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=0.5 .*(uf+circshift(uf,(-1,0,0)));v=0.5 .*(vf+circshift(vf,(0,-1,0)));w=0.5 .*(wf[:,:,1:end-1]+wf[:,:,2:end])
    rho=Array(interior(density_field(m)));T=Array(interior(m.temperature));θ=Array(interior(m.formulation.potential_temperature))
    dye=Array(interior(m.tracers.dye))./rho
    (;uf,vf,wf,u,v,w,rho,T,θ,dye)
end
function constraint(m,p)
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    d=(circshift(mx,(-1,0,0))-mx)/p.dx+(circshift(my,(0,-1,0))-my)/p.dx+diff(mz;dims=3)/p.dz
    maximum(abs,d)/max(maximum(abs,mx)/p.dx,maximum(abs,my)/p.dx,maximum(abs,mz)/p.dz,1e-12)
end

function build(c,ratio;arch=CPU(),at_rest=false)
    nx=c["nx"];ny=c["ny"];nz=c["nz"];dz=c["lz"]/nz;dx=c["lx"]/nx
    grid=RectilinearGrid(arch,Float64;size=(nx,ny,nz),x=(0,c["lx"]),y=(0,c["lx"]),z=(0,c["lz"]),halo=(6,6,6),topology=(Periodic,Periodic,Bounded))
    tc=ThermodynamicConstants(Float64;molar_gas_constant=1.,dry_air_molar_mass=1.,dry_air_heat_capacity=2.5,gravitational_acceleration=c["gravity"],energy_reference_temperature=0.)
    seed=horizontal_seed(c)
    neutral=reference(grid,tc,c,0.);r0=vec(Array(interior(neutral.density)))
    shape=[capshape(k*dz,c) for k=1:nz-1]
    target_c=c["f"]*seed.L/ratio
    N2=(target_c/internal_speed(r0,shape,dz))^2
    for it=1:8
        ref=reference(grid,tc,c,N2)
        cr=internal_speed(vec(Array(interior(ref.density))),N2.*shape,dz)
        N2*=(target_c/cr)^2
    end
    ref=reference(grid,tc,c,N2);rho=vec(Array(interior(ref.density)))
    cr=internal_speed(rho,N2.*shape,dz);@assert abs(cr/target_c-1)<1e-6
    F=ones(nz);shape_error=0.
    closure=ScalarDiffusivity(ν=seed.ν,κ=seed.ν/c["prandtl"])
    rot=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=c["f"]/2,latitude=90.))
    Bflux=seed.U^3/c["cap_base"]
    flux=Bflux*c["bottom_density"]*2.5c["bottom_temperature"]/c["gravity"]
    hp=(;dz,flux,thickness=max(1,round(Int,c["heating_thickness"]/dz))*dz,lz=c["lz"],heating=at_rest ? 0. : c["heating"])
    m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64),thermodynamic_constants=tc,coriolis=rot,closure,tracers=(:dye,),momentum_advection=WENO(order=9),scalar_advection=(ρθ=WENO(order=9),ρqᵛ=WENO(order=9,bounds=(0,1)),dye=WENO(order=9,bounds=(0,1))),forcing=(ρE=Forcing(thermal_source;discrete_form=true,parameters=hp),))
    # This Breeze pin stores a scalar representative potential_temperature in
    # ReferenceState; recover the actual hydrostatic profile from p and rho.
    θ0=vec(Array(interior(ref.pressure))).^.6 ./rho
    # Nonlinear horizontal balance pressure; first-order hydrostatic buoyancy.
    P=c["f"].*reshape(seed.pc,nx,ny,1).*reshape(F,1,1,nz)+reshape(seed.nl,nx,ny,1).*reshape(F.^2,1,1,nz)
    b=zeros(nx,ny,nz)
    b[:,:,2:end-1].=(P[:,:,3:end]-P[:,:,1:end-2])/(2dz)
    b[:,:,1].=(P[:,:,2]-P[:,:,1])/(2dz);b[:,:,end].=(P[:,:,end]-P[:,:,end-1])/(2dz)
    at_rest && (b.=0.;P.=0.)
    theta=reshape(θ0,1,1,nz)./(1 .- b/c["gravity"])
    if !at_rest
        for k=1:nz,j=1:ny,i=1:nx
            z=(k-.5)*dz
            if z<c["cap_base"]-.5c["transition"]
                theta[i,j,k]*=1+c["thermal_seed"]*sinpi(2*6*(i-.5)/nx)*sinpi(2*5*(j-.5)/ny)*sinpi(z/(c["cap_base"]-.5c["transition"]))^2
            end
        end
    end
    minN2=minimum(c["gravity"].*diff(log.(theta);dims=3)/dz)
    capidx=findall(z->c["cap_base"]+c["transition"]<z<c["lz"]-c["heating_thickness"],[(k-.5)*dz for k=1:nz-1])
    capmin=minimum((c["gravity"].*diff(log.(theta);dims=3)/dz)[:,:,capidx]);@assert capmin>.9N2
    uf=reshape(seed.u,nx,ny,1).*reshape(F,1,1,nz);vf=reshape(seed.v,nx,ny,1).*reshape(F,1,1,nz)
    at_rest && (uf.=0.;vf.=0.)
    set!(m;θ=theta,qᵗ=0.,u=uf,v=vf,w=0.)
    set!(m.dynamics.pressure_anomaly,P)
    dye=[.1+.05cospi(2*(i-.5)/nx)*cospi(2*(j-.5)/ny) for i=1:nx,j=1:ny,k=1:nz]
    set!(m.tracers.dye,dye.*reshape(rho,1,1,nz))
    Oceananigans.TimeSteppers.update_state!(m)
    p=(;nx,ny,nz,dx,dz,z=[(k-.5)*dz for k=1:nz],rho,θ0,N2,F,seed,f=c["f"],g=c["gravity"],ratio,cr,flux,Bflux,capidx)
    @assert constraint(m,p)<1e-9 || at_rest
    a=fields_cpu(m);ζ=(a.vf-circshift(a.vf,(1,0,0)))/dx-(a.uf-circshift(a.uf,(0,1,0)))/dx
    skew=at_rest ? 0. : mean(ζ[:,:,end].^3)/mean(ζ[:,:,end].^2)^1.5
    @assert abs(skew)<1e-10
    hydroerr=maximum(abs,0.5 .*(b[:,:,1:end-1]+b[:,:,2:end])-diff(P;dims=3)/dz)/max(maximum(abs,b),1e-12)
    @assert hydroerr<.15
    info=Dict("requested_L_over_Ld"=>ratio,"measured_L_over_Ld"=>seed.L*c["f"]/cr,"internal_phase_speed"=>cr,"N2"=>N2,"minimum_initial_N2"=>minN2,"L"=>seed.L,"U"=>seed.U,"Ro"=>seed.U/(c["f"]*seed.L),"turnover_time"=>seed.turnover,"viscosity"=>seed.ν,"scalar_diffusivity"=>seed.ν/c["prandtl"],"enthalpy_flux"=>flux,"bottom_buoyancy_flux"=>Bflux,"cap_base"=>c["cap_base"],"minimum_initial_interior_cap_N2"=>capmin,"initial_vorticity_skewness"=>skew,"common_mode_relative_error"=>shape_error,"hydrostatic_interpolation_relative_error"=>hydroerr,"density_contrast"=>rho[1]/rho[end],"theta_departure_max"=>maximum(abs,theta./reshape(θ0,1,1,nz).-1),"maximum_initial_zeta_over_f"=>maximum(abs,ζ)/c["f"],"mode_shape"=>F,"reference_density"=>rho,"reference_theta"=>θ0)
    m,p,info
end

mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)
