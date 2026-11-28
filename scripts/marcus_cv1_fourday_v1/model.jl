using Breeze, Oceananigans, CUDA, JLD2, TOML, Statistics, LinearAlgebra, Printf
using Oceananigans: interior
using Oceananigans.Operators: ℑzᵃᵃᶠ
using Oceananigans.Coriolis: ConstantCartesianCoriolis
const FFT=Oceananigans.Solvers.FFTW
struct MomentumCoriolis{C} <: Oceananigans.Coriolis.AbstractRotation{Oceananigans.Coriolis.EnergyConserving}
    rotation::C
end
@inline aliases(M)=(u=M.ρu,v=M.ρv,w=M.ρw)
@inline Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,c.rotation,aliases(M))
@inline Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,c.rotation,aliases(M))
@inline Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c::MomentumCoriolis,M)=Oceananigans.Coriolis.z_f_cross_U(i,j,k,g,c.rotation,aliases(M))
write_toml(path,d)=open(io->TOML.print(io,d),path,"w")
# Dry native b = g δ/(1+δ); Bannon b = g δ. Match at cell centres,
# then use the SAME vertical interpolation as the native buoyancy force.
@inline function buoyancy_correction(i,j,k,grid,F,p)
    @inbounds δ=F.T[i,j,k]/p.Tref[i,j,k]-1
    @inbounds return p.g*p.rho[i,j,k]*δ^2/(1+δ)
end
@inline source_x(i,j,k,g,c,F,p)=(@inbounds -p.rate[i,j,k]*p.rho[i,j,k]*(F.u[i,j,k]-p.ubg[i,j,k]))
@inline source_y(i,j,k,g,c,F,p)=(@inbounds -.5*(p.rate[i,j,k]+p.rate[i,j-1,k])*p.rho[i,j,k]*F.v[i,j,k])
@inline function source_z(i,j,k,g,c,F,p)
    correction=ℑzᵃᵃᶠ(i,j,k,g,buoyancy_correction,F,p)
    @inbounds sponge=-.25*(p.rate[i,j,k]+p.rate[i,j,k-1])*(p.rho[i,j,k]+p.rho[i,j,k-1])*F.w[i,j,k]
    correction+sponge
end
@inline source_theta(i,j,k,g,c,F,p)=(@inbounds -p.rate[i,j,k]*(F.ρθ[i,j,k]-p.rho[i,j,k]*p.thbg[i,j,k]))
# This fixture exercises the benchmark machinery, not Jovian physics.
const REFCSV=joinpath(@__DIR__, "reference.csv")
function read_reference()
    rows=readlines(REFCSV)[2:end]
    reduce(vcat,permutedims.(parse.(Float64,split(row,',')) for row in rows))
end
function refinterp(z,col,table)
    @assert table[1,1]<=z<=table[end,1] "Reference interpolation outside documented range"
    j=clamp(searchsortedlast(view(table,:,1),z),1,size(table,1)-1)
    a=(z-table[j,1])/(table[j+1,1]-table[j,1])
    (1-a)*table[j,col]+a*table[j+1,col]
end
function fixture()
    Dict{String,Any}("nx"=>24,"ny"=>24,"nz"=>64,"lx"=>105e6,"ly"=>70e6,"zmin"=>-220e3,"zmax"=>155e3,
        "g"=>25.,"R"=>3600.,"cp"=>12600.,"T0"=>162.71551724137933,"p0"=>1e5,"rotation_rate"=>1.75e-4,"latitude"=>-23.,
        "shear_vorticity"=>9.1e-6,"epsilon"=>1.5,"Rv"=>4.6e6,"Lr"=>2e6,"z0"=>8200.,"Htop"=>50e3,"Hbot"=>120e3,
        "sponge_seconds"=>2e4,"dt"=>.5,"stop_time"=>2.,"fixture_only"=>true)
end
@inline smoothwindow(a)=a<=.3 ? 1. : a>=.5 ? 0. : (q=(a-.3)/.2;1-10q^3+15q^4-6q^5)
@inline background_u(y,z,c)=begin q=mod(y-z*cotd(c["latitude"])+c["ly"]/2,c["ly"])-c["ly"]/2;-c["shear_vorticity"]*q*smoothwindow(abs(q)/c["ly"]) end
@inline function reference_vorticity(x,y,c)
    e=c["epsilon"];η=hypot(x/e,y);ω=c["shear_vorticity"]*(e+1)/(e^2-e)
    η<=c["Rv"] && return ω
    W=e/(e+1)*ω*c["Rv"]^2/η^2*exp((c["Rv"]-η)/c["Lr"])
    W*(1+1/e)+W*(-2/η-1/c["Lr"])/η*(x^2/e^3+y^2)
end
function poisson(rhs,dx,dy)
    nx,ny=size(rhs);a=FFT.fft(rhs)
    for j=1:ny,i=1:nx
        λ=-4*(sinpi((i-1)/nx)^2/dx^2+sinpi((j-1)/ny)^2/dy^2)
        a[i,j]=λ<0 ? a[i,j]/λ : 0
    end
    real.(FFT.ifft(a))
end
function initialized_flow(c,z;vortex=true)
    nx=c["nx"];ny=c["ny"];nz=length(z);dx=c["lx"]/nx;dy=c["ly"]/ny;dz=z[2]-z[1]
    uf=zeros(nx,ny,nz);vf=similar(uf);fill!(vf,0);P=similar(uf);ubg=similar(uf);Pbg=similar(uf)
    f=2c["rotation_rate"]*sind(c["latitude"]);fh=2c["rotation_rate"]*cosd(c["latitude"])
    removed=zeros(nz)
    ddx(a)=(circshift(a,(-1,0))-circshift(a,(1,0)))/(2dx)
    ddy(a)=(circshift(a,(0,-1))-circshift(a,(0,1)))/(2dy)
    for k=1:nz
        H=z[k]>=c["z0"] ? c["Htop"] : c["Hbot"]
        Hzero=0>=c["z0"] ? c["Htop"] : c["Hbot"]
        scale=exp(-((z[k]-c["z0"])/H)^2+(c["z0"]/Hzero)^2)
        vort=[vortex ? reference_vorticity((-c["lx"]/2+(i-1)*dx)/scale,(-c["ly"]/2+(j-1)*dy-z[k]*cotd(c["latitude"]))/scale,c) : 0. for i=1:nx,j=1:ny]
        removed[k]=mean(vort);ψ=poisson(vort.-mean(vort),dx,dy)
        ubg[:,:,k]=[background_u(-c["ly"]/2+(j-.5)*dy,z[k],c) for i=1:nx,j=1:ny]
        uf[:,:,k]=-(circshift(ψ,(0,-1))-ψ)/dy+ubg[:,:,k]
        vf[:,:,k]=(circshift(ψ,(-1,0))-ψ)/dx
        for (U,V,target) in ((uf[:,:,k],vf[:,:,k],P),(ubg[:,:,k],zeros(nx,ny),Pbg))
            uc=.5*(U+circshift(U,(-1,0)));vc=.5*(V+circshift(V,(0,-1)))
            ζ=(V-circshift(V,(1,0)))/dx-(U-circshift(U,(0,1)))/dy
            ζc=.25*(ζ+circshift(ζ,(-1,0))+circshift(ζ,(0,-1))+circshift(ζ,(-1,-1)))
            target[:,:,k]=poisson(-ddx(uc.*ddx(uc)+vc.*ddy(uc))-ddy(uc.*ddx(vc)+vc.*ddy(vc))+f*ζc,dx,dy)
        end
    end
    function buoyancy(p,u)
        b=similar(p);b[:,:,2:end-1]=(p[:,:,3:end]-p[:,:,1:end-2])/(2dz)
        b[:,:,1]=(p[:,:,2]-p[:,:,1])/dz;b[:,:,end]=(p[:,:,end]-p[:,:,end-1])/dz
        b-fh*.5*(u+circshift(u,(-1,0,0)))
    end
    (;uf,vf,P,ubg,b=buoyancy(P,uf),bbg=buoyancy(Pbg,ubg),removed)
end
function build(c;arch=CPU(),vortex=true,rest=false)
    @assert c["reconstruction_authorized"] "Requires documented physical reconstruction"
    nx=c["nx"];ny=c["ny"];nz=c["nz"]
    grid=RectilinearGrid(arch,Float64;size=(nx,ny,nz),x=(-c["lx"]/2,c["lx"]/2),y=(-c["ly"]/2,c["ly"]/2),z=(c["zmin"],c["zmax"]),halo=(6,6,6),topology=(Periodic,Periodic,Bounded))
    constants=ThermodynamicConstants(Float64;molar_gas_constant=c["R"],dry_air_molar_mass=1.,dry_air_heat_capacity=c["cp"],gravitational_acceleration=c["g"],energy_reference_temperature=0.)
    table=read_reference()
    pref=refinterp(c["zmin"],2,table)
    θfun(z)=refinterp(z,5,table)
    ref=ReferenceState(grid,constants;base_pressure=pref,standard_pressure=c["p0"],potential_temperature=θfun)
    rho=vec(Array(interior(ref.density)));pr=vec(Array(interior(ref.pressure)));Tref=vec(Array(interior(ref.temperature)))
    θref=Tref.*(c["p0"]./pr).^(c["R"]/c["cp"])
    dz=(c["zmax"]-c["zmin"])/nz;z=[c["zmin"]+(k-.5)*dz for k=1:nz]
    flow=initialized_flow(c,z;vortex=vortex&&!rest)
    if rest;flow.uf.=0.;flow.vf.=0.;flow.P.=0.;flow.ubg.=0.;flow.b.=0.;flow.bbg.=0.;end
    θ=reshape(θref,1,1,nz).*(1 .+flow.b/c["g"]);θbg=reshape(θref,1,1,nz).*(1 .+flow.bbg/c["g"])
    @assert minimum(θ)>0
    rate=CenterField(grid);ubg=XFaceField(grid);thbg=CenterField(grid)
    set!(ubg,flow.ubg);set!(thbg,θbg)
    set!(rate,(x,y,z)->(max(0.,(abs(y)/c["ly"]-.4)/.1)^2+max(0.,1-min(z-c["zmin"],c["zmax"]-z)/c["sponge_depth"])^2)/c["sponge_seconds"])
    for a in (rate,ubg,thbg);Oceananigans.BoundaryConditions.fill_halo_regions!(a);end
    p=(;rho=ref.density,Tref=ref.temperature,rate,ubg,thbg,g=c["g"])
    force(fun)=Forcing(fun;discrete_form=true,parameters=p)
    rot=MomentumCoriolis(ConstantCartesianCoriolis(rotation_rate=c["rotation_rate"],latitude=c["latitude"]))
    m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref;pressure_solver_float_type=Float64),thermodynamic_constants=constants,coriolis=rot,closure=nothing,tracers=(:dye,),momentum_advection=WENO(order=9),scalar_advection=(ρθ=WENO(order=9),ρqᵛ=WENO(order=9,bounds=(0,1)),dye=WENO(order=9,bounds=(0,1))),forcing=(ρu=force(source_x),ρv=force(source_y),ρw=force(source_z),ρθ=force(source_theta)))
    set!(m;θ,qᵗ=0.,u=flow.uf,v=flow.vf,w=0.);set!(m.dynamics.pressure_anomaly,flow.P)
    dye=[.5+.4sinpi(16*(i-.5)/nx) for i=1:nx,j=1:ny,k=1:nz]
    set!(m.tracers.dye,dye.*reshape(rho,1,1,nz));Oceananigans.TimeSteppers.update_state!(m)
    m,p,(;z,rho,θref,removed=flow.removed,flow)
end
