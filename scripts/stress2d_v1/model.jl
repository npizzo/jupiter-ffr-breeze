# Isolated strain-aligned stress hypothesis. Never modify previous campaign drivers.
include("../viscous2d_fast512_v1/model.jl")
const CASES=("decay","symmetric","biased_plus","biased_minus")
const GRID=256
const DTMAX=.005
const VISCOSITY=.0001
const KCUT=16.
function coefficients(case)
    case in CASES || error("Unknown stress case")
    off=case=="decay"
    (;nu0=off ? 0. : -.0002,nu1=case in ("biased_plus","biased_minus") ? -.001 : 0.,
      nu2=off ? 0. : .01,f=case=="biased_minus" ? -100. : 100.)
end
mutable struct StressTerms{PF,PB}
    op::BandOperator
    velocity::Matrix{ComplexF64}
    stress::Matrix{ComplexF64}
    zeta::Matrix{Float64}
    mask::Matrix{Float64}
    forward::PF
    backward::PB
    nu0::Float64
    nu1::Float64
    nu2::Float64
    f::Float64
    dt::Float64
    rate::Float64
    max_rossby::Float64
end
@noinline function pack_velocity!(a,u,v,n)
    @inbounds for j=1:n,i=1:n;a[i,j]=complex(u[i,j,1],v[i,j,1]);end
end
function project!(a,b)
    mul!(a,b.forward,a)
    @inbounds for k in eachindex(a);a[k]*=b.mask[k];end
    mul!(a,b.backward,a)
end
@noinline function stress_buffer!(r,z,a,n,ih,n0,n1,n2,f)
    maxro=0.;maxrate=0.
    @inbounds for j=1:n,i=1:n
        im=i==1 ? n : i-1;jm=j==1 ? n : j-1
        z[i,j]=((imag(a[i,j])-imag(a[im,j]))-(real(a[i,j])-real(a[i,jm])))*ih
    end
    @inbounds for j=1:n,i=1:n
        ip=i==n ? 1 : i+1;jp=j==n ? 1 : j+1
        im=i==1 ? n : i-1;jm=j==1 ? n : j-1
        sx=((real(a[ip,j])-real(a[i,j]))-(imag(a[i,jp])-imag(a[i,j])))*ih/2
        sy=((real(a[i,j])-real(a[i,jm]))+(imag(a[i,j])-imag(a[im,j])))*ih/2
        rc=z[i,j]/f;rn=(z[i,j]+z[ip,j]+z[i,jp]+z[ip,jp])/(4f)
        ntx=n0+n1*rn+n2*rn*rn;nty=n0+n1*rc+n2*rc*rc
        r[i,j]=complex(-2ntx*sx,-2nty*sy) # Rxx at centers, Rxy at corners; Ryy=-Rxx.
        ar=max(abs(rc),abs(rn));maxro=max(maxro,ar)
        bound=abs(n0)+abs(n1)*ar+abs(n2)*ar*ar+(abs(n1)+2abs(n2)*ar)*2(abs(sx)+abs(sy))/abs(f)
        maxrate=max(maxrate,bound*KCUT^2)
    end
    maxrate,maxro
end
@noinline function add_stress_divergence!(gu,gv,r,n,ih)
    @inbounds for j=1:n,i=1:n
        ip=i==n ? 1 : i+1;jp=j==n ? 1 : j+1
        im=i==1 ? n : i-1;jm=j==1 ? n : j-1
        gu[i,j,1]-=((real(r[i,j])-real(r[im,j]))+(imag(r[i,jp])-imag(r[i,j])))*ih
        gv[i,j,1]-=((imag(r[ip,j])-imag(r[i,j]))-(real(r[i,j])-real(r[i,jm])))*ih
    end
end
function (b::StressTerms)(m)
    n=b.op.nx;ih=n/(2pi);u=interior(m.velocities.u);v=interior(m.velocities.v)
    guard_cfl(u,v,n,n,ih,ih,b.dt)
    b.nu0==b.nu1==b.nu2==0 && return nothing
    pack_velocity!(b.velocity,u,v,n);project!(b.velocity,b)
    b.rate,b.max_rossby=stress_buffer!(b.stress,b.zeta,b.velocity,n,ih,b.nu0,b.nu1,b.nu2,b.f)
    isfinite(b.rate) && b.dt*b.rate<.15 || error("Stress timestep guard failed")
    b.max_rossby<.5 || error("Filtered |Z/f| reached 0.5; stop for review of expansion validity")
    project!(b.stress,b)
    add_stress_divergence!(interior(m.timestepper.Gⁿ.u),interior(m.timestepper.Gⁿ.v),b.stress,n,ih)
    nothing
end
function setup_stress(u,v;case,t=0.,dt=DTMAX,override=nothing,coriolis=nothing)
    c=isnothing(override) ? coefficients(case) : override;n=size(u,1)
    o=BandOperator(n,n;kmin=1.,kmax=KCUT,allow_zonal=true)
    m=build_model(o;f=isnothing(coriolis) ? c.f : coriolis,nu=VISCOSITY)
    seed!(m,u,v);m.clock.time=t
    a=zeros(ComplexF64,n,n);r=similar(a);z=zeros(n,n)
    pf=FFTW.plan_fft!(a;flags=FFTW.MEASURE|FFTW.UNALIGNED)
    pb=FFTW.plan_bfft!(a;flags=FFTW.MEASURE|FFTW.UNALIGNED)
    mask=Float64.(o.k2.<=KCUT^2)./(n*n)
    m,StressTerms(o,a,r,z,mask,pf,pb,c.nu0,c.nu1,c.nu2,c.f,dt,0.,0.)
end
@noinline function velocity_rate(u,v,n)
    a=0.;b=0.
    @inbounds for j=1:n,i=1:n;a=max(a,abs(u[i,j,1]));b=max(b,abs(v[i,j,1]));end
    (a+b)*n/(2pi)
end
function stable!(sim,b;ceiling=DTMAX)
    rate=velocity_rate(interior(sim.model.velocities.u),interior(sim.model.velocities.v),b.op.nx)
    dt=min(ceiling,.3/max(rate,eps()),.05/max(b.rate,eps()))
    dt>=.00025 && isfinite(dt) || error("Timestep below bounded-cost floor; stop for review")
    sim.Δt=dt;b.dt=dt
end
function source_hashes()
    paths=["model.jl","run.jl","reference.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in paths)
end
