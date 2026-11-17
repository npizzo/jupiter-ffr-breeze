# Finite-volume/C-grid polynomial closure plus the existing linear band forcing.
include("../viscous2d_fast512_v1/model.jl")
const CASES=("cubic_control","biased_plus","biased_minus")
function coefficients(case)
    case in CASES || error("Unsupported case")
    (;a=.0002,b=case=="cubic_control" ? 0. : .002,c=.0002,f=case=="biased_minus" ? -1. : 1.)
end
mutable struct PolynomialTerms{PF,PB}
    op::BandOperator
    q::Matrix{ComplexF64}
    reaction::Matrix{ComplexF64}
    invlap::Matrix{Float64}
    linear::Matrix{Float64}
    forward::PF
    backward::PB
    b_f::Float64
    c::Float64
    dt::Float64
end
@noinline function polynomial_buffer!(r,q,bf,c,dt)
    maxrate=0.
    @inbounds for k in eachindex(q)
        z=real(q[k]);r[k]=complex(bf*z*z-c*z*z*z,0.)
        maxrate=max(maxrate,abs(2bf*z-3c*z*z))
    end
    isfinite(maxrate) && dt*maxrate<.2 || error("Reaction step guard: $(dt*maxrate)")
    nothing
end
@noinline function combine_sources!(q,r,lin,invlap)
    @inbounds for k in eachindex(q)
        q[k]=(lin[k]*q[k]+r[k])*invlap[k]
    end
    q[1,1]=0. # Exact zero Fourier mode; implements mean removal for the polynomial.
    nothing
end
function (b::PolynomialTerms)(m)
    o=b.op;nx,ny=o.nx,o.ny;idx=nx/o.lx;idy=ny/o.ly
    u=interior(m.velocities.u);v=interior(m.velocities.v)
    guard_cfl(u,v,nx,ny,idx,idy,b.dt)
    curl_to_buffer!(b.q,u,v,nx,ny,idx,idy)
    polynomial_buffer!(b.reaction,b.q,b.b_f,b.c,b.dt)
    mul!(b.q,b.forward,b.q);mul!(b.reaction,b.forward,b.reaction)
    combine_sources!(b.q,b.reaction,b.linear,b.invlap)
    mul!(b.q,b.backward,b.q)
    add_streamfunction_tendency!(interior(m.timestepper.Gⁿ.u),interior(m.timestepper.Gⁿ.v),b.q,nx,ny,idx,idy)
    nothing
end
function setup_polynomial(u,v;case,t=0.,dt=.0025,b_override=nothing,c_override=nothing)
    cf=coefficients(case);bf=isnothing(b_override) ? cf.b*cf.f : b_override*cf.f
    c=isnothing(c_override) ? cf.c : c_override
    o=BandOperator(size(u,1),size(u,2);kmin=12.,kmax=16.,allow_zonal=false)
    m=build_model(o;f=cf.f,nu=.0001);seed!(m,u,v);m.clock.time=t;m.clock.iteration=0
    q=zeros(ComplexF64,o.nx,o.ny);react=similar(q);fill!(react,0.)
    # FFTW.UNALIGNED permits the same plan on either separately allocated buffer.
    pf=FFTW.plan_fft!(q;flags=FFTW.MEASURE|FFTW.UNALIGNED)
    pb=FFTW.plan_bfft!(q;flags=FFTW.MEASURE|FFTW.UNALIGNED)
    invlap=map(x->x==0 ? 0. : -1/(x*o.nx*o.ny),o.lambda)
    lin=cf.a .*o.k2.*o.mask
    m,PolynomialTerms(o,q,react,invlap,lin,pf,pb,bf,c,dt)
end
@noinline function stable_step(u,v,nx,ny,idx,idy,bf,c,ceiling)
    umax=0.;vmax=0.;rate=0.
    @inbounds for j=1:ny,i=1:nx
        im=i==1 ? nx : i-1;jm=j==1 ? ny : j-1
        z=(v[i,j,1]-v[im,j,1])*idx-(u[i,j,1]-u[i,jm,1])*idy
        umax=max(umax,abs(u[i,j,1]));vmax=max(vmax,abs(v[i,j,1]))
        rate=max(rate,abs(2bf*z-3c*z*z))
    end
    dt=min(ceiling,.3/max(umax*idx+vmax*idy,eps()),.1/max(rate,eps()))
    isfinite(dt) && dt>=1e-5 || error("Step too small; preserve checkpoint and review")
    dt
end
function set_stable_step!(sim,b;ceiling=.0025)
    o=b.op;m=sim.model
    dt=stable_step(interior(m.velocities.u),interior(m.velocities.v),o.nx,o.ny,o.nx/o.lx,o.ny/o.ly,b.b_f,b.c,ceiling)
    sim.Δt=dt;b.dt=dt
    nothing
end
function source_hashes()
    paths=["model.jl","run.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(x=>bytes2hex(sha256(read(joinpath(@__DIR__,x)))) for x in paths)
end
