# Native viscosity and a reusable, two-transform instability operator.
if !isdefined(Main, :Instability2D)
    include("../instability2d_v1/core.jl")
end
using .Instability2D, Oceananigans, Statistics, LinearAlgebra, Printf, TOML, SHA, JLD2
using .Instability2D: FFTW
using Oceananigans: TendencyCallsite

struct FastTerms{PF,PB}
    op::BandOperator
    buffer::Matrix{ComplexF64}
    multiplier::Matrix{Float64}
    forward::PF
    backward::PB
    active::Bool
    dt::Float64
end
function FastTerms(o; a=0., dt=.0025)
    buffer=zeros(ComplexF64,o.nx,o.ny)
    forward=FFTW.plan_fft!(buffer;flags=FFTW.MEASURE)
    backward=FFTW.plan_bfft!(buffer;flags=FFTW.MEASURE)
    # Combine the forcing multiplier, inverse discrete Laplacian and inverse-FFT normalization.
    multiplier=[o.lambda[i,j]==0 ? 0. : -a*o.k2[i,j]*o.mask[i,j]/(o.lambda[i,j]*o.nx*o.ny) for i=1:o.nx,j=1:o.ny]
    FastTerms(o,buffer,multiplier,forward,backward,a!=0.,dt)
end
# Function barriers specialize hot loops on the actual interior-view types.
# Without these barriers, inferred field-view unions caused scalar boxing.
@noinline function guard_cfl(u,v,nx,ny,idx,idy,dt)
    umax=0.;vmax=0.
    @inbounds for j=1:ny,i=1:nx
        umax=max(umax,abs(u[i,j,1]));vmax=max(vmax,abs(v[i,j,1]))
    end
    cfl=dt*(umax*idx+vmax*idy)
    isfinite(cfl) && cfl<.5 || error("CFL/finite-state guard failed: $cfl")
    nothing
end
@noinline function curl_to_buffer!(q,u,v,nx,ny,idx,idy)
    @inbounds for j=1:ny,i=1:nx
        im=i==1 ? nx : i-1;jm=j==1 ? ny : j-1
        q[i,j]=complex((v[i,j,1]-v[im,j,1])*idx-(u[i,j,1]-u[i,jm,1])*idy,0.)
    end
    nothing
end
@noinline function add_streamfunction_tendency!(gu,gv,q,nx,ny,idx,idy)
    @inbounds for j=1:ny,i=1:nx
        ip=i==nx ? 1 : i+1;jp=j==ny ? 1 : j+1
        gu[i,j,1]-=(real(q[i,jp])-real(q[i,j]))*idy
        gv[i,j,1]+=(real(q[ip,j])-real(q[i,j]))*idx
    end
    nothing
end
function (b::FastTerms)(m)
    o=b.op;nx,ny=o.nx,o.ny;idx=nx/o.lx;idy=ny/o.ly
    u=interior(m.velocities.u);v=interior(m.velocities.v)
    guard_cfl(u,v,nx,ny,idx,idy,b.dt)
    b.active || return nothing
    q=b.buffer
    curl_to_buffer!(q,u,v,nx,ny,idx,idy)
    mul!(q,b.forward,q)
    @inbounds for k in eachindex(q)
        q[k]*=b.multiplier[k]
    end
    mul!(q,b.backward,q)
    gu=interior(m.timestepper.Gⁿ.u);gv=interior(m.timestepper.Gⁿ.v)
    add_streamfunction_tendency!(gu,gv,q,nx,ny,idx,idy)
    nothing
end
function setup_fast(u,v;case,t=0.,iteration=0)
    case in ("forced","unforced") || error("Only the authorized pair is supported")
    o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    m=build_model(o;f=1.,nu=.0001)
    seed!(m,u,v);m.clock.time=t;m.clock.iteration=0
    b=FastTerms(o;a=case=="forced" ? .0002 : 0.)
    m,b
end
