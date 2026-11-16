include("../instability2d_v1/core.jl")
using .Instability2D, Oceananigans, Statistics, LinearAlgebra, Printf, TOML, SHA, JLD2
using Oceananigans: TendencyCallsite

struct LinearTerms
    op::BandOperator
    multiplier::Matrix{Float64}
    dt::Float64
end
LinearTerms(o;a=0.,nu=.0001,dt=.0025)=LinearTerms(o,a.*o.k2.*o.mask.-nu.*o.lambda,dt)
function (b::LinearTerms)(m)
    o=b.op;u,v=state(m)
    cfl=b.dt*(maximum(abs,u)/(o.lx/o.nx)+maximum(abs,v)/(o.ly/o.ny))
    isfinite(cfl) && cfl<.5 || error("CFL/finite-state guard failed: $cfl")
    q=vorticity(o,u,v)
    # Identical linear instability and native Laplacian viscosity to Audit(alpha=0).
    fu,fv=velocity_from_curl(o,spectral_curl(o,q,b.multiplier))
    interior(m.timestepper.Gⁿ.u)[:,:,1].+=fu
    interior(m.timestepper.Gⁿ.v)[:,:,1].+=fv
    nothing
end
function setup(u,v;case,t=0.,iteration=0)
    case in ("forced","unforced") || error("Only forced/unforced authorized")
    o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    m=build_model(o;f=1.);seed!(m,u,v);m.clock.time=t;m.clock.iteration=0 # New model must initialize RK tendency caches; output keeps a separate global offset.
    b=LinearTerms(o;a=case=="forced" ? .0002 : 0.)
    m,b
end
