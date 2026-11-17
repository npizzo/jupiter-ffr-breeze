# Reuse the frozen efficient implementation; only the coefficient differs.
include("../viscous2d_fast512_v1/model.jl")
function setup_double(u,v;t=0.)
    o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    m=build_model(o;f=1.,nu=.0001)
    seed!(m,u,v);m.clock.time=t;m.clock.iteration=0
    m,FastTerms(o;a=.0004,dt=.0025)
end
function source_hashes()
    paths=["model.jl","run.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in paths)
end
