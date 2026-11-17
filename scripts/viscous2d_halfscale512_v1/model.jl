# Reuse the frozen efficient implementation; only the coefficient differs.
include("../viscous2d_fast512_v1/model.jl")
function setup_halfscale(u,v;t=0.)
    o=BandOperator(512,512;kmin=24.,kmax=32.,allow_zonal=false)
    m=build_model(o;f=1.,nu=.0001)
    seed!(m,u,v);m.clock.time=t;m.clock.iteration=0
    m,FastTerms(o;a=.000175,dt=.00125)
end
function source_hashes()
    paths=["model.jl","run.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in paths)
end
