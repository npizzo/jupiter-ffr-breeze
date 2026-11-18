# Frozen polynomial model, with only native momentum advection disabled.
include("../feedback_maintenance_v1/model.jl")
function setup_noadvection(u,v;t=0.,dt=.0025)
    _,p=setup_polynomial(u,v;case="biased_plus",t,dt)
    m=build_model(p.op;f=1.,nu=.0001,advection=false)
    seed!(m,u,v);m.clock.time=t;m.clock.iteration=0
    m,MaintenanceTerms(p,zeros(ComplexF64,size(u)),BitMatrix((p.op.k2.>0).&(p.op.k2.<=4.)),false)
end
function formation_hashes()
    paths=["model.jl","run.jl","check.jl","../feedback_maintenance_v1/model.jl",
           "../viscous2d_feedback512_v1/model.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(x=>bytes2hex(sha256(read(joinpath(@__DIR__,x)))) for x in paths)
end
function potential_diagnostics(o,q)
    qh=FFTW.fft(q);n=length(q)^2
    g=real.(FFTW.ifft((.0001 .*o.lambda .-.0002 .*o.k2.*o.mask).*qh)) .-.002 .*q.^2 .+.0002 .*q.^3
    g.-=mean(g)
    potential=.5sum((.0001 .*o.lambda .-.0002 .*o.k2.*o.mask).*abs2.(qh))/n- .002mean(q.^3)/3 + .0002mean(q.^4)/4
    (;potential,gradient_norm_squared=mean(abs2,g))
end
