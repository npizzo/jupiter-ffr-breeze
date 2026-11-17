include("../viscous2d_long512_v1/model.jl")
include("model.jl")
using Test
source,baseline,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"]
initial_hash=bytes2hex(sha256(read(source)))
@assert initial_hash=="35cb85c5c7f179ac97dce1658f6779cd88e1d39c40b6193292412d250a7a8c6d"
@assert d["t"]==0 && size(u0)==size(v0)==(512,512)
function evolve(u,v;reference=false,steps=30,t=0.)
    if reference
        o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
        m=build_model(o;f=1.);seed!(m,u,v);m.clock.time=t
        b=LinearTerms(o;a=.0004,nu=.0001,dt=.0025)
    else
        m,b=setup_double(u,v;t)
    end
    sim=Simulation(m;Δt=.0025,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
    run!(sim);(;m,b,sim)
end
rel(a,b)=sqrt(sum(abs2,a.-b)/sum(abs2,b))
metrics=Dict{String,Any}()
@testset "Matched IC, doubled forcing, restart, efficient startup" begin
    o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    saved=read!(baseline,Array{Float32}(undef,512,512))
    @test Float32.(vorticity(o,u0,v0))==saved
    b1=FastTerms(o;a=.0002);b2=FastTerms(o;a=.0004)
    @test b2.multiplier==2 .* b1.multiplier
    @test b2.op.mask==b1.op.mask
    ref=evolve(u0,v0;reference=true);fast=evolve(u0,v0)
    ur,vr=state(ref.m);uf,vf=state(fast.m)
    err=max(rel(uf,ur),rel(vf,vr));@test err<1e-10
    half=evolve(u0,v0;steps=13);u,v=state(half.m)
    p=joinpath(out,"restart.jld2");jldsave(p;u,v,t=Float64(half.m.clock.time))
    dd=load(p);res=evolve(dd["u"],dd["v"];steps=17,t=dd["t"])
    us,vs=state(res.m);restart_err=max(rel(us,uf),rel(vs,vf));@test restart_err<1e-10
    fast.sim.stop_iteration=130;GC.gc();perf=@timed run!(fast.sim)
    @test perf.bytes/100<65536
    @test perf.time/100<.15
    fast.sim.stop_iteration=800;run!(fast.sim)
    u,v=state(fast.m);cfl=.0025*(maximum(abs,u)+maximum(abs,v))/(2pi/512)
    @test isfinite(cfl) && cfl<.5
    @test abs(fast.m.clock.time-2)<1e-10
    merge!(metrics,Dict("relative_velocity_error_30_steps"=>err,"relative_restart_error"=>restart_err,
        "seconds_per_step"=>perf.time/100,"bytes_per_step"=>perf.bytes/100,"pilot_end_time"=>fast.m.clock.time,
        "pilot_end_cfl"=>cfl,"identical_initial_saved_vorticity"=>true,"forcing_multiplier_exactly_doubled"=>true))
end
merge!(metrics,Dict("status"=>"PASS","forcing_coefficient"=>.0004,"dt"=>.0025,"viscosity"=>.0001,
    "initial_sha256"=>initial_hash,"code_sha256"=>source_hashes(),"node"=>get(ENV,"SLURMD_NODENAME","")))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)
