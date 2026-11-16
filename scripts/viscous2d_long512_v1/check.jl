# The legacy audit is used only for this brief numerical-equivalence check.
include("../instability2d_long_v1/audit.jl")
# Avoid including the shared core module twice.
include_string(Main,replace(read(joinpath(@__DIR__,"model.jl"),String),"include(\"../instability2d_v1/core.jl\")"=>""),joinpath(@__DIR__,"model.jl"))
using Test
source,out=ARGS;ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];metrics=Dict{String,Any}()
function evolve(u,v,case;legacy=false,steps=12,t=0.,iteration=0)
 m,b=setup(u,v;case,t,iteration)
 if legacy;b=Audit(b.op;a=case=="forced" ? .0002 : 0.,nu=.0001,alpha=0.,dt=.0025);end
 sim=Simulation(m;Δt=.0025,stop_iteration=iteration+steps,verbose=false)
 add_callback!(sim,b,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
 legacy && add_callback!(sim,s->finish!(b,s.model),IterationInterval(1);name=:finish)
 started=time();run!(sim);elapsed=time()-started
 (;m,elapsed)
end
@testset "Same equations, consistent restart" begin
 for case in ("unforced","forced")
  ref=evolve(u0,v0,case;legacy=true);lean=evolve(u0,v0,case)
  ur,vr=state(ref.m);u,v=state(lean.m)
  error=sqrt(mean((u.-ur).^2 .+(v.-vr).^2))/sqrt(mean(ur.^2 .+vr.^2))
  @test error<1e-10
  half=evolve(u0,v0,case;steps=6);uh,vh=state(half.m)
  path=joinpath(out,case*"_restart.jld2");jldsave(path;u=uh,v=vh,t=Float64(half.m.clock.time),iteration=half.m.clock.iteration)
  dd=load(path);res=evolve(dd["u"],dd["v"],case;steps=6,t=dd["t"],iteration=dd["iteration"])
  us,vs=state(res.m);restart_error=sqrt(mean((u.-us).^2 .+(v.-vs).^2))/sqrt(mean(u.^2 .+v.^2))
  @test restart_error<1e-10
  @test isapprox(res.m.clock.time,lean.m.clock.time;atol=1e-14)
  bench=evolve(u0,v0,case;steps=100)
  metrics[case]=Dict("relative_reference_error"=>error,"relative_restart_error"=>restart_error,"seconds_per_step"=>bench.elapsed/100,"estimated_40_time_unit_seconds"=>bench.elapsed/100*16000)
 end
end
metrics["status"]="PASS";metrics["source_initial_sha256"]=bytes2hex(sha256(read(source)))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true)
