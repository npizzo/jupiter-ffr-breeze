include("../viscous2d_long512_v1/model.jl")
include("model.jl")
using Test
root,out=ARGS;ispath(out)&&error("Output exists");mkpath(out)
function evolve(u,v,case;fast=true,steps=30,t=200.)
    m,b=fast ? setup_fast(u,v;case,t) : setup(u,v;case,t)
    sim=Simulation(m;Δt=.0025,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
    run!(sim)
    (;m,b,sim)
end
rel(a,b)=sqrt(sum(abs2,a.-b)/sum(abs2,b))
metrics=Dict{String,Any}()
@testset "Fast driver: same equations, restart, measured efficiency" begin
    for case in ("unforced","forced")
        source=joinpath(root,case,"to_200/checkpoint.jld2");d=load(source)
        @test abs(d["t"]-200)<1e-7
        ref=evolve(d["u"],d["v"],case;fast=false)
        fast=evolve(d["u"],d["v"],case)
        ur,vr=state(ref.m);uf,vf=state(fast.m)
        err=max(rel(uf,ur),rel(vf,vr));@test err<1e-10
        half=evolve(d["u"],d["v"],case;steps=13)
        u,v=state(half.m)
        # Exercise serialization, then reconstruct plans and RK state exactly as production does.
        p=joinpath(out,case*"_restart.jld2");jldsave(p;u,v,t=Float64(half.m.clock.time))
        dd=load(p);res=evolve(dd["u"],dd["v"],case;steps=17,t=dd["t"])
        us,vs=state(res.m);restart_err=max(rel(us,uf),rel(vs,vf));@test restart_err<1e-10
        ref.sim.stop_iteration=ref.m.clock.iteration+100;GC.gc();old=@timed run!(ref.sim)
        fast.sim.stop_iteration=fast.m.clock.iteration+100;GC.gc();new=@timed run!(fast.sim)
        speedup=old.time/new.time;@test speedup>=3
        @test new.bytes<old.bytes/10
        callback_bytes=@allocated fast.b(fast.m);@test callback_bytes<65536
        metrics[case]=Dict("relative_velocity_error_30_steps"=>err,"relative_restart_error"=>restart_err,
            "old_seconds_per_step"=>old.time/100,"fast_seconds_per_step"=>new.time/100,"speedup"=>speedup,
            "old_bytes_per_step"=>old.bytes/100,"fast_bytes_per_step"=>new.bytes/100,
            "fast_callback_bytes"=>callback_bytes,"estimated_t200_to400_compute_seconds"=>new.time/100*80000,
            "restart_source_sha256"=>bytes2hex(sha256(read(source))))
        println(case);TOML.print(stdout,metrics[case];sorted=true);flush(stdout)
    end
end
metrics["status"]="PASS";metrics["node"]=get(ENV,"SLURMD_NODENAME","")
metrics["dt"]=.0025;metrics["viscosity"]=.0001;metrics["forcing_coefficient"]=.0002
metrics["code_sha256"]=Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in [joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"../instability2d_v1/core.jl")])
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
