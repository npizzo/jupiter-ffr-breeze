include("model.jl")
using Test
source,savedframe,out=ARGS;ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];t0=Float64(d["t"])
@assert abs(t0-160)<1e-7 && d["config"]["dt"]==.0025
for (k,v) in [("case","double_forced"),("instability_coefficient",.0004),("viscosity",.0001),("drag",0.),("f",1.),("beta",0.),("forcing_band",[12.,16.]),("allow_zonal_forcing",false)]
    @assert d["config"][k]==v
end
function evolve(u,v;steps=30,t=t0,dt=.00125)
    m,b=setup_restart(u,v;t,dt)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
    run!(sim);(;m,b,sim)
end
rel(a,b)=sqrt(sum(abs2,a.-b)/sum(abs2,b))
metrics=Dict{String,Any}()
@testset "t=160 restart with half timestep" begin
    o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    @test Float32.(vorticity(o,u0,v0))==read!(savedframe,Array{Float32}(undef,512,512))
    @test all(isfinite,u0) && all(isfinite,v0)
    new=evolve(u0,v0);un,vn=state(new.m)
    half=evolve(u0,v0;steps=13);u,v=state(half.m)
    p=joinpath(out,"restart.jld2");jldsave(p;u,v,t=Float64(half.m.clock.time));dd=load(p)
    res=evolve(dd["u"],dd["v"];steps=17,t=dd["t"]);ur,vr=state(res.m)
    restart_err=max(rel(un,ur),rel(vn,vr));@test restart_err<1e-10
    # Same physical interval at half of the new dt: short local timestep check.
    fine=evolve(u0,v0;steps=60,dt=.000625);uf,vf=state(fine.m)
    dt_err=max(rel(un,uf),rel(vn,vf));@test dt_err<1e-4
    new.sim.stop_iteration=130;GC.gc();perf=@timed run!(new.sim)
    @test perf.bytes/100<65536 && perf.time/100<.15
    new.sim.stop_iteration=1600;run!(new.sim)
    u,v=state(new.m);cfl=.00125*(maximum(abs,u)+maximum(abs,v))/(2pi/512)
    @test isfinite(cfl) && cfl<.5
    @test abs(new.m.clock.time-162)<1e-7
    merge!(metrics,Dict("checkpoint_frame_matches"=>true,"relative_restart_error"=>restart_err,
        "relative_short_timestep_error"=>dt_err,"seconds_per_step"=>perf.time/100,"bytes_per_step"=>perf.bytes/100,
        "pilot_end_time"=>new.m.clock.time,"pilot_end_cfl"=>cfl))
end
merge!(metrics,Dict("status"=>"PASS","forcing_coefficient"=>.0004,"dt"=>.00125,"previous_dt"=>.0025,
    "initial_sha256"=>d["initial_sha256"],"restart_source_sha256"=>bytes2hex(sha256(read(source))),
    "code_sha256"=>source_hashes(),"node"=>get(ENV,"SLURMD_NODENAME","")))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)
