include("model.jl")
out,source,target_str,validation=ARGS;target=parse(Float64,target_str);case="double_forced"
gate=TOML.parsefile(validation)
@assert gate["status"]=="PASS" && gate["forcing_coefficient"]==.0004
@assert gate["code_sha256"]==source_hashes()
ispath(out)&&error("Refusing existing output");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];t0=Float64(d["t"]);dt=.00125
@assert size(u0)==size(v0)==(512,512) && 160-1e-7<=t0<target<=400
it0=round(Int,t0/dt);itend=round(Int,target/dt)
@assert abs(it0*dt-t0)<1e-7 && abs(itend*dt-target)<1e-10
if haskey(d,"config")
    c=d["config"]
    for (k,val) in [("case",case),("viscosity",.0001),("drag",0.),("f",1.),("beta",0.),("instability_coefficient",.0004),("forcing_band",[12.,16.]),("allow_zonal_forcing",false)]
        @assert c[k]==val "Restart configuration mismatch: $k"
    end
end
initial_hash=haskey(d,"initial_sha256") ? d["initial_sha256"] : bytes2hex(sha256(read(source)))
@assert initial_hash==gate["initial_sha256"]
@assert haskey(d,"config")
if abs(t0-160)<1e-7
    @assert d["config"]["dt"]==.0025
    @assert bytes2hex(sha256(read(source)))==gate["restart_source_sha256"]
else
    @assert d["config"]["dt"]==dt && d["config"]["code_sha256"]==source_hashes()
end
m,b=setup_restart(u0,v0;t=t0,dt)
config=Dict("case"=>case,"grid"=>[512,512],"domain"=>[2pi,2pi],"f"=>1.,"beta"=>0.,"dt"=>dt,
 "viscosity"=>.0001,"drag"=>0.,"instability_coefficient"=>(.0004),
 "forcing_band"=>[12.,16.],"allow_zonal_forcing"=>false,"start_time"=>t0,"target_time"=>target,
 "final_campaign_time"=>400.,"initial_sha256"=>initial_hash,"source_sha256"=>bytes2hex(sha256(read(source))),
 "compute_location"=>"Unity","previous_dt"=>.0025,"restart_reason"=>"Original dt reached CFL guard near t=161; dt halved from saved t=160 checkpoint, unchanged physics","implementation"=>"Native ScalarDiffusivity; two planned in-place FFTs per forced stage, none for the unforced case","outputs"=>"Checkpoints and raw vorticity; downstream movie rendering on Unity; no scientific diagnostic analysis",
 "frame_interval"=>.5,"checkpoint_interval"=>10.,"internal_wall_seconds"=>5400.,
 "code_sha256"=>source_hashes())
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,file");end
started=time();lastframe=Ref(-1);lastcheckpoint=Ref(-1)
function checkpoint!()
    u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),iteration=it0+m.clock.iteration,config,initial_sha256=initial_hash)
    mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true)
end
function output!(sim;force=false)
    it=it0+m.clock.iteration
    if (it%400==0 || force) && lastframe[]!=it
        u,v=state(m);q=vorticity(b.op,u,v);all(isfinite,q)||error("Nonfinite saved field")
        name=@sprintf("zeta_%08d.f32",it)
        open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
        open(joinpath(out,"frames.csv"),"a") do io;println(io,join((it,m.clock.time,name),','));end
        lastframe[]=it
    end
    if (it%8000==0 || force) && lastcheckpoint[]!=it
        checkpoint!();lastcheckpoint[]=it
    end
    if it>it0 && it%2000==0
        println(@sprintf("PROGRESS case=%s t=%.3f segment_target=%.1f elapsed=%.1fs",case,m.clock.time,target,time()-started));flush(stdout)
    end
end
sim=Simulation(m;Δt=dt,stop_iteration=itend-it0,wall_time_limit=5400.,verbose=false)
add_callback!(sim,b,IterationInterval(1);name=:linear_terms,callsite=TendencyCallsite())
add_callback!(sim,output!,IterationInterval(1);name=:saved_fields)
output!(sim;force=true);run!(sim);output!(sim;force=true)
complete=m.clock.iteration>=itend-it0 && abs(m.clock.time-target)<1e-7
open(joinpath(out,"status.toml"),"w") do io
 TOML.print(io,Dict("status"=>(complete ? "COMPLETED" : "PARTIAL_WALL_TIME_LIMIT"),"t"=>m.clock.time,"target"=>target,"elapsed_seconds"=>time()-started);sorted=true)
end
complete||error("Partial output retained; downstream jobs stay stopped")
println("SEGMENT_COMPLETED case=$case t=$(m.clock.time)");flush(stdout)

