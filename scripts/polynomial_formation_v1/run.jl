include("model.jl")
out,source,target_arg,validation,case=ARGS;target=parse(Float64,target_arg)
gate=TOML.parsefile(validation);@assert gate["status"]=="PASS" && gate["code_sha256"]==formation_hashes()
ispath(out)&&error("Refusing existing output");mkpath(out)
d=load(source);u,v=d["u"],d["v"];t0=Float64(d["t"])
@assert 0.0 <= t0 < target <= gate["target_limit"] && size(u)==size(v)==(512,512)
if t0==0.
    @assert !haskey(d,"config") && case=="no_advection";origin_hash=bytes2hex(sha256(read(source)))
    @assert origin_hash==gate["initial_sha256"]
else
    cfg=d["config"];@assert cfg["case"]==case && cfg["code_sha256"]==formation_hashes()
    origin_hash=cfg["initial_sha256"];@assert origin_hash==gate["initial_sha256"]
end
m,b=setup_noadvection(u,v;t=t0)
config=Dict("case"=>case,"grid"=>[512,512],"domain"=>[2pi,2pi],"f"=>1.,"beta"=>0.,"a"=>.0002,"b"=>.002,"c"=>.0002,
 "viscosity"=>.0001,"forcing_band"=>[12.,16.],"allow_zonal_forcing"=>false,"velocity_drag"=>0.,"hyperviscosity"=>0.,
 "momentum_advection"=>false,"cubic_unchanged"=>true,"dt_ceiling"=>.0025,"start_time"=>t0,
 "target_time"=>target,"final_campaign_time"=>400.,"frame_interval"=>.5,"diagnostic_interval"=>.25,"checkpoint_interval"=>5.,
 "initial_sha256"=>origin_hash,"code_sha256"=>formation_hashes(),"compute_location"=>"Unity",
 "interpretation"=>"Matched formation intervention: disable only momentum advection",
 "budget_semantics"=>"Native instantaneous low-mode tendency work; trapezoidal cumulative samples are not exact every-stage budget integrals")
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,file");end
started=time();lastsample=Ref(-Inf);lastcheckpoint=Ref(-Inf);lastfields=Ref(-Inf)
qstart,startmetrics=maintenance_diagnostics(m,b);merge!(startmetrics,Dict(string(k)=>v for (k,v) in pairs(potential_diagnostics(b.base.op,qstart))));startenergy=startmetrics["energy_k_le_2"]
keys_sorted=sort(collect(keys(startmetrics)))
terms=("advection","linear","quadratic","cubic","viscous")
integrals=Dict(k=>0. for k in terms);previous=Ref(startmetrics);previoust=Ref(t0)
headers=vcat(["t","dt"],keys_sorted,["integral_"*k*"_lowk_work" for k in terms],["lowk_energy_change","sampled_lowk_budget_residual"])
open(joinpath(out,"scalars.csv"),"w") do io;println(io,join(headers,','));end
function checkpoint!(sim)
    abs(m.clock.time-lastcheckpoint[])<1e-8 && return
    u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),config)
    mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true);lastcheckpoint[]=m.clock.time
end
function output!(sim)
    abs(m.clock.time-lastsample[])<1e-8 && return
    q,di=maintenance_diagnostics(m,b);merge!(di,Dict(string(k)=>v for (k,v) in pairs(potential_diagnostics(b.base.op,q))));@assert abs(di["advection_lowk_energy_rate"])<1e-14;all(isfinite,q)&&all(isfinite,values(di))||error("Nonfinite state")
    t=Float64(m.clock.time);h=t-previoust[]
    for k in terms;integrals[k]+=.5h*(previous[][k*"_lowk_energy_rate"]+di[k*"_lowk_energy_rate"]);end
    delta=di["energy_k_le_2"]-startenergy;residual=delta-sum(values(integrals))
    vals=vcat([t,sim.Δt],[di[k] for k in keys_sorted],[integrals[k] for k in terms],[delta,residual])
    open(joinpath(out,"scalars.csv"),"a") do io;println(io,join(vals,','));end
    ix=round(Int,2t)
    if abs(ix*.5-t)<1e-7 && abs(t-lastfields[])>1e-8
        name=@sprintf("zeta_%08d.f32",ix);open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
        open(joinpath(out,"frames.csv"),"a") do io;println(io,join((m.clock.iteration,t,name),','));end
        lastfields[]=t
    end
    previous[]=di;previoust[]=t;lastsample[]=t
end
function progress!(sim)
 println(@sprintf("PROGRESS case=%s t=%.3f target=%.1f dt=%.6g elapsed=%.1fs",case,m.clock.time,target,sim.Δt,time()-started));flush(stdout)
end
sim=Simulation(m;Δt=.0025,stop_time=target,wall_time_limit=2050.,verbose=false)
maintenance_step!(sim,b)
add_callback!(sim,b,IterationInterval(1);name=:maintenance_terms,callsite=TendencyCallsite())
add_callback!(sim,s->maintenance_step!(s,b),IterationInterval(10);name=:adaptive_step)
add_callback!(sim,output!,TimeInterval(.25);name=:saved_diagnostics)
add_callback!(sim,checkpoint!,TimeInterval(5.);name=:checkpoint)
add_callback!(sim,progress!,TimeInterval(2.5);name=:progress)
output!(sim);checkpoint!(sim);run!(sim);checkpoint!(sim)
complete=abs(m.clock.time-target)<1e-7
open(joinpath(out,"status.toml"),"w") do io;TOML.print(io,Dict("status"=>(complete ? "COMPLETED" : "PARTIAL_WALL_LIMIT"),"t"=>m.clock.time,"target"=>target,"elapsed_seconds"=>time()-started);sorted=true);end
complete||error("Partial checkpoint retained; no automatic continuation")
output!(sim);println("COMPLETED case=$case t=$(m.clock.time)");flush(stdout)
