include("model.jl")
include("reference.jl")
out,source,target_arg,validation,case=ARGS;target=parse(Float64,target_arg)
gate=TOML.parsefile(validation);@assert gate["status"]=="PASS" && gate["code_sha256"]==source_hashes()
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u,v=d["u"],d["v"];t0=Float64(d["t"]);c=coefficients(case)
@assert size(u)==size(v)==(GRID,GRID) && 0<=t0<target<=400
initial_hash=d["initial_sha256"];@assert initial_hash==gate["initial_sha256"]
if t0==0
    @assert d["prepared_initial"] && bytes2hex(sha256(read(source)))==gate["prepared_initial_sha256"]
    case=="biased_minus" && ((u,v)=reflect_state(u,v))
else
    @assert d["config"]["case"]==case && d["config"]["code_sha256"]==source_hashes()
end
m,b=setup_stress(u,v;case,t=t0)
cfg=Dict("case"=>case,"grid"=>[GRID,GRID],"domain"=>[2pi,2pi],"f"=>c.f,"beta"=>0.,
 "nu0"=>c.nu0,"nu1"=>c.nu1,"nu2"=>c.nu2,"viscosity"=>VISCOSITY,"cutoff"=>KCUT,
 "dt_ceiling"=>DTMAX,"frame_interval"=>1.,"start_time"=>t0,"target_time"=>target,"initial_sha256"=>initial_hash,
 "code_sha256"=>source_hashes(),"closure"=>"Rdev=P[-2(nu0+nu1*Zbar/f+nu2*(Zbar/f)^2) Sbar]; force=-div Rdev",
 "projection"=>"orthogonal radial Fourier low-pass, all modes |k|<=16 including zonal; input and output projected",
 "center_vorticity"=>"adjacent four-corner arithmetic mean; shear stress at corners; normal stresses at centers",
 "linear_drag"=>0.,"hyperviscosity"=>0.,"additional_instability_forcing"=>0.,"compute_location"=>"Unity",
 "interpretation"=>"Minimal constitutive hypothesis; coefficients not derived from or fitted to Breeze",
 "reflected_initial_condition"=>(case=="biased_minus"))
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,cfg;sorted=true);end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,file");end
started=time();lastframe=Ref(-Inf);lastcheckpoint=Ref(-Inf);keysout=sort(collect(keys(diagnostics(u,v,c))))
open(joinpath(out,"scalars.csv"),"w") do io;println(io,join(vcat(["t","dt"],keysout),','));end
function checkpoint!(sim)
    abs(m.clock.time-lastcheckpoint[])<1e-8 && return
    u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),config=cfg,initial_sha256=initial_hash)
    mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true);lastcheckpoint[]=m.clock.time
end
function output!(sim)
    abs(m.clock.time-lastframe[])<1e-8 && return
    u,v=state(m);q=vorticity(b.op,u,v);all(isfinite,q)||error("Nonfinite field")
    ix=round(Int,m.clock.time);@assert abs(ix-m.clock.time)<1e-7
    name=@sprintf("zeta_%08d.f32",ix);open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
    open(joinpath(out,"frames.csv"),"a") do io;println(io,join((m.clock.iteration,m.clock.time,name),','));end
    ds=diagnostics(u,v,c)
    open(joinpath(out,"scalars.csv"),"a") do io;println(io,join(vcat([m.clock.time,sim.Δt],[ds[k] for k in keysout]),','));end
    lastframe[]=m.clock.time
end
function progress!(sim)
    println(@sprintf("PROGRESS case=%s t=%.2f target=%.1f dt=%.6g elapsed=%.1fs filtered_max_Ro=%.4f",case,m.clock.time,target,sim.Δt,time()-started,b.max_rossby));flush(stdout)
end
sim=Simulation(m;Δt=DTMAX,stop_time=target,wall_time_limit=3300.,verbose=false);stable!(sim,b)
add_callback!(sim,b,IterationInterval(1);name=:stress,callsite=TendencyCallsite())
add_callback!(sim,s->stable!(s,b),IterationInterval(10);name=:adaptive)
add_callback!(sim,output!,TimeInterval(1.);name=:output)
add_callback!(sim,checkpoint!,TimeInterval(10.);name=:checkpoint)
add_callback!(sim,progress!,TimeInterval(5.);name=:progress)
output!(sim);checkpoint!(sim)
try
    run!(sim)
catch e
    checkpoint!(sim)
    rethrow()
end
checkpoint!(sim);complete=abs(m.clock.time-target)<1e-7
open(joinpath(out,"status.toml"),"w") do io;TOML.print(io,Dict("status"=>(complete ? "COMPLETED" : "PARTIAL_WALL_LIMIT"),"t"=>m.clock.time,"target"=>target,"elapsed_seconds"=>time()-started);sorted=true);end
complete||error("Partial output retained; no automatic continuation")
output!(sim);println("COMPLETED case=$case t=$(m.clock.time)");flush(stdout)
