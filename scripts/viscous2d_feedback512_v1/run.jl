include("model.jl")
out,source,target_arg,validation,case=ARGS;target=parse(Float64,target_arg)
gate=TOML.parsefile(validation);@assert gate["status"]=="PASS" && gate["code_sha256"]==source_hashes()
ispath(out)&&error("Refusing existing output");mkpath(out)
d=load(source);u,v=d["u"],d["v"];t0=Float64(d["t"]);cf=coefficients(case)
@assert size(u)==size(v)==(512,512) && 0<=t0<target<=400
initial_hash=haskey(d,"initial_sha256") ? d["initial_sha256"] : bytes2hex(sha256(read(source)))
@assert initial_hash==gate["initial_sha256"]
if t0==0
    @assert !haskey(d,"config") && bytes2hex(sha256(read(source)))==gate["initial_sha256"]
    case=="biased_minus" && ((u,v)=reflect_state(u,v))
else
    cfg=d["config"]
    @assert cfg["case"]==case && cfg["code_sha256"]==source_hashes()
    for (k,val) in [("a",cf.a),("b",cf.b),("c",cf.c),("f",cf.f),("viscosity",.0001),("dt_ceiling",.0025)]
        @assert cfg[k]==val
    end
end
m,b=setup_polynomial(u,v;case,t=t0)
config=Dict("case"=>case,"grid"=>[512,512],"domain"=>[2pi,2pi],"f"=>cf.f,"beta"=>0.,"a"=>cf.a,"b"=>cf.b,"c"=>cf.c,
 "viscosity"=>.0001,"forcing_band"=>[12.,16.],"allow_zonal_forcing"=>false,"selective_viscosity"=>0.,
 "velocity_drag"=>0.,"closure"=>"P0(b*f*zeta^2-c*zeta^3)","dt_ceiling"=>.0025,"advective_cfl_target"=>.3,
 "reaction_step_target"=>.1,"dt_floor"=>1e-5,"start_time"=>t0,"target_time"=>target,"final_campaign_time"=>400.,
 "frame_interval"=>.5,"checkpoint_interval"=>10.,"initial_sha256"=>initial_hash,"code_sha256"=>source_hashes(),
 "reflected_initial_condition"=>(case=="biased_minus"),"initial_mirror_axis"=>"x", "compute_location"=>"Unity",
 "interpretation"=>"Imposed polynomial closure; neither a convection model nor an established resonance mechanism")
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,file");end
open(joinpath(out,"scalars.csv"),"w") do io
 println(io,"t,dt,energy,enstrophy,mean_zeta,third_moment,skewness,energy_k_le_2,enstrophy_k_ge_128,quadratic_m3_rate,quadratic_m3_identity,cubic_m3_rate,quadratic_energy_rate,cubic_energy_rate,quadratic_lowk_energy_rate,cubic_lowk_energy_rate,linear_energy_rate")
end
started=time();lastframe=Ref(-Inf);lastcheckpoint=Ref(-Inf)
function checkpoint!(sim)
    abs(m.clock.time-lastcheckpoint[])<1e-8 && return
    u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),config,initial_sha256=initial_hash)
    mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true);lastcheckpoint[]=m.clock.time
end
function output!(sim)
    abs(m.clock.time-lastframe[])<1e-8 && return
    u,v=state(m);q=vorticity(b.op,u,v);all(isfinite,q)||error("Nonfinite field")
    ix=round(Int,2m.clock.time);@assert abs(ix*.5-m.clock.time)<1e-7
    name=@sprintf("zeta_%08d.f32",ix);open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
    open(joinpath(out,"frames.csv"),"a") do io;println(io,join((m.clock.iteration,m.clock.time,name),','));end
    q2=q.^2;q3=q.^3;qb=b.b_f.*q2;qb.-=mean(qb);qc=-b.c.*q3;qc.-=mean(qc)
    m2=mean(q2);m3=mean(q3);mb=3mean(q2.*qb);identity=3b.b_f*(mean(q2.^2)-m2^2)
    @assert isapprox(mb,identity;rtol=1e-10,atol=1e-12)
    qhat=FFTW.fft(q);bh=FFTW.fft(qb);ch=FFTW.fft(qc);n=length(q)^2
    invlambda=map(x->x==0 ? 0. : 1/x,b.op.lambda);low=(b.op.k2.>0).&(b.op.k2.<=4.)
    ek=0.5 .*abs2.(qhat).*invlambda./n
    workb=real.(conj.(qhat).*bh).*invlambda./n;workc=real.(conj.(qhat).*ch).*invlambda./n
    linwork=sum(abs2.(qhat).*b.linear.*invlambda)/n
    hi=.5*sum(abs2.(qhat)[b.op.k2.>=128^2])/n
    vals=(m.clock.time,sim.Δt,.5*(mean(abs2,u)+mean(abs2,v)),.5m2,mean(q),m3,m3/m2^1.5,
          sum(ek[low]),hi,mb,identity,3mean(q2.*qc),sum(workb),sum(workc),sum(workb[low]),sum(workc[low]),linwork)
    open(joinpath(out,"scalars.csv"),"a") do io;println(io,join(vals,','));end
    lastframe[]=m.clock.time
end
function progress!(sim)
 println(@sprintf("PROGRESS case=%s t=%.3f target=%.1f dt=%.6g elapsed=%.1fs",case,m.clock.time,target,sim.Δt,time()-started));flush(stdout)
end
sim=Simulation(m;Δt=.0025,stop_time=target,wall_time_limit=5400.,verbose=false)
set_stable_step!(sim,b)
add_callback!(sim,b,IterationInterval(1);name=:polynomial_terms,callsite=TendencyCallsite())
add_callback!(sim,s->set_stable_step!(s,b),IterationInterval(10);name=:adaptive_step)
add_callback!(sim,output!,TimeInterval(.5);name=:saved_fields)
add_callback!(sim,checkpoint!,TimeInterval(10.);name=:checkpoint)
add_callback!(sim,progress!,TimeInterval(2.5);name=:progress)
output!(sim);checkpoint!(sim);run!(sim);checkpoint!(sim)
complete=abs(m.clock.time-target)<1e-7
open(joinpath(out,"status.toml"),"w") do io;TOML.print(io,Dict("status"=>(complete ? "COMPLETED" : "PARTIAL_WALL_LIMIT"),"t"=>m.clock.time,"target"=>target,"elapsed_seconds"=>time()-started);sorted=true);end
complete||error("Partial checkpoint retained; dependent jobs remain held")
output!(sim);println("SEGMENT_COMPLETED case=$case t=$(m.clock.time)");flush(stdout)
