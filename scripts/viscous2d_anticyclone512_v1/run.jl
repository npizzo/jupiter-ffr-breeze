include("model.jl")
out,source,target_str,validation=ARGS;target=parse(Float64,target_str);case="anticyclone_diffusion"
gate=TOML.parsefile(validation)
@assert gate["status"]=="PASS" && gate["selective_viscosity"]==.0004
@assert gate["code_sha256"]==source_hashes()
ispath(out)&&error("Refusing existing output");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];t0=Float64(d["t"]);dt=.0025
@assert size(u0)==size(v0)==(512,512) && 0<=t0<target<=400
it0=round(Int,t0/dt);itend=round(Int,target/dt)
@assert abs(it0*dt-t0)<1e-7 && abs(itend*dt-target)<1e-10
if haskey(d,"config")
    c=d["config"]
    for (k,val) in [("case",case),("dt",dt),("viscosity",.0001),("drag",0.),("f",1.),("beta",0.),("instability_coefficient",0.),("selective_viscosity",.0004)]
        @assert c[k]==val "Restart configuration mismatch: $k"
    end
end
initial_hash=haskey(d,"initial_sha256") ? d["initial_sha256"] : bytes2hex(sha256(read(source)))
@assert initial_hash==gate["initial_sha256"]
if t0==0
    @assert bytes2hex(sha256(read(source)))==gate["initial_sha256"]
else
    @assert haskey(d,"config") && d["config"]["code_sha256"]==source_hashes()
end
m,b=setup_anticyclone(u0,v0;t=t0)
config=Dict("case"=>case,"grid"=>[512,512],"domain"=>[2pi,2pi],"f"=>1.,"beta"=>0.,"dt"=>dt,
 "viscosity"=>.0001,"drag"=>0.,"instability_coefficient"=>0.,"selective_viscosity"=>.0004,
 "selective_operator"=>"nu_A Laplacian(sign(f) min(sign(f) zeta, 0))","start_time"=>t0,"target_time"=>target,
 "final_campaign_time"=>400.,"initial_sha256"=>initial_hash,"source_sha256"=>bytes2hex(sha256(read(source))),
 "compute_location"=>"Unity","implementation"=>"Native ScalarDiffusivity plus local C-grid rotated gradient of clipped vorticity; no forcing FFTs","outputs"=>"Checkpoints, raw vorticity and inexpensive scalar dissipation diagnostics; downstream movie rendering on Unity",
 "frame_interval"=>.5,"checkpoint_interval"=>10.,"internal_wall_seconds"=>5400.,
 "code_sha256"=>source_hashes())
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,file");end
open(joinpath(out,"scalars.csv"),"w") do io
    println(io,"t,energy,enstrophy,mean_zeta,skewness,cyclonic_enstrophy,anticyclonic_enstrophy,selective_energy_rate,selective_enstrophy_rate")
end
started=time();lastframe=Ref(-1);lastcheckpoint=Ref(-1)
function checkpoint!()
    u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),iteration=it0+m.clock.iteration,config,initial_sha256=initial_hash)
    mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true)
end
function output!(sim;force=false)
    it=it0+m.clock.iteration
    if (it%200==0 || force) && lastframe[]!=it
        u,v=state(m);q=vorticity(b.op,u,v);all(isfinite,q)||error("Nonfinite saved field")
        name=@sprintf("zeta_%08d.f32",it)
        open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
        open(joinpath(out,"frames.csv"),"a") do io;println(io,join((it,m.clock.time,name),','));end
        qa=min.(q,0.);qc=max.(q,0.);dx=b.op.lx/b.op.nx;dy=b.op.ly/b.op.ny
        # Periodic summation by parts: exact semi-discrete selective work.
        zwork=-b.nu_A*(mean((circshift(q,(-1,0)).-q).*(circshift(qa,(-1,0)).-qa))/dx^2+
                      mean((circshift(q,(0,-1)).-q).*(circshift(qa,(0,-1)).-qa))/dy^2)
        ework=-b.nu_A*mean(abs2,qa)
        @assert ework<=0 && zwork<=0
        vals=(m.clock.time,.5*(mean(abs2,u)+mean(abs2,v)),.5*mean(abs2,q),mean(q),
              mean(q.^3)/mean(abs2,q)^1.5,.5*mean(abs2,qc),.5*mean(abs2,qa),ework,zwork)
        open(joinpath(out,"scalars.csv"),"a") do io;println(io,join(vals,','));end
        lastframe[]=it
    end
    if (it%4000==0 || force) && lastcheckpoint[]!=it
        checkpoint!();lastcheckpoint[]=it
    end
    if it>it0 && it%1000==0
        println(@sprintf("PROGRESS case=%s t=%.3f segment_target=%.1f elapsed=%.1fs",case,m.clock.time,target,time()-started));flush(stdout)
    end
end
sim=Simulation(m;Δt=dt,stop_iteration=itend-it0,wall_time_limit=5400.,verbose=false)
add_callback!(sim,b,IterationInterval(1);name=:selective_diffusion,callsite=TendencyCallsite())
add_callback!(sim,output!,IterationInterval(1);name=:saved_fields)
output!(sim;force=true);run!(sim);output!(sim;force=true)
complete=m.clock.iteration>=itend-it0 && abs(m.clock.time-target)<1e-7
open(joinpath(out,"status.toml"),"w") do io
 TOML.print(io,Dict("status"=>(complete ? "COMPLETED" : "PARTIAL_WALL_TIME_LIMIT"),"t"=>m.clock.time,"target"=>target,"elapsed_seconds"=>time()-started);sorted=true)
end
complete||error("Partial output retained; downstream jobs stay stopped")
println("SEGMENT_COMPLETED case=$case t=$(m.clock.time)");flush(stdout)
