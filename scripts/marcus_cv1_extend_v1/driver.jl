include("model.jl")

mutable struct SavedStep
    dt::Float64
end
(s::SavedStep)(sim)=(s.dt=sim.Δt;nothing)
Oceananigans.prognostic_state(s::SavedStep)=(;dt=s.dt)
Oceananigans.restore_prognostic_state!(s::SavedStep,v::NamedTuple)=(s.dt=v.dt;s)

function config(;small=false)
    c=fixture()
    merge!(c,Dict("nx"=>small ? 24 : 384,"ny"=>small ? 24 : 384,"nz"=>small ? 64 : 256,
        "fixture_only"=>small,"reconstruction_authorized"=>true,"sponge_depth"=>31e3,
        "dt"=>small ? .5 : 2.,"max_dt"=>small ? .5 : 4.,"stop_time"=>small ? 4. : 4*86400.,
        "diagnostic_interval"=>small ? 1. : 2160.,"checkpoint_interval"=>small ? 2. : 21600.,
        "wall_seconds"=>small ? 600. : 42000.,"max_mach"=>.8))
    c
end
sync(m)=m.architecture isa GPU ? CUDA.synchronize() : nothing
function geometry_checks(c,info,p)
    dz=(c["zmax"]-c["zmin"])/c["nz"]
    rate=Array(interior(p.rate));rho=reshape(info.rho,1,1,:)
    uv2=info.flow.uf .^ 2 .+ info.flow.vf .^ 2
    mask=rate .> 0
    overlap=sum(rho .* ((info.flow.uf .- info.flow.ubg) .^ 2 .+ info.flow.vf .^ 2) .* mask)/max(sum(rho .* ((info.flow.uf .- info.flow.ubg) .^ 2 .+ info.flow.vf .^ 2)),1e-30)
    N2=c["g"] .* diff(log.(info.θref)) ./ dz
    Dict("vortex_energy_fraction_in_sponge"=>overlap,"max_N_squared"=>maximum(N2),"min_reference_density"=>minimum(info.rho),
         "max_reference_density"=>maximum(info.rho),"removed_mean_vorticity_max"=>maximum(abs,info.removed),"sponge_depth_m"=>c["sponge_depth"])
end
function sample(m,c,p,info,out;write_frame=true,full=false)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5 .* (uf .+ circshift(uf,(-1,0,0)));v=.5 .* (vf .+ circshift(vf,(0,-1,0)));w=.5 .* (wf[:,:,1:end-1] .+ wf[:,:,2:end])
    ρ=reshape(info.rho,1,1,:);θ=Array(interior(m.formulation.potential_temperature));T=Array(interior(m.temperature));dye=Array(interior(m.tracers.dye)) ./ ρ
    dx=c["lx"]/c["nx"];dy=c["ly"]/c["ny"];dz=(c["zmax"]-c["zmin"])/c["nz"];dv=dx*dy*dz
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv));mz=Array(interior(m.momentum.ρw))
    div=(circshift(mx,(-1,0,0)) .- mx)/dx+(circshift(my,(0,-1,0)) .- my)/dy+diff(mz;dims=3)/dz
    constraint=maximum(abs,div)/max(maximum(abs,mx)/dx,maximum(abs,my)/dy,maximum(abs,mz)/dz,1e-12)
    f=2c["rotation_rate"]*sind(c["latitude"])
    ζf=(vf .- circshift(vf,(1,0,0))) ./ dx .- (uf .- circshift(uf,(0,1,0))) ./ dy
    ζ=.25 .* (ζf .+ circshift(ζf,(-1,0,0)) .+ circshift(ζf,(0,-1,0)) .+ circshift(ζf,(-1,-1,0)))
    rate=Array(interior(p.rate));Ub=Array(interior(p.ubg));Thb=Array(interior(p.thbg));Pi=reshape(vec(Array(interior(m.dynamics.reference_state.pressure))),1,1,:) ./ c["p0"]
    Pi=Pi .^ (c["R"]/c["cp"])
    Eh=.5dv*sum(ρ .* (u .^ 2 .+ v .^ 2));Ev=.5dv*sum(ρ .* w .^ 2)
    E=dv*sum(ρ .* (.5 .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .+ c["cp"] .* Pi .* θ))
    sponge_power=-dv*sum(rate .* ρ .* (u .* (u .- Ub) .+ v .^ 2 .+ w .^ 2 .+ c["cp"] .* Pi .* (θ .- Thb)))
    mach=maximum(sqrt.(u .^ 2 .+ v .^ 2 .+ w .^ 2) ./ sqrt.((c["cp"]/(c["cp"]-c["R"]))*c["R"] .* T))
    delta=maximum(abs,θ ./ reshape(info.θref,1,1,:) .- 1)
    k=[argmin(abs.(info.z .- z0)) for z0 in (-40e3,c["z0"],50e3)]
    # Core circulation is measured in a fixed ellipse at the centre height, including background shear.
    kk=k[2];mask=[(((-c["lx"]/2+(i-.5)*dx)/(c["epsilon"]*c["Rv"]))^2+((-c["ly"]/2+(j-.5)*dy-info.z[kk]*cotd(c["latitude"]))/c["Rv"])^2)<=1 for i=1:c["nx"],j=1:c["ny"]]
    circulation=sum(ζ[:,:,kk] .* mask)*dx*dy
    vals=(Float64(m.clock.time),maximum(abs,w),maximum(abs,u),constraint,minimum(T),E,sponge_power,Eh,Ev,mach,delta,minimum(ζ[:,:,kk]),maximum(ζ[:,:,kk]),circulation,minimum(dye),maximum(dye))
    @assert all(isfinite,vals) && minimum(T)>0 && constraint<1e-8
    @assert mach<c["max_mach"] "Mach guard"
    open(joinpath(out,"diagnostics.csv"),"a") do io;println(io,join(vals,','));end
    if write_frame
        jldsave(joinpath(out,@sprintf("frame_%012d.jld2",round(Int,m.clock.time)));t=Float64(m.clock.time),z=info.z,lx=c["lx"],ly=c["ly"],f,
            video_heights=info.z[k],height_coordinate="metres relative to one-bar level",vorticity_slices=Float32.(ζ[:,:,k]),tracer_slices=Float32.(dye[:,:,k]),
            temperature_slices=Float32.(T[:,:,k]),vertical_velocity_slices=Float32.(w[:,:,k]),theta_slices=Float32.(θ[:,:,k]),
            u_slices=Float32.(u[:,:,k]),v_slices=Float32.(v[:,:,k]),
            vorticity_xz=Float32.(ζ[:,c["ny"]÷2,:]),theta_xz=Float32.(θ[:,c["ny"]÷2,:]),
            temperature_xz=Float32.(T[:,c["ny"]÷2,:]),vertical_velocity_xz=Float32.(w[:,c["ny"]÷2,:]),
            section_y=-c["ly"]/2+(c["ny"]÷2-.5)*dy)
    end
    if full
        jldsave(joinpath(out,@sprintf("full_%012d.jld2",round(Int,m.clock.time)));t=Float64(m.clock.time),z=info.z,u_face=Float32.(uf),v_face=Float32.(vf),w_face=Float32.(wf),temperature=Float32.(T),theta=Float32.(θ))
    end
    @printf("CV1 t=%.1fs days=%.6f max_w=%.5g Mach=%.4g constraint=%.3g vertical_KE_fraction=%.4g\n",m.clock.time,m.clock.time/86400,vals[2],mach,constraint,Ev/max(Eh,1e-30));flush(stdout)
    vals
end
function run_case(c,out;arch=CPU(),vortex=true,rest=false,pickup=false,pilot=false)
    @assert !ispath(out);mkpath(out)
    started=time();m,p,info=build(c;arch,vortex,rest);sync(m)
    audit=geometry_checks(c,info,p);N=sqrt(max(audit["max_N_squared"],1e-12));cap=min(c["max_dt"],.1/N)
    sim=Simulation(m;Δt=min(c["dt"],cap),stop_time=c["stop_time"],wall_time_limit=c["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpoint]=Checkpointer(m;schedule=TimeInterval(c["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt
        Oceananigans.TimeSteppers.update_state!(m)
        @assert m.clock.time==c["restart_time"] && m.clock.iteration==c["restart_iteration"]
    end
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    open(joinpath(out,"diagnostics.csv"),"w") do io;println(io,"time,max_w,max_u,mass_flux_divergence_relative,min_temperature,diagnostic_energy,sponge_power,horizontal_ke,vertical_ke,max_mach,max_relative_theta_departure,zeta_min_centre,zeta_max_centre,core_ellipse_circulation,tracer_min,tracer_max");end
    initial=sample(m,c,p,info,out;full=c["fixture_only"])
    # Explicit reset after build/output: timings exclude initialization and first40steps.
    rows=Vector{Vector{Float64}}();lw=Ref(0.);lt=Ref(0.);li=Ref(0)
    function timer(s)
        sync(m);now=time();it=m.clock.iteration;t=Float64(m.clock.time)
        if lw[]>0;push!(rows,[it,t,now-lw[],t-lt[],it-li[],s.Δt]);end
        lw[]=now;lt[]=t;li[]=it
    end
    add_callback!(sim,timer,IterationInterval(20);name=:timing)
    last=Ref(initial)
    add_callback!(sim,s->(last[]=sample(m,c,p,info,out)),TimeInterval(c["diagnostic_interval"]);name=:output)
    meta=Dict{String,Any}("status"=>"RUNNING","parameters"=>c,"audit"=>audit,"start_time"=>Float64(m.clock.time),"max_timestep"=>cap,
        "reconstruction"=>true,"exact_author_input"=>false,"energy_note"=>"Cell-centred anelastic diagnostic proxy; not an exact discrete conservation identity. Instantaneous sponge power recorded.",
        "initial_core_circulation"=>initial[14],"initial_KE_horizontal"=>initial[8],"scientific_question"=>"Persistence of an initialized CV1-like anticyclone, not spontaneous sign selection")
    write_toml(joinpath(out,"metadata.toml"),meta)
    try
        run!(sim;checkpoint_at_end=true)
        final=sample(m,c,p,info,out;full=c["fixture_only"])
        meta["status"]=m.clock.time>=c["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
        meta["final_core_circulation"]=final[14];meta["final_KE_horizontal"]=final[8];meta["final_vertical_ke_fraction"]=final[9]/max(final[8],1e-30)
        meta["final_mach"]=final[10];meta["final_theta_departure"]=final[11]
    catch e
        meta["status"]="FAILED";meta["reason"]=sprint(showerror,e);rethrow()
    finally
        meta["final_time"]=Float64(m.clock.time);meta["wall_seconds_actual"]=time()-started
        good=[r for r in rows if r[1]-r[5]>=40 && r[4]>0]
        if !isempty(good)
            meta["timed_seconds_per_model_second"]=sum(r[3] for r in good)/sum(r[4] for r in good)
            meta["estimated_four_day_wall_hours"]=meta["timed_seconds_per_model_second"]*4*86400/3600
        end
        write_toml(joinpath(out,"metadata.toml"),meta)
        open(joinpath(out,"timing.csv"),"w") do io;println(io,"iteration,time,wall_seconds,model_seconds,steps,dt");for r in rows;println(io,join(r,','));end;end
    end
    sim,meta
end
function readback(out;full=true)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(s,','))) for s in rows[2:end])
    for f in filter(x->startswith(x,"frame_"),readdir(out))
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        for key in ("vorticity_slices","tracer_slices","temperature_slices","vertical_velocity_slices","u_slices","v_slices","vorticity_xz","theta_xz","temperature_xz","vertical_velocity_xz");@assert all(isfinite,a[key]);end
    end
    if full
        for f in filter(x->startswith(x,"full_"),readdir(out));a=load(joinpath(out,f));@assert maximum(abs,a["w_face"][:,:,[1,end]])==0;end
    end
end
function exercise(out,arch)
    c=config(small=true);results=Dict{String,Any}()
    for (name,vortex,rest) in (("rest",false,true),("shear",false,false),("vortex",true,false))
        a,meta=run_case(c,joinpath(out,name);arch,vortex,rest);@assert meta["status"]=="COMPLETE";readback(joinpath(out,name))
        if rest
            speed=maximum(maximum(abs,Array(interior(f))) for f in values(a.model.velocities));@assert speed<1e-8;results["rest_max_speed"]=speed
        end
        if vortex
            first,_=run_case(merge(c,Dict("stop_time"=>2.)),joinpath(out,"first");arch)
            cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
            b,_=run_case(merge(c,Dict("restart_time"=>2.,"restart_iteration"=>first.model.clock.iteration)),joinpath(out,"restart");arch,pickup=cp)
            fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
            err=maximum(maximum(abs,Array(interior(fa[k])) .- Array(interior(fb[k])))/max(maximum(abs,Array(interior(fa[k]))),1.) for k in keys(fa));@assert err<1e-9
            results["restart_relative_error"]=err;readback(joinpath(out,"first"));readback(joinpath(out,"restart"))
        end
        a=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    end
    results
end
function native_buoyancy_check(c)
    m,p,info=build(c;arch=CPU(),rest=true);errs=Float64[];diffs=Float64[]
    for δ in (-.05,-.01,.01,.05)
        set!(m;θ=reshape(info.θref,1,1,:) .* ones(c["nx"],c["ny"],1) .* (1+δ),qᵗ=0.)
        Oceananigans.TimeSteppers.update_state!(m)
        F=Oceananigans.fields(m);k=c["nz"]÷2
        bn=Breeze.AtmosphereModels.buoyancy_forceᶜᶜᶜ(1,1,k,m.grid,m.dynamics,m.temperature,zeros(c["nx"],c["ny"],c["nz"]),m.microphysics,m.microphysical_fields,m.thermodynamic_constants)
        bc=buoyancy_correction(1,1,k,m.grid,F,p);target=c["g"]*info.rho[k]*δ
        push!(errs,abs(bn+bc-target)/abs(target));push!(diffs,abs(bn-target)/abs(target))
    end
    @assert maximum(errs)<1e-10
    Dict("corrected_buoyancy_max_relative_error"=>maximum(errs),"native_buoyancy_max_relative_difference"=>maximum(diffs),"anomaly_samples"=>[-.05,-.01,.01,.05])
end


function case_config(which)
    @assert which in ("baseline", "fine")
    c=config();c["stop_time"]=(which=="baseline" ? 12. : 2.)*86400
    if which=="fine";c["nx"]=512;c["ny"]=512;end
    c["case"]=which;c
end
function resume_segment(c,out,request;arch)
    r=TOML.parsefile(request)
    @assert r["target_time"]==c["stop_time"] && r["restart_time"]<c["stop_time"]
    cfg=merge(c,Dict("restart_time"=>r["restart_time"],"restart_iteration"=>r["restart_iteration"]))
    sim,meta=run_case(cfg,joinpath(out,"parent");arch,pickup=r["checkpoint"])
    @assert meta["final_time"]>r["restart_time"]
    cp=only(filter(p->endswith(p,".jld2"),readdir(joinpath(out,"parent/checkpoints");join=true)))
    write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],
        "final_iteration"=>sim.model.clock.iteration,"checkpoint"=>abspath(cp),"start_time"=>r["restart_time"],"target_time"=>c["stop_time"]))
    sim,meta
end
function continuation_check(out,arch)
    c=config(small=true)
    a,_=run_case(merge(c,Dict("stop_time"=>2.)),joinpath(out,"seed");arch)
    cp=only(readdir(joinpath(out,"seed/checkpoints");join=true));req=joinpath(out,"request.toml")
    write_toml(req,Dict("checkpoint"=>abspath(cp),"restart_time"=>2.,"restart_iteration"=>a.model.clock.iteration,"target_time"=>4.))
    b,_=resume_segment(c,joinpath(out,"resumed"),req;arch)
    direct,_=run_case(c,joinpath(out,"direct");arch)
    fa=Oceananigans.prognostic_fields(direct.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])) .- Array(interior(fb[k])))/max(maximum(abs,Array(interior(fa[k]))),1.) for k in keys(fa))
    @assert err<1e-9;readback(joinpath(out,"resumed/parent"));readback(joinpath(out,"direct"))
    @assert TOML.parsefile(joinpath(out,"resumed/OUTCOME.toml"))["status"]=="COMPLETE"
    err
end
function main()
    mode,out=ARGS[1:2];@assert mode in ("preflight","resume","finepilot");@assert !ispath(out);mkpath(out)
    if mode=="preflight"
        r=exercise(out,CPU());r["continuation_restart_error"]=continuation_check(joinpath(out,"continuation"),CPU())
        r["buoyancy_check"]=native_buoyancy_check(config(small=true));r["status"]="PASS"
        @assert case_config("baseline")["stop_time"]==12*86400 && case_config("fine")["nx"]==512
        write_toml(joinpath(out,"DRIVER_PASS.toml"),r)
    else
        @assert CUDA.functional()
        err=continuation_check(joinpath(out,"gpu_checks"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        GC.gc();CUDA.reclaim()
        jobend=parse(Float64,strip(read(`python $(joinpath(@__DIR__,"allocation.py")) deadline`,String)))
        write_toml(joinpath(out,"ALLOCATION.toml"),Dict("live_deadline_epoch"=>jobend,"checked_at_epoch"=>time(),"environment_end_time"=>get(ENV,"SLURM_JOB_END_TIME","missing")))
        if mode=="resume"
            c=case_config(ARGS[4]);c["wall_seconds"]=jobend-time()-1200;@assert c["wall_seconds"]>60
            resume_segment(c,out,ARGS[3];arch=GPU())
        else
            c=case_config("fine");merge!(c,Dict("stop_time"=>1800.,"diagnostic_interval"=>300.,"checkpoint_interval"=>1800.,"wall_seconds"=>jobend-time()-600))
            sim,meta=run_case(c,joinpath(out,"parent");arch=GPU());@assert meta["status"]=="COMPLETE"
            readback(joinpath(out,"parent");full=false)
            @assert meta["final_vertical_ke_fraction"]<.3 && .25<meta["final_KE_horizontal"]/meta["initial_KE_horizontal"]<2 && meta["final_theta_departure"]<.25
            hours=meta["timed_seconds_per_model_second"]*2*86400/3600;@assert hours<16 "Two-day fine-grid projection exceeds16GPUhours"
            cp=only(readdir(joinpath(out,"parent/checkpoints");join=true))
            write_toml(joinpath(out,"FINE_RELEASE.toml"),Dict("status"=>"PASS","estimated_two_day_hours"=>hours,"checkpoint"=>abspath(cp),"restart_time"=>meta["final_time"],"restart_iteration"=>sim.model.clock.iteration))
        end
    end
    println("MARCUS_EXTENSION_DONE ",mode);flush(stdout)
end
main()
