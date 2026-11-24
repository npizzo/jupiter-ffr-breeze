include("model.jl")
include("transfer.jl")
include("diagnostics.jl")
include("run.jl")
function config()
    c=merge(TOML.parsefile(joinpath(@__DIR__,"parent_config.toml")),Dict("equations"=>"compressible",
        "duration"=>10.,"wall_seconds"=>4800.,"diagnostic_interval"=>.1,"sparse_interval"=>1.,
        "movie_interval"=>.5,"checkpoint_interval"=>5.))
    @assert c["omega"]==3 && c["beta"]==0 && c["jet_speed"]==0 && c["closure"]=="none"
    c
end
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    fs=filter(f->startswith(f,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        @assert all(isfinite,a["vorticity_slices"]) && all(isfinite,a["tracer_slices"])
    end
    @assert !isempty(readdir(joinpath(out,"checkpoints")))
end
function exercise(out,arch)
    mkpath(out)
    c=merge(config(),Dict("nx"=>16,"ny"=>16,"nz"=>24,"max_dt"=>.001,"stop_time"=>.04,
        "diagnostic_interval"=>.01,"sparse_interval"=>.02,"movie_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>600.))
    m,p=build_model(c;arch=CPU())
    # Smooth nonzero flow exercises staggered momentum transfer, rather than only rest.
    set!(m;u=(x,y,z)->.001sinpi(x/2)*cospi(y/2),v=(x,y,z)->-.001cospi(x/2)*sinpi(y/2),
         w=(x,y,z)->.0001sinpi(z)*cospi(x/2))
    s=extract_state(m);sf=joinpath(out,"state.jld2");jldsave(sf;state=s)
    ident=state_error(s,mapped_state(s,(16,16,24)));@assert ident<1e-12
    # Refine all dimensions by 3/2, exactly as in production.
    target=(24,24,36);f=mapped_state(s,target)
    for key in (:rho,:rhotheta,:mx,:my,:dye)
        err=abs(mean(getproperty(f,key))-mean(getproperty(s,key)))
        @assert err<1e-11
    end
    ff=joinpath(out,"fine.jld2");jldsave(ff;state=f)
    fc=merge(c,Dict("nx"=>24,"ny"=>24,"nz"=>36,"max_dt"=>.001*2/3,"stop_time"=>.02))
    run_branch(fc,ff,joinpath(out,"fine");arch);readback(joinpath(out,"fine"))
    run_branch(merge(c,Dict("stop_time"=>.02)),sf,joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,bb=run_branch(c,sf,joinpath(out,"restart");arch,pickup=cp)
    a,ab=run_branch(c,sf,joinpath(out,"direct");arch)
    err=state_error(extract_state(a.model),extract_state(b.model));@assert err<1e-7
    herr=abs(ab.state.heat-bb.state.heat);@assert herr<1e-9
    @assert ab.state.E0==bb.state.E0
    for name in ("first","restart","direct");readback(joinpath(out,name));end
    Dict("status"=>"PASS","identity_transfer_error"=>ident,"restart_error"=>err,"restart_heat_error"=>herr)
end
function prepare(checkpoint,out)
    @assert !ispath(out);mkpath(out);c=config()
    m,p=build_model(c;arch=CPU());sim=Simulation(m;Δt=.001)
    saved=SavedStep(.001);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    set!(sim;checkpoint);Oceananigans.TimeSteppers.update_state!(m);s=extract_state(m);e=energy(m,p)
    donorE=e.E;donorKE=e.ke;donorT=e.a.T;donoru=e.a.u
    jldsave(joinpath(out,"baseline.jld2");state=s)
    # Verify exact field reinstallation before defining the common branch start.
    install_state!(m,s);ea=energy(m,p)
    installerr=state_error(s,extract_state(m));Terr=maximum(abs,ea.a.T.-donorT);uerr=maximum(abs,ea.a.u.-donoru)
    println("REINSTALL_AUDIT ", (;installerr,Terr,uerr,energy_change=ea.E-donorE));flush(stdout)
    @assert installerr<1e-11 && Terr<1e-11 && uerr<1e-11 && abs(ea.E-donorE)<1e-8
    m=nothing;sim=nothing;e=nothing;ea=nothing;donorT=nothing;donoru=nothing;GC.gc()
    target=(384,384,120);f=mapped_state(s,target)
    audit=Dict{String,Any}("source_time"=>s.t,"source_energy"=>donorE,"source_ke"=>donorKE,
        "same_grid_prognostic_error"=>installerr,"same_grid_temperature_error"=>Terr,"same_grid_velocity_error"=>uerr,
        "method"=>"Limited conservative piecewise-linear density/rhotheta/dye and horizontal-momentum dual cells; linear vertical face momentum; no thermal or velocity rebalance")
    for key in (:rho,:rhotheta,:mx,:my,:dye)
        x=mean(getproperty(f,key))-mean(getproperty(s,key));audit[string(key)*"_mean_change"]=x
        @assert abs(x)<1e-10*max(1.,maximum(abs,getproperty(s,key)))
    end
    @assert minimum(f.rho)>0 && minimum(f.rhotheta)>0 && maximum(abs,f.mz[:,:,[1,end]])==0
    jldsave(joinpath(out,"fine.jld2");state=f)
    s=nothing;GC.gc()
    fc=merge(c,Dict("nx"=>384,"ny"=>384,"nz"=>120))
    fm,fp=build_model(fc;arch=CPU());install_state!(fm,f);fe=energy(fm,fp)
    audit["fine_energy"]=fe.E;audit["transfer_energy_change"]=fe.E-donorE
    audit["fine_ke"]=fe.ke;audit["transfer_ke_fraction"]=(fe.ke-donorKE)/max(abs(donorKE),1e-30)
    audit["fine_max_mach"]=maximum(sqrt.(fe.a.u.^2 .+fe.a.v.^2 .+fe.a.w.^2)./sqrt.((5/3).*fe.a.T))
    @assert audit["fine_max_mach"]<2 && abs(audit["transfer_ke_fraction"])<.1
    # A large interpolation jump would invalidate this short-window comparison.
    @assert abs(audit["transfer_energy_change"])<.1*(16c["bottom_flux"]*c["duration"])
    audit["status"]="PASS";write_toml(joinpath(out,"TRANSFER_PASS.toml"),audit)
    println("PRODUCTION_TRANSFER_PASS ",audit);flush(stdout)
end
function main()
    mode,out=ARGS[1:2]
    if mode=="preflight"
        @assert !ispath(out);mkpath(out)
        result=exercise(joinpath(out,"exercise"),CPU());write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
        prepare(ARGS[3],joinpath(out,"prepared"))
    elseif mode=="run"
        @assert !ispath(out);mkpath(out);CUDA.functional()||error("GPU required")
        check=exercise(joinpath(out,"gpu_check"),GPU());check["gpu"]=CUDA.name(CUDA.device())
        write_toml(joinpath(out,"GPU_PASS.toml"),check);GC.gc();CUDA.reclaim()
        c=config();idx=parse(Int,ARGS[3]);prepared=ARGS[4]
        names=["baseline","refined","smaller_dt"];@assert idx in 0:2
        statefile=joinpath(prepared,idx==1 ? "fine.jld2" : "baseline.jld2")
        audit=TOML.parsefile(joinpath(prepared,"TRANSFER_PASS.toml"));@assert audit["status"]=="PASS"
        c["case"]=names[idx+1];c["stop_time"]=audit["source_time"]+c["duration"]
        if idx==1;c["nx"]=384;c["ny"]=384;c["nz"]=120;end
        if idx!=0;c["max_dt"]*=2/3;end
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+5100)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-180);@assert c["wall_seconds"]>60
        run_branch(c,statefile,joinpath(out,"branch");arch=GPU())
    else
        error("Unknown mode")
    end
end
main()
