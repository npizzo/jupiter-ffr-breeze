include("model.jl")
include("run.jl")
config()=TOML.parsefile(joinpath(@__DIR__,"config.toml"))
function small(c)
    τ=1/(c["rossby"]*c["f"])
    merge(c,Dict("nx"=>32,"ny"=>32,"nz"=>32,"turnovers"=>.04/τ,"initial_dt"=>.001,"diagnostic_turnovers"=>.01/τ,"movie_turnovers"=>.01/τ,"full_turnovers"=>.02/τ,"checkpoint_turnovers"=>.02/τ,"wall_seconds"=>300.,"test"=>true))
end
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>=2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    fs=filter(x->startswith(x,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        for k in ("vorticity_slices","tracer_slices","vertical_velocity_slices","buoyancy_slices");@assert all(isfinite,a[k]);end
    end
    fs=filter(x->startswith(x,"full_"),readdir(out));@assert !isempty(fs)
    a=load(joinpath(out,last(sort(fs))));@assert maximum(abs,a["w_face"][:,:,[1,end]])==0
end
function exercise(c,ratio,out,arch)
    τ=1/(c["rossby"]*c["f"])
    a,ma=run_case(c,ratio,joinpath(out,"direct");arch)
    _,mf=run_case(merge(c,Dict("turnovers"=>.02/τ)),ratio,joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first","checkpoints");join=true))
    b,mb=run_case(c,ratio,joinpath(out,"restart");arch,pickup=cp)
    @assert ma["status"]==mf["status"]==mb["status"]=="COMPLETE"
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k]))-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-9
    for part in ("direct","first","restart");readback(joinpath(out,part));end
    err
end
function initial_checks(c,out)
    result=Dict{String,Any}();common=nothing
    for (label,ratio) in (("weak",.5),("strong",2.))
        m,p,info=build(c,ratio;arch=CPU());a=fields_cpu(m)
        if common===nothing
            common=(copy(a.uf),copy(a.vf))
            info["velocity_match_relative_error"]=0.
        else
            # Native velocities are reconstructed from rho*u / rho, so allow
            # roundoff from the two slightly different reference densities.
            err=max(maximum(abs,a.uf-common[1]),maximum(abs,a.vf-common[2]))/p.seed.U
            info["velocity_match_relative_error"]=err
            println("VELOCITY_MATCH relative_error=$err");flush(stdout)
            @assert err<1e-12
        end
        result[label]=info
        # A short actual-model step quantifies imbalance using the projected velocity.
        sim=Simulation(m;Δt=.001,stop_time=.01);run!(sim);aa=fields_cpu(m)
        info["initial_adjustment_w_over_U"]=maximum(abs,aa.w)/p.seed.U
        @assert info["initial_adjustment_w_over_U"]<.01
        @assert constraint(m,p)<c["constraint_tolerance"]
        GC.gc()
    end
    write_toml(joinpath(out,"INITIAL_PAIR_PASS.toml"),result)
    result
end
function continuation(c,ratio,out,arch;test=false)
    label=ratio==.5 ? "weak" : "strong"
    τ=1/(c["rossby"]*c["f"])
    if test
        # Same production continuation function, with a freshly serialized short parent.
        sc=small(c)
        _,meta=run_case(merge(sc,Dict("turnovers"=>.02/τ)),ratio,joinpath(out,"seed_checkpoint");arch)
        cp=only(readdir(joinpath(out,"seed_checkpoint","checkpoints");join=true))
        start=.02;pilotend=.03;finish=.04;active=sc
    else
        cp=c[label*"_checkpoint"];@assert isfile(cp)
        start=25τ;pilotend=30τ;finish=100τ;active=c
    end
    # Five-turnover pilot measures actual continuation cost before the longer segment.
    _,pilot=run_case(merge(active,Dict("expected_start_time"=>start,"turnovers"=>pilotend/τ)),ratio,joinpath(out,"pilot");arch,pickup=cp)
    @assert pilot["status"]=="COMPLETE"
    readback(joinpath(out,"pilot"))
    projected=pilot["wall_seconds_actual"]*(finish-pilotend)/(pilotend-start)
    write_toml(joinpath(out,"PACE.toml"),Dict("pilot_wall_seconds"=>pilot["wall_seconds_actual"],"projected_remaining_seconds"=>projected,"allowed_remaining_seconds"=>2400,"production_released"=>test || projected<2400))
    @assert test || projected<2400 "Continuation unexpectedly expensive; hold for review"
    cp2=only(readdir(joinpath(out,"pilot","checkpoints");join=true))
    _,meta=run_case(merge(active,Dict("expected_start_time"=>pilotend,"turnovers"=>finish/τ)),ratio,joinpath(out,"parent");arch,pickup=cp2)
    @assert meta["status"]=="COMPLETE"
    readback(joinpath(out,"parent"))
    write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"final_turnovers"=>meta["final_turnovers"],"extension_to_200"=>"Requires assessment of completed 100-turnover pair"))
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out);c=config()
    @assert c["nx"]==c["ny"]==128 && c["nz"]==64 && c["turnovers"]==100
    if mode=="preflight"
        result=Dict{String,Any}("status"=>"PASS")
        for (label,ratio) in (("weak",.5),("strong",2.))
            result[label*"_restart_error"]=exercise(small(c),ratio,joinpath(out,label*"_equivalence"),CPU())
            continuation(c,ratio,joinpath(out,label*"_path"),CPU();test=true)
        end
        write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    elseif mode=="run"
        @assert CUDA.functional();ratio=parse(Float64,ARGS[3]);@assert ratio in (.5,2.)
        err=exercise(small(c),ratio,joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        GC.gc();CUDA.reclaim()
        continuation(c,ratio,out,GPU())
    else
        error("Unknown mode")
    end
end
main()
