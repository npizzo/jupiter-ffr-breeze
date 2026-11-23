include("model.jl")
include("initialization.jl")
include("cap_initialization.jl")
include("guards.jl")
include("capture.jl")
include("native_audit.jl")
include("stages.jl")
include("eddy_and_source.jl")
include("probe.jl")
include("checks.jl")
include("extended.jl")
include("waves.jl")
include("hotspots.jl")
include("run.jl")
function readback(out)
    for name in ["diagnostics.csv","extrema.csv"]
        lines=readlines(joinpath(out,name));@assert length(lines)>2
        @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    end
    @assert length(readlines(joinpath(out,"hotspots.csv")))>2
    @assert length(readlines(joinpath(out,"vortex_metrics.csv")))>2
    @assert isfile(joinpath(out,"budget_reference.jld2"))
    for prefix in ["frame_","yz_","full_","profiles_"]
        fs=sort(filter(f->startswith(f,prefix),readdir(out)));@assert !isempty(fs)
        a=load(joinpath(out,last(fs)))
        if prefix=="frame_";@assert size(a["vorticity_slices"],3)==3 && all(isfinite,a["tracer_slices"]);end
        if prefix=="full_";@assert maximum(abs,a["v_face"][:,[1,end],:])==0 && maximum(abs,a["w_face"][:,:,[1,end]])==0;end
    end
end
function exercise(c,out,arch)
    run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,_=run_case(merge(c,Dict("replay_duration"=>.02,"health_reference_path"=>joinpath(out,"first","health_reference.jld2"),"budget_reference_path"=>joinpath(out,"first","budget_reference.jld2"))),joinpath(out,"restart");arch,pickup=cp)
    a,_=run_case(c,joinpath(out,"direct");arch)
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-7
    for x in ["first","restart","direct"];readback(joinpath(out,x));end
    err
end
function base_config()
    merge(TOML.parsefile(joinpath(@__DIR__,"pilot.toml")),TOML.parsefile(joinpath(@__DIR__,"candidate.toml")),Dict(
        "nx"=>32,"ny"=>64,"nz"=>80,"lx"=>4.,"initial_rms"=>.1,"perturbation"=>0.,"initial_dt"=>.001,"max_dt"=>1/120,
        "wall_seconds"=>180.,"diagnostic_interval"=>.05,"movie_interval"=>1.,"checkpoint_interval"=>10.,"save_3d"=>true,
        "mach_stop"=>.8,"diagnostic_extension"=>true,"temperature_stop"=>.5,"profile_drift_stop"=>.5,"theta_log_departure_stop"=>log(4.),
        "heat"=>1.,"conductivity"=>.0004*.03125/(GRAV/2.5),"cap_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"bottom_factor"=>1.,
        "jet_design"=>"cap","jet_speed"=>0.,"allow_guard_stop"=>true,"neutral_coriolis"=>true,"closure"=>"lilly","entropy_diffusion_factor"=>1.,
        "momentum_scheme"=>"weno5","scalar_scheme"=>"weno5","stop_time"=>20.,"enforce_energy"=>false,
        "wave_diagnostics"=>false,"eddy_tilt_correction"=>true,"mean_restoration"=>true,"sponge"=>false,"horizontal_sponge"=>false,
        "purpose"=>"Matched modest-shear response pilot; no long-time LSV claim from t100; physical parameters fixed throughout"))
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
    arch=mode=="preflight" ? CPU() : GPU();mode=="preflight" || CUDA.functional() || error("GPU required")
    base=base_config();proofs=Dict{String,Any}()
    for (label,speed) in [("nojet",0.),("jet",.015625)]
        c=merge(base,Dict("jet_speed"=>speed,"nx"=>8,"ny"=>16,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,
            "diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"allow_guard_stop"=>false))
        err=exercise(c,joinpath(out,"preflight",label),arch)
        proofs[label]=Dict("status"=>"PASS","restart_error"=>err)
        GC.gc();mode=="preflight" || CUDA.reclaim()
    end
    proofs["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),proofs)
    mode=="preflight" && return
    outcomes=Dict{String,Any}();passed=true
    for (label,speed) in [("nojet",0.),("jet",.015625)]
        println("BEGIN_CASE ",label," pilot target20");flush(stdout)
        sim,meta=run_case(merge(base,Dict("jet_speed"=>speed)),joinpath(out,label*"_pilot");arch)
        outcomes[label*"_pilot"]=meta;passed=(meta["status"]=="COMPLETE")&&passed
        write_toml(joinpath(out,"OUTCOMES.toml"),outcomes);sim=nothing;GC.gc();CUDA.reclaim()
    end
    write_toml(joinpath(out,"PAIR_GATE.toml"),Dict("pass"=>passed,"continuation_target"=>100.,"both_required"=>true,
        "energy_residual"=>"Recorded as exploratory fidelity limitation; no1percent energy gate; no compensating heat"))
    if passed
        for (label,speed) in [("nojet",0.),("jet",.015625)]
            parent=joinpath(out,label*"_pilot");cp=only(readdir(joinpath(parent,"checkpoints");join=true))
            cfg=merge(base,Dict("jet_speed"=>speed,"stop_time"=>100.,"replay_duration"=>80.,"wall_seconds"=>480.,
                "health_reference_path"=>joinpath(parent,"health_reference.jld2"),"budget_reference_path"=>joinpath(parent,"budget_reference.jld2")))
            println("BEGIN_CASE ",label," continuation target100");flush(stdout)
            sim,meta=run_case(cfg,joinpath(out,label*"_to100");arch,pickup=cp)
            outcomes[label*"_to100"]=meta;write_toml(joinpath(out,"OUTCOMES.toml"),outcomes);sim=nothing;GC.gc();CUDA.reclaim()
        end
    end
    write_toml(joinpath(out,"RELEASE.toml"),Dict("production_release"=>false,"pilot_pair_passed"=>passed,"status"=>"EXPLORATORY_PHYSICS_PAIR_REVIEW_REQUIRED"))
end
main()
