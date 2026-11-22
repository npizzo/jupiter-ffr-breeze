function config()
    c=TOML.parsefile("configs/omega3-long-v1.toml")
    @assert c["omega"]==3.0
    c
end
function small_config()
    merge(config(),Dict("nx"=>16,"ny"=>16,"nz"=>24,"restart_time"=>.2,"stop_time"=>.4,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"field_interval"=>.1,"wall_seconds"=>240.))
end
function difference(a,b)
    fa=Oceananigans.prognostic_fields(a);fb=Oceananigans.prognostic_fields(b)
    @assert keys(fa)==keys(fb)
    maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
end
function readback(out)
    ls=readlines(joinpath(out,"diagnostics.csv"));@assert length(ls)>2
    rows=[parse.(Float64,split(l,',')) for l in ls[2:end]]
    @assert all(all(isfinite,r) for r in rows)
    @assert all(diff(first.(rows)).>=0)
    frames=sort(filter(f->startswith(f,"frame_"),readdir(out)))
    @assert !isempty(frames)
    for name in frames
        f=load(joinpath(out,name));@assert size(f["vorticity_slices"],3)==3
        @assert all(isfinite,f["vorticity_slices"])&&all(isfinite,f["tracer_slices"])
    end
    fs=readdir(joinpath(out,"fields3d");join=true);@assert !isempty(fs)
    f=load(last(sort(fs)));@assert minimum(f["density"])>0&&all(isfinite,f["temperature"])
    return rows
end
function restart_pair(c,out,arch)
    first,_=run_case(merge(c,Dict("stop_time"=>.2)),joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    resumed,meta=run_case(c,joinpath(out,"restart");arch,pickup=cp)
    for name in ["first","restart"];readback(joinpath(out,name));end
    @assert meta["status"]=="COMPLETE"
    first,resumed,cp
end
function main()
    mode,out=ARGS[1:2];ispath(out)&&error("Refuse existing output");mkpath(out)
    if mode=="preflight"
        c=small_config()
        _,resumed,cp=restart_pair(c,out,CPU())
        direct,_=run_case(c,joinpath(out,"direct");arch=CPU())
        err=difference(resumed.model,direct.model);@assert err<1e-7
        old=merge(c,Dict("target_omega"=>3.,"ramp_start"=>.2,"ramp_duration"=>20.))
        ref,_=Previous.run_case(old,joinpath(out,"reference");arch=CPU(),pickup=cp)
        parity=difference(resumed.model,ref.model);@assert parity<1e-10
        readback(joinpath(out,"direct"))
        # Compare refreshed saved fields at matched times, excluding the old manual sample.
        frame_error=0.
        for name in filter(f->startswith(f,"frame_"),readdir(joinpath(out,"restart")))
            a=load(joinpath(out,"restart",name));b=load(joinpath(out,a["t"]<=.2 ? "first" : "reference",name))
            for k in ["vorticity_slices","tracer_slices","temperature_slices","ubar","vbar"]
                frame_error=max(frame_error,maximum(abs,a[k].-b[k]))
            end
        end
        @assert frame_error<1e-7
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","omega_constant"=>3.,"restart_max_abs_error"=>err,"previous_model_max_abs_error"=>parity,"refreshed_output_max_abs_error"=>frame_error,"scope"=>"Actual continuation driver advances, serializes/restarts, exercises movies/full fields/diagnostics; compares prognostic state to previous fixed-Omega=3 implementation"))
    elseif mode=="run"
        CUDA.functional()||error("GPU required")
        origin=TOML.parsefile("runs/omega3-long-v1/ORIGIN.toml")
        @assert origin["status"]=="PASS"&&origin["restart_time"]==500.
        cp=origin["checkpoint"];@assert bytes2hex(open(sha256,cp))==origin["checkpoint_sha256"]
        restart_pair(small_config(),joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"omega_constant"=>3.,"scope"=>"GPU serialized restart and output readback"))
        GC.gc();CUDA.reclaim()
        _,meta=run_case(config(),joinpath(out,"parent");pickup=cp)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_continuation"=>false))
    else
        error("Unknown mode")
    end
end
include("model.jl")
include("initialization.jl")
include("checks.jl")
module Previous
include("../chan_rotation_response_v1/model.jl")
include("../chan_rotation_response_v1/initialization.jl")
end
main()
