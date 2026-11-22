function config(label)
    TOML.parsefile("configs/chan-rotation-response-v1/"*label*".toml")
end
function small_config(label)
    merge(config(label),Dict("nx"=>16,"ny"=>16,"nz"=>24,"restart_time"=>.2,"ramp_start"=>.2,"ramp_duration"=>.2,"stop_time"=>.6,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"field_interval"=>.2,"wall_seconds"=>180.))
end
function state_difference(a,b)
    fa=Oceananigans.prognostic_fields(a);fb=Oceananigans.prognostic_fields(b)
    @assert keys(fa)==keys(fb)
    maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
end
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>2
    rows=[parse.(Float64,split(l,',')) for l in lines[2:end]]
    @assert all(all(isfinite,r) for r in rows)
    @assert all(diff(first.(rows)).>=0)
    files=sort(filter(f->startswith(f,"frame_"),readdir(out)))
    @assert length(files)>1
    for file in (first(files),last(files))
        v=load(joinpath(out,file));@assert size(v["vorticity_slices"],3)==3
        @assert all(isfinite,v["tracer_slices"])
        @assert isfinite(v["omega_actual"])&&v["f"]>0
    end
    return rows
end
function ramp_operator_check(cfg)
    m,p=build_model(cfg;arch=CPU(),perturb=false)
    rng=MersenneTwister(97)
    M=(ρu=rand(rng,24,24,32),ρv=rand(rng,24,24,32),ρw=rand(rng,24,24,32))
    err=0.
    for q in [-.1,0.,.25,.5,.75,1.,1.1]
        t=p.ramp_start+q*p.ramp_duration;om=omega_at(t,p)
        r=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=om,latitude=cfg["latitude"]))
        for (native,extra) in [(Oceananigans.Coriolis.x_f_cross_U,rotation_x),(Oceananigans.Coriolis.y_f_cross_U,rotation_y),(Oceananigans.Coriolis.z_f_cross_U,rotation_z)]
            want=-native(6,6,6,m.grid,r,M)
            got=-native(6,6,6,m.grid,p.rotation,M)+extra(6,6,6,m.grid,(;time=t),M,p)
            err=max(err,abs(got-want))
        end
    end
    @assert err<1e-12 "Ramp differs from full native Coriolis operator"
    @assert omega_at(p.ramp_start,p)==3.
    @assert omega_at(p.ramp_start+.5p.ramp_duration,p)≈3.75
    @assert omega_at(p.ramp_start+p.ramp_duration,p)==4.5
    @assert omega_at(p.ramp_start+2p.ramp_duration,p)==4.5
    return err
end
function main()
    mode,out=ARGS[1:2];ispath(out)&&error("Refuse existing output");mkpath(out)
    oldcp=only(readdir("data/chan-eddies-v3-cpu-2089424/first/checkpoints";join=true))
    if mode=="preflight"
        energy_budget_notice(.11,1.1)
        c=small_config("control")
        operator_checks(merge(c,Dict("initial_rms"=>0.,"check_output"=>out)))
        op_error=ramp_operator_check(small_config("ramp"))
        base,mb=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"control");arch=CPU(),pickup=oldcp)
        ref,mref=ReferenceV4.run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"reference_v4");arch=CPU(),pickup=oldcp)
        control_error=state_difference(base.model,ref.model);@assert control_error<1e-10
        ramp=small_config("ramp")
        first,mfirst=run_case(merge(ramp,Dict("stop_time"=>.4)),joinpath(out,"ramp_first");arch=CPU(),pickup=oldcp)
        cp=only(readdir(joinpath(out,"ramp_first/checkpoints");join=true))
        resumed,mres=run_case(merge(ramp,Dict("restart_time"=>.4)),joinpath(out,"ramp_restart");arch=CPU(),pickup=cp)
        direct,mdirect=run_case(ramp,joinpath(out,"ramp_direct");arch=CPU(),pickup=oldcp)
        restart_error=state_difference(resumed.model,direct.model);@assert restart_error<1e-7
        initial_errors=Dict{String,Float64}()
        for key in ["vorticity_slices","tracer_slices","temperature_slices","ubar","vbar","density_profile"]
            a=load(joinpath(out,"control/frame_0000020.jld2"),key)
            b=load(joinpath(out,"ramp_direct/frame_0000020.jld2"),key)
            initial_errors[key]=maximum(abs,a.-b);@assert initial_errors[key]==0.
        end
        for name in ["control","ramp_first","ramp_restart","ramp_direct"]
            rr=readback(joinpath(out,name))
            if startswith(name,"ramp");@assert last(rr)[end-1]==4.5;end
        end
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","control_matches_v4_max_abs_error"=>control_error,"ramp_operator_max_abs_error"=>op_error,"ramp_restart_max_abs_error"=>restart_error,"identical_initial_field_errors"=>initial_errors,"scope"=>"Actual driver: full-vector ramp, exact shared initial fields, v4 control comparison, outputs and serialized restart"))
    elseif mode in ["control","ramp"]
        CUDA.functional()||error("GPU required")
        origin=TOML.parsefile(ARGS[3]);@assert origin["status"]=="PASS"
        checkpoint=origin["checkpoint"];@assert bytes2hex(open(sha256,checkpoint))==origin["checkpoint_sha256"]
        # GPU compilation and actual ramp/restart/output path, inside this allocation.
        _,sm=run_case(small_config(mode),joinpath(out,"gpu_check");pickup=oldcp)
        @assert sm["status"]=="COMPLETE";rr=readback(joinpath(out,"gpu_check"))
        @assert last(rr)[end-1]==config(mode)["target_omega"]
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"target_rotation"=>last(rr)[end-1]))
        GC.gc();CUDA.reclaim()
        _,meta=run_case(config(mode),joinpath(out,"parent");pickup=checkpoint)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_continuation"=>false))
    else
        error("Unknown mode")
    end
end
include("model.jl")
include("initialization.jl")
include("checks.jl")
module ReferenceV4
include("../chan_eddies_v4/model.jl")
include("../chan_eddies_v4/initialization.jl")
end
main()
