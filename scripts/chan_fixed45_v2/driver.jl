function config()
    TOML.parsefile("configs/chan-fixed45-v2.toml")
end
function small_config()
    merge(config(),Dict("nx"=>32,"ny"=>32,"nz"=>40,"stop_time"=>.1,"diagnostic_interval"=>.05,"movie_interval"=>.05,"checkpoint_interval"=>.1,"field_interval"=>.1,"wall_seconds"=>240.))
end
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
    frames=sort(filter(f->startswith(f,"frame_"),readdir(out)));@assert length(frames)>1
    a=load(joinpath(out,last(frames)));@assert size(a["vorticity_slices"],3)==3
    @assert all(isfinite,a["tracer_slices"])
    fields=sort(readdir(joinpath(out,"fields3d");join=true));@assert !isempty(fields)
    f=load(last(fields));@assert all(isfinite,f["density"])&&minimum(f["density"])>0
end
function main()
    mode,out=ARGS[1:2];ispath(out)&&error("Refuse existing output");mkpath(out)
    if mode=="preflight"
        c=config();energy_budget_notice(.11,1.1)
        operator_checks(merge(c,Dict("initial_rms"=>0.,"check_output"=>out)))
        cc=merge(c,Dict("nx"=>64,"ny"=>64,"nz"=>80))
        m,p=build_model(cc;arch=CPU(),perturb=false)
        balance=merge(Dict{String,Any}(initial_balance(m,p)),p.seed_info)
        @assert balance["horizontal_imbalance_rms_over_fU"]<.15
        @assert balance["vertical_imbalance_max_over_fU"]<1e-8
        @assert abs(balance["coriolis_work_over_bottom_power"])<.01
        write_toml(joinpath(out,"INITIAL_BALANCE.toml"),balance)
        p3=parameters(merge(cc,Dict("omega"=>3.)));p45=parameters(cc)
        s3=seed_velocity(p3,cc);s45=seed_velocity(p45,cc)
        velocity_difference=maximum(maximum(abs,s3[i].-s45[i]) for i=2:3);@assert velocity_difference==0
        m=nothing;GC.gc()
        c=small_config()
        a,_=run_case(c,joinpath(out,"first");arch=CPU())
        cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
        b,_=run_case(merge(c,Dict("stop_time"=>.2)),joinpath(out,"restart");arch=CPU(),pickup=cp)
        e,_=run_case(merge(c,Dict("stop_time"=>.2)),joinpath(out,"direct");arch=CPU())
        fa=Oceananigans.prognostic_fields(b.model);fb=Oceananigans.prognostic_fields(e.model)
        err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-7
        for name in ["first","restart","direct"];readback(joinpath(out,name));end
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","omega_constant"=>4.5,"restart_max_abs_error"=>err,"initial_velocity_difference_from_omega3"=>velocity_difference,"balance"=>balance,"density_departure_gate"=>"Old amplitude heuristic reported, not enforced; positivity, stable cap, near-neutral bulk, discrete force balance and exact restart remain mandatory","scope"=>"Original mixed-sign velocity field at constant Omega=4.5; balanced thermodynamic initialization, full-vector Coriolis, actual driver output and restart"))
    elseif mode=="run"
        CUDA.functional()||error("GPU required")
        # Same executable and constant rotation in GPU smoke check and production.
        _,sm=run_case(small_config(),joinpath(out,"gpu_check"));@assert sm["status"]=="COMPLETE"
        readback(joinpath(out,"gpu_check"))
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"omega_constant"=>4.5))
        GC.gc();CUDA.reclaim()
        _,meta=run_case(config(),joinpath(out,"parent"))
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_continuation"=>false))
    else
        error("Unknown mode")
    end
end
include("model.jl")
include("initialization.jl")
include("checks.jl")
main()
