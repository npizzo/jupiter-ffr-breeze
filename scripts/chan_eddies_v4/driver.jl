function config(case)
    c=TOML.parsefile("configs/chan-eddies-v4/"*case*".toml")
    c["stop_time"]=c["initial_turnovers"]*4/(sqrt(5)*c["initial_rms"])
    c
end
function main()
    mode,out=ARGS[1:2]
    ispath(out)&&error("Refuse existing output");mkpath(out)
    if mode=="preflight"
        energy_budget_notice(.11,1.1) # Exercise the warning-only policy through the actual driver.
        c=config("case4");operator_checks(merge(c,Dict("initial_rms"=>0.,"check_output"=>out)))
        reports=Dict{String,Any}();states=[]
        for label in ["case3","case4"]
            cc=merge(config(label),Dict("nx"=>64,"ny"=>64,"nz"=>80))
            m,p=build_model(cc;arch=CPU(),perturb=false)
            balance=merge(Dict{String,Any}(initial_balance(m,p)),p.seed_info)
            write_toml(joinpath(out,label*"_INITIAL_BALANCE.toml"),balance)
            @assert balance["horizontal_imbalance_rms_over_fU"]<.15 "Initial horizontal balance inadequate"
            @assert balance["vertical_imbalance_max_over_fU"]<1e-8 "Initial vertical balance inadequate"
            @assert abs(balance["coriolis_work_over_bottom_power"])<.01 "Coriolis injects appreciable energy"
            reports[label]=balance;push!(states,(arrays(m).uf,arrays(m).vf));m=nothing;GC.gc()
        end
        @assert maximum(maximum(abs,states[1][i].-states[2][i]) for i=1:2)<1e-13
        for key in ["mass","kinetic_energy"]
            @assert abs(reports["case3"][key]/reports["case4"][key]-1)<1e-10
        end
        write_toml(joinpath(out,"BALANCE_PASS.toml"),reports)
        c=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>24,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"field_interval"=>.2,"wall_seconds"=>180.))
        a,_=run_case(c,joinpath(out,"first");arch=CPU())
        cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
        b,_=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"restart");arch=CPU(),pickup=cp)
        e,_=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"direct");arch=CPU())
        aa=arrays(b.model);bb=arrays(e.model)
        err=maximum(maximum(abs,getproperty(aa,k).-getproperty(bb,k)) for k in (:uf,:vf,:wf,:rho,:T))
        @assert err<1e-7
        oldcp=only(readdir("data/chan-eddies-v3-cpu-2089424/first/checkpoints";join=true))
        cross,_=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"cross_version_restart");arch=CPU(),pickup=oldcp)
        ccross=arrays(cross.model)
        cross_error=maximum(maximum(abs,getproperty(ccross,k).-getproperty(bb,k)) for k in (:uf,:vf,:wf,:rho,:T))
        @assert cross_error<1e-7 "Diagnostic-only change altered restarted physics"

        run_case(merge(c,Dict("omega"=>1.5)),joinpath(out,"lower_rotation");arch=CPU())
        for name in ["first","restart","direct","lower_rotation","cross_version_restart"]
            rows=readlines(joinpath(out,name,"diagnostics.csv"));@assert length(rows)>1
            @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
            frames=filter(f->startswith(f,"frame_"),readdir(joinpath(out,name)))
            @assert !isempty(frames)
            v=load(joinpath(out,name,last(frames)));@assert size(v["vorticity_slices"],3)==3
        end
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"cross_version_restart_error"=>cross_error,"policy"=>"Warning-only energy residual; no changes to physical equations or parameters","scope"=>"Actual driver; matched velocity, mass and KE; stable initialized cap; balance; output/restart"))
    elseif mode=="screen"
        CUDA.functional()||error("GPU required")
        for label in ["case3","case4"]
            c=merge(config(label),Dict("stop_time"=>5.,"diagnostic_interval"=>.1,"movie_interval"=>.5,"checkpoint_interval"=>5.,"field_interval"=>5.,"wall_seconds"=>360.))
            _,mt=run_case(c,joinpath(out,label));@assert mt["status"]=="COMPLETE" "Screen partial; no production release"
            GC.gc();CUDA.reclaim()
        end
        write_toml(joinpath(out,"SCREEN_COMPLETE.toml"),Dict("status"=>"COMPLETE_REVIEW_REQUIRED","gpu"=>CUDA.name(CUDA.device())))
    elseif mode in ["case3","case4"]
        CUDA.functional()||error("GPU required")
        gate=TOML.parsefile(ARGS[3]);@assert gate["status"]=="PASS"
        c=config(mode)
        cp=gate[mode*"_checkpoint"];@assert bytes2hex(open(sha256,cp))==gate[mode*"_checkpoint_sha256"]
        _,meta=run_case(c,joinpath(out,"parent");pickup=cp)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_continuation"=>false))
    else
        error("Unknown mode")
    end
end
include("model.jl")
include("initialization.jl")
include("checks.jl")
main()
