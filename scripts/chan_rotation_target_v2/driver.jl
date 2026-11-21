function main()
    mode,out=ARGS[1:2]
    ispath(out)&&error("Refuse overwrite");mkpath(out)
    cfg=TOML.parsefile("configs/chan-rotation-v2/case4.toml")
    if mode=="preflight"
        cfg["check_output"]=out;operator_checks(cfg)
        c=merge(cfg,Dict("nx"=>8,"ny"=>8,"nz"=>24,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"field_interval"=>.2,"wall_seconds"=>180.,"audit_only"=>false))
        a,_=run_case(c,joinpath(out,"heated");arch=CPU())
        cp=only(readdir(joinpath(out,"heated/checkpoints");join=true))
        b,_=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"restart");arch=CPU(),pickup=cp)
        e,_=run_case(merge(c,Dict("stop_time"=>.4)),joinpath(out,"direct");arch=CPU())
        aa=arrays(b.model);bb=arrays(e.model)
        err=maximum(maximum(abs,getproperty(aa,k).-getproperty(bb,k)) for k in (:uf,:vf,:wf,:rho,:T))
        @assert err<1e-7
        for name in ["heated","restart","direct"]
            rows=readlines(joinpath(out,name,"diagnostics.csv"));@assert length(rows)>1
            @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
            frames=filter(f->startswith(f,"frame_"),readdir(joinpath(out,name);join=false))
            @assert !isempty(frames)
            v=load(joinpath(out,name,last(frames)));@assert size(v["vorticity_slices"],3)==3
        end
        # Exercise second rotation and rest through identical output path.
        run_case(merge(c,Dict("omega"=>1.5,"heat"=>0.,"perturbation"=>0.)),joinpath(out,"rest");arch=CPU())
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"scope"=>"actual driver, operators, all output, native-face energy and restart"))
    elseif mode=="target"
        CUDA.functional()||error("GPU required")
        for (label,om) in [("case3",1.5),("case4",3.)]
            c=merge(cfg,Dict("omega"=>om,"stop_time"=>60.,"diagnostic_interval"=>.5,"movie_interval"=>5.,"checkpoint_interval"=>60.,"field_interval"=>60.,"wall_seconds"=>360.,"audit_only"=>true))
            _,mt=run_case(c,joinpath(out,label));@assert mt["status"]=="COMPLETE" "Target-grid check incomplete; no production release"
            GC.gc();CUDA.reclaim()
        end
        write_toml(joinpath(out,"TARGET_COMPLETE.toml"),Dict("status"=>"TARGET_COMPLETE_REVIEW_REQUIRED","production_authorized_by_receipt"=>false,"gpu"=>CUDA.name(CUDA.device()),"warning"=>"Final bounded target-grid audit. Not production. No physics or energy threshold changes."))
    else
        error("Unknown mode; no production route before energy review")
    end
end
include("../chan_rotation_v2/model.jl")
include("../chan_rotation_v2/checks.jl")
main()
