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
    elseif mode=="energy"
        CUDA.functional()||error("GPU required")
        cases=[("omega1p5",1.5,32,80,.02),("omega3",3.,32,80,.02),("omega1p5_dt2",1.5,32,80,.008333333333333333),("omega3_dt2",3.,32,80,.008333333333333333),("omega1p5_grid2",1.5,64,80,.02),("omega3_grid2",3.,64,80,.02)]
        for (label,om,n,nz,dt) in cases
            c=merge(cfg,Dict("omega"=>om,"nx"=>n,"ny"=>n,"nz"=>nz,"max_dt"=>dt,"stop_time"=>40.,"diagnostic_interval"=>.1,"movie_interval"=>5.,"checkpoint_interval"=>40.,"field_interval"=>40.,"wall_seconds"=>150.,"audit_only"=>true))
            _,mt=run_case(c,joinpath(out,label));@assert mt["status"]=="COMPLETE" "Audit incomplete; no production release"
            GC.gc();CUDA.reclaim()
        end
        write_toml(joinpath(out,"AUDIT_COMPLETE.toml"),Dict("status"=>"AUDIT_COMPLETE_REVIEW_REQUIRED","production_authorized_by_receipt"=>false,"gpu"=>CUDA.name(CUDA.device()),"warning"=>"Energy guard disabled ONLY for bounded t40 audit to locate discrepancy. This is not a production pass."))
    else
        error("Unknown mode; no production route before energy review")
    end
end
include("model.jl")
include("checks.jl")
main()
