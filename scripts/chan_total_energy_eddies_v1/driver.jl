include("model.jl")
include("initialization.jl")
include("run.jl")
function config(eq)
    c=merge(TOML.parsefile(joinpath(@__DIR__,"config.toml")),Dict("equations"=>eq))
    @assert c["omega"]==3 && c["beta"]==0 && c["jet_speed"]==0 && c["initial_rms"]==.136
    @assert c["closure"]=="none" && c["advection_order"]==9
    c
end
small(c)=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>24,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"field_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>300.,"budget_gate_start"=>.01))
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
    fs=filter(f->startswith(f,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        @assert all(isfinite,a["vorticity_slices"]) && all(isfinite,a["tracer_slices"])
    end
    fs=filter(f->startswith(f,"full_"),readdir(out));@assert !isempty(fs)
    a=load(joinpath(out,last(sort(fs))))
    @assert size(a["u_face"])==size(a["v_face"])==size(a["temperature"])
    @assert maximum(abs,a["w_face"][:,:,[1,end]])==0 && all(isfinite,a["theta"])
end
function exercise(c,out,arch)
    _,meta=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch)
    @assert meta["status"]=="COMPLETE"
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,bmeta=run_case(merge(c,Dict("restart_time"=>.02)),joinpath(out,"restart");arch,pickup=cp)
    a,ameta=run_case(c,joinpath(out,"direct");arch)
    @assert ameta["status"]==bmeta["status"]=="COMPLETE"
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
    @assert err<1e-7
    for part in ("first","restart","direct");readback(joinpath(out,part));end
    err
end
function energy_checks(c,out,arch)
    # Identical physical initialization on the same revised solver.
    ct=merge(c,Dict("formulation"=>"PotentialTemperature"))
    ce=merge(c,Dict("formulation"=>"TotalEnergy"))
    mt,pt=build_model(ct;arch);me,pe=build_model(ce;arch)
    at=arrays(mt);ae=arrays(me)
    initial_error=maximum(maximum(abs,getproperty(at,k).-getproperty(ae,k)) for k in (:rho,:T,:uf,:vf,:wf))
    @assert initial_error<1e-10
    mt=nothing;me=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    me,pe=build_model(merge(ce,Dict("initial_rms"=>0.));arch)
    # At rest, forcing is the entire energy tendency (all advective energy fluxes vanish).
    Oceananigans.TimeSteppers.update_state!(me)
    F=Oceananigans.fields(me)
    a=arrays(me)
    expected=[conduction(i,j,k,a.T,pe)-a.rho[i,j,k]*(a.T[i,j,k]-1)*cool_rate(k,pe) for i=1:pe.nx,j=1:pe.ny,k=1:pe.nz]
    source_error=maximum(abs,Array(interior(me.timestepper.Gⁿ.ρE)).-expected)
    @assert source_error<1e-9
    mt=nothing;me=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    function ledger(path)
        rows=readlines(joinpath(path,"diagnostics.csv"));names=split(rows[1],',')
        v=[parse.(Float64,split(r,',')) for r in rows[2:end]]
        return Dict(k=>[r[i] for r in v] for (i,k) in enumerate(names))
    end
    checks=Dict{String,Any}("initial_physical_state_max_error"=>initial_error,"native_energy_source_max_error"=>source_error)
    for name in ("closed_moving","constant_heat","original_heat")
        q=merge(ce,Dict("heat"=>name=="closed_moving" ? 0. : 1., "heat_mode"=>name=="constant_heat" ? "constant" : "original", "smooth_velocity"=>false))
        sim,meta=run_case(q,joinpath(out,name);arch);@assert meta["status"]=="COMPLETE"
        readback(joinpath(out,name));d=ledger(joinpath(out,name));E=d["energy_proxy"]
        dt=last(d["time"])-first(d["time"])
        expected=name=="constant_heat" ? 16*.001*dt : name=="closed_moving" ? 0. : last(d["cumulative_heat"])
        err=abs(last(E)-first(E)-expected)
        if name=="original_heat"
            @assert err<max(1e-9,1e-3*.5*dt)
        else
            @assert err/max(abs(first(E)),1.)<1e-11
        end
        @assert maximum(abs,d["reconstruction_minus_native"])/max(abs(first(E)),1.)<1e-11
        checks[name]=Dict("absolute_energy_residual"=>err,"relative_to_initial_energy"=>err/max(abs(first(E)),1.),"expected_energy_input"=>expected,"max_mach"=>maximum(d["max_mach"]))
        sim=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    end
    checks
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
    @assert mode in ("preflight","gpu_validate","run")
    arch=mode=="preflight" ? CPU() : GPU()
    mode!="preflight" && @assert CUDA.functional()
    if mode=="run"
        c=config("compressible");@assert c["formulation"]=="TotalEnergy" && c["stop_time"]==2000.
        seed_checks=energy_checks(small(c),joinpath(out,"gpu_energy_check"),GPU())
        err=exercise(small(c),joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device()),"energy_and_seed_checks"=>seed_checks))
        GC.gc();CUDA.reclaim()
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+c["wall_seconds"]+300)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-300);@assert c["wall_seconds"]>60
        _,meta=run_case(c,joinpath(out,"parent");arch=GPU())
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_extension"=>false))
        return
    end
    result=energy_checks(small(config("compressible")),joinpath(out,"energy"),arch)
    for form in ("PotentialTemperature","TotalEnergy")
        c=merge(small(config("compressible")),Dict("formulation"=>form))
        result[form]=Dict("restart_error"=>exercise(c,joinpath(out,form),arch))
        GC.gc();arch isa GPU && CUDA.reclaim()
    end
    result["status"]="PASS"
    result["revision"]="8273be89d9357ca2ac572347c45db8c647e3c45c"
    write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    println("TOTAL_ENERGY_DRIVER_PASS");flush(stdout)
end
main()
