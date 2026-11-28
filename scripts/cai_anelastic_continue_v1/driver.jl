include("model.jl")
include("initialization.jl")
include("run.jl")

function config()
    c=TOML.parsefile(joinpath(@__DIR__,"B7.toml"))
    @assert c["nx"]==c["ny"]==512 && c["nz"]==101
    @assert c["omega"]==2.56 && c["latitude"]==90 && c["initial_rms"]==c["beta"]==c["jet_speed"]==0
    @assert c["bottom_flux"]==.001 && c["stop_time"]==5000
    c
end
small(c)=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>16,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"field_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>600.))
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
    fs=filter(f->startswith(f,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        for key in ("vorticity_slices","tracer_slices","temperature_slices","vertical_velocity_slices");@assert all(isfinite,a[key]);end
    end
    a=load(joinpath(out,last(sort(filter(f->startswith(f,"full_"),readdir(out))))))
    @assert size(a["u_face"])==size(a["v_face"])==size(a["temperature"])
    @assert maximum(abs,a["w_face"][:,:,[1,end]])==0
end
function physics_checks(c,out)
    mkpath(out);p=parameters(c);T=[reference(zcell(k,p)).T for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    cond=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz];@assert maximum(abs,cond)<1e-12
    m,mp=build_science(merge(c,Dict("heat"=>0.,"perturbation"=>0.));arch=CPU())
    set!(m;u=(x,y,z)->.1sinpi(y/2),v=0.,w=0.);Oceananigans.TimeSteppers.update_state!(m)
    F=Oceananigans.fields(m);nu=Array(interior(mp.nu))
    expected=[c["deardorff_coefficient"]^2*sqrt(mp.dx*mp.dy)*mp.dz*sqrt(2ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,m.grid,F.u,F.v,F.w)) for i=1:mp.nx,j=1:mp.ny,k=1:mp.nz]
    err=maximum(abs,nu.-expected);@assert err<1e-12
    q=[heat_source(i,j,k,m.grid,m.clock,F,mp) for i=1:mp.nx,j=1:mp.ny,k=1:mp.nz]
    loss=-sum(Array(interior(m.velocities.u)).*Array(interior(m.timestepper.Gⁿ.ρu)))*mp.dx*mp.dy*mp.dz
    heat=sum(q)*mp.dx*mp.dy*mp.dz;@assert loss>0 && abs(heat/loss-1)<1e-8
    # Compare initialization against the actual compressible constructor, not analytic assumptions.
    am,ap=build_science(c;arch=CPU());cm,cp=build_model(c;arch=CPU())
    rhoerr=maximum(abs,parent(am.dynamics.reference_state.density).-parent(cm.dynamics.reference_state.density))
    perr=maximum(abs,parent(am.dynamics.reference_state.pressure).-parent(cm.dynamics.reference_state.pressure))
    thetaerr=maximum(abs,Array(interior(am.formulation.potential_temperature_density))./Array(interior(ap.rho0)).-arrays(cm).theta)
    @assert rhoerr==0 && perr==0 && thetaerr<1e-10
    checks=Dict{String,Any}("conductive_equilibrium_max_source"=>maximum(abs,cond),"viscosity_formula_error"=>err,"viscous_ke_loss"=>loss,"viscous_heat_return"=>heat,
        "reference_density_error"=>rhoerr,"reference_pressure_error"=>perr,"initial_theta_error"=>thetaerr)
    write_toml(joinpath(out,"PHYSICS_PASS.toml"),checks);checks
end
function resume_segment(c,out,request;arch)
    r=TOML.parsefile(request)
    @assert r["target_time"]==c["stop_time"] && r["restart_time"]<c["stop_time"]
    @assert isfile(r["checkpoint"])
    cfg=merge(c,Dict("restart_time"=>r["restart_time"],"restart_iteration"=>r["restart_iteration"],"restart_checkpoint"=>r["checkpoint"],
        "budget_scope"=>"Per segment since restart; do not interpret as cumulative since t=0",
        "profile_drift_scope"=>"Relative to segment-start profile; theta departure retains original fixed reference"))
    sim,meta=run_case(cfg,joinpath(out,"parent");arch,pickup=r["checkpoint"])
    @assert meta["final_time"]>r["restart_time"]
    checkpoints=filter(p->endswith(p,".jld2"),readdir(joinpath(out,"parent/checkpoints");join=true));@assert length(checkpoints)==1
    outcome=Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"final_iteration"=>sim.model.clock.iteration,
        "checkpoint"=>abspath(only(checkpoints)),"start_time"=>r["restart_time"],"target_time"=>c["stop_time"],"automatic_continuation"=>true)
    write_toml(joinpath(out,"OUTCOME.toml"),outcome)
    sim,meta
end
function continuation_check(c,out,arch)
    # Exercise the same request, resume, output and outcome path as production.
    c=merge(c,Dict("stop_time"=>.04))
    a,_=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"seed");arch)
    cp=only(readdir(joinpath(out,"seed/checkpoints");join=true))
    request=joinpath(out,"request.toml")
    write_toml(request,Dict("checkpoint"=>abspath(cp),"restart_time"=>.02,"restart_iteration"=>a.model.clock.iteration,"target_time"=>.04))
    b,_=resume_segment(c,joinpath(out,"resumed"),request;arch)
    direct,_=run_case(c,joinpath(out,"direct");arch)
    fa=Oceananigans.prognostic_fields(direct.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-7
    readback(joinpath(out,"resumed/parent"));readback(joinpath(out,"direct"))
    result=TOML.parsefile(joinpath(out,"resumed/OUTCOME.toml"));@assert result["status"]=="COMPLETE"
    err
end
function main()
    mode,out=ARGS[1:2];@assert mode in ("preflight","resume");@assert !ispath(out);mkpath(out);c=config()
    if mode=="preflight"
        result=physics_checks(small(c),joinpath(out,"physics"))
        result["restart_error"]=continuation_check(small(c),joinpath(out,"continuation"),CPU())
        warn=merge(small(c),Dict("perturbation"=>.8))
        _,meta=run_case(warn,joinpath(out,"warning_branch");arch=CPU());@assert meta["status"]=="COMPLETE"
        @assert isfile(joinpath(out,"warning_branch/APPROXIMATION_WARNING.toml"));readback(joinpath(out,"warning_branch"))
        result["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    else
        @assert CUDA.functional()
        err=continuation_check(small(c),joinpath(out,"gpu_validation"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"restart_error"=>err))
        GC.gc();CUDA.reclaim()
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+43200)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-900);@assert c["wall_seconds"]>60
        resume_segment(c,out,ARGS[3];arch=GPU())
    end
    println("CAI_ANELASTIC_CONTINUATION_DONE ",mode);flush(stdout)
end
main()
