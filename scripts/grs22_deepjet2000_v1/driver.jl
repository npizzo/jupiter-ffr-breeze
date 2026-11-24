include("model.jl")
include("initialization.jl")
include("cap_initialization.jl")
include("eddy_and_source.jl")
include("guards.jl")
include("config.jl")
include("run.jl")
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    spatial=readlines(joinpath(out,"mach_spatial.csv"));@assert length(spatial)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in spatial[2:end])
    for pre in ("frame_","full_")
        fs=sort(filter(f->startswith(f,pre),readdir(out)));@assert !isempty(fs)
        x=load(joinpath(out,last(fs)))
        if pre=="frame_";@assert size(x["vorticity_slices"],3)==3 && all(isfinite,x["tracer_slices"]) && all(isfinite,x["ubar"]);end
        if pre=="full_";@assert maximum(abs,x["v_face"][:,[1,end],:])==0 && maximum(abs,x["w_face"][:,:,[1,end]])==0;end
    end
end
function small(c)
    merge(c,Dict("nx"=>16,"ny"=>32,"nz"=>40,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"field_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>240.))
end
function exercise(c,out,arch)
    _,meta=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch);@assert meta["status"]=="COMPLETE"
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    rc=merge(c,Dict("restart_time"=>.02,"health_reference_path"=>joinpath(out,"first/health_reference.jld2"),"budget_reference_path"=>joinpath(out,"first/budget_reference.jld2")))
    b,bmeta=run_case(rc,joinpath(out,"restart");arch,pickup=cp)
    a,ameta=run_case(c,joinpath(out,"direct");arch)
    @assert bmeta["status"]==ameta["status"]=="COMPLETE"
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-7
    for part in ("first","restart","direct");readback(joinpath(out,part));end
    err
end
function operators(c)
    m,p=build_model(c;arch=CPU());@assert m.closure===nothing
    set!(m;u=0.,v=0.,w=0.);Oceananigans.TimeSteppers.update_state!(m)
    F=Oceananigans.fields(m);T=Array(interior(m.temperature));dv=p.dx*p.dy*p.dz
    Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    qtop=p.top_factor*mean(2p.kc.*(T[:,:,end].-p.tt)./p.dz)
    fluxerr=abs(sum(Q)*dv-p.lx*p.ly*(p.bottom_factor*p.flux-qtop));@assert fluxerr<1e-10
    expected=[heat_source(i,j,k,m.grid,m.clock,F,p)*F.ρθ[i,j,k]/mass_density(i,j,k,F,p)/(2.5F.T[i,j,k]) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    herr=maximum(abs,Array(interior(m.timestepper.Gⁿ.ρθ)).-expected);@assert herr<1e-8
    # Resting reference without external sources tests hydrostatic support and wall geometry.
    restcfg=merge(c,Dict("initial_rms"=>0.,"jet_speed"=>0.,"bottom_factor"=>0.,"top_factor"=>0.,"conduction_factor"=>0.,"cap_factor"=>0.))
    r,rp=build_model(restcfg;arch=CPU());sim=Simulation(r;Δt=.001,stop_time=.01);run!(sim)
    ra=arrays(r);resterr=max(maximum(abs,ra.uf),maximum(abs,ra.vf),maximum(abs,ra.wf));@assert resterr<1e-8
    Dict("heat_flux_telescoping_error"=>fluxerr,"native_theta_source_error"=>herr,"rest_speed"=>resterr,"closure"=>"nothing")
end
function science_checks(c)
    p=parameters(c)
    @assert p.f<0 && p.fh>0 && p.beta>0 && !p.cap_jet
    @assert isapprox(p.beta/abs(p.f),250/71492*cotd(22.);rtol=1e-12)
    @assert jet_target(-p.jet_width,.2,p)>0 && jet_target(p.jet_width,.2,p)<0
    @assert jet_target(p.jet_width,.2,p)==jet_target(p.jet_width,.9,p)
    # Zeta=-dU/dy is positive at the center, hence anticyclonic where f<0.
    shear=(jet_profile(1e-5,p)-jet_profile(-1e-5,p))/2e-5
    @assert -shear/p.f<0
    # An x-uniform mean jet may have zero material acceleration, but tilted
    # Coriolis requires both horizontal pressure and vertical density balance.
    cjet=merge(c,Dict("initial_rms"=>0.))
    m,mp=build_model(cjet;arch=CPU());r=spatial_screen(m,mp)
    @assert r["cap_N2_min"]>0
    @assert r["vertical_balance_max"]<1e-8
    @assert r["horizontal_balance_rms_over_fU"]<.15
    merge(r,Dict("f"=>p.f,"fh"=>p.fh,"beta"=>p.beta,"jet_peak"=>abs(p.jet_speed),"central_shear"=>shear,"jet_to_beta_width2"=>abs(p.jet_speed)/(p.beta*p.jet_width^2)))
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
    c=base_config("compressible",true)
    if mode=="preflight"
        sc=small(c);guard_checks(sc);op=operators(sc)
        op["restart_error"]=exercise(sc,joinpath(out,"compressible_deepjet"),CPU())
        op["science_checks"]=science_checks(sc);op["status"]="PASS"
        write_toml(joinpath(out,"DRIVER_PASS.toml"),op)
    elseif mode=="run"
        CUDA.functional() || error("GPU required")
        err=exercise(small(c),joinpath(out,"gpu_check"),GPU());GC.gc();CUDA.reclaim()
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        pc=merge(c,Dict("stop_time"=>10.,"wall_seconds"=>1800.,"diagnostic_interval"=>.5,"movie_interval"=>1.,"field_interval"=>5.,"checkpoint_interval"=>5.))
        _,pm=run_case(pc,joinpath(out,"pilot");arch=GPU())
        @assert pm["status"]=="COMPLETE" && pm["final_time"]==10.
        readback(joinpath(out,"pilot"))
        lines=readlines(joinpath(out,"pilot","diagnostics.csv"));lastrow=Dict(zip(split(lines[1],','),parse.(Float64,split(lines[end],','))))
        @assert lastrow["cap_mean_N2_min"]>0 && lastrow["mach"]<.8 && lastrow["mass_error"]<1e-7
        write_toml(joinpath(out,"PILOT_PASS.toml"),Dict("status"=>"PASS","time"=>10.,"diagnostics"=>lastrow))
        cp=only(readdir(joinpath(out,"pilot/checkpoints");join=true))
        rc=merge(c,Dict("restart_time"=>10.,"health_reference_path"=>joinpath(out,"pilot/health_reference.jld2"),"budget_reference_path"=>joinpath(out,"pilot/budget_reference.jld2")))
        GC.gc();CUDA.reclaim()
        _,meta=run_case(rc,joinpath(out,"parent");arch=GPU(),pickup=cp)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_further_extension"=>false))
    else
        error("Unknown mode")
    end
end
main()
