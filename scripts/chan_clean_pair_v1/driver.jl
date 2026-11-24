include("model.jl")
include("run.jl")
function config(eq)
    c=merge(TOML.parsefile(joinpath(@__DIR__,"config.toml")),Dict("equations"=>eq))
    @assert c["omega"]==3 && c["beta"]==0 && c["jet_speed"]==0 && c["initial_rms"]==0
    @assert c["closure"]=="none" && c["advection_order"]==9
    c
end
small(c)=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>24,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"field_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>300.))
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
function operators(c)
    m,p=build_model(c;arch=CPU());@assert m.closure===nothing
    @assert Oceananigans.Grids.topology(m.grid)==(Periodic,Periodic,Bounded)
    F=Oceananigans.fields(m);a=arrays(m);dv=p.dx*p.dy*p.dz
    Q=[conduction(i,j,k,a.T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    qt=mean(2p.kc.*(a.T[:,:,end].-1)./p.dz)
    fluxerr=abs(sum(Q)*dv-16*(p.flux-qt));@assert fluxerr<1e-10
    Oceananigans.TimeSteppers.update_state!(m);F=Oceananigans.fields(m)
    expected=[heat_source(i,j,k,m.grid,m.clock,F,p)*F.ρθ[i,j,k]/mass_density(i,j,k,F,p)/(2.5F.T[i,j,k]) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    herr=maximum(abs,Array(interior(m.timestepper.Gⁿ.ρθ)).-expected);@assert herr<1e-8
    totalQ=sum(heat_source(i,j,k,m.grid,m.clock,F,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz)*dv
    cooling=dv*sum(a.rho.*(a.T.-1).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
    powererr=abs(totalQ-p.heat*(16*(p.flux-qt)-cooling));@assert powererr<1e-10
    restcfg=merge(c,Dict("heat"=>0.,"perturbation"=>0.));r,rp=build_model(restcfg;arch=CPU())
    sim=Simulation(r;Δt=.001,stop_time=.01);run!(sim);ra=arrays(r)
    resterr=max(maximum(abs,ra.uf),maximum(abs,ra.vf),maximum(abs,ra.wf));@assert resterr<1e-8
    # Audit Coriolis orientation directly, without pressure-projection contamination.
    set!(r;u=1.,v=2.,w=0.);M=Oceananigans.prognostic_fields(r)
    k=div(rp.nz,2);rho=Array(interior(density_field(r)))[1,1,k]
    xx=Oceananigans.Coriolis.x_f_cross_U(5,5,k,r.grid,r.coriolis,M)
    yy=Oceananigans.Coriolis.y_f_cross_U(5,5,k,r.grid,r.coriolis,M)
    corerr=max(abs(xx+2rp.f*rho),abs(yy-rp.f*rho));@assert corerr<1e-8
    Dict("heat_flux_telescoping_error"=>fluxerr,"native_theta_source_error"=>herr,"local_cooling_power_error"=>powererr,"rest_speed"=>resterr,"coriolis_horizontal_error"=>corerr)
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
    if mode=="preflight"
        result=Dict{String,Any}("status"=>"PASS")
        for eq in ["compressible","anelastic"]
            c=small(config(eq));op=operators(c);op["restart_error"]=exercise(c,joinpath(out,eq),CPU());result[eq]=op;GC.gc()
        end
        write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    elseif mode=="run"
        CUDA.functional()||error("GPU required")
        eq=["compressible","anelastic"][parse(Int,ARGS[3])+1];c=config(eq)
        err=exercise(small(c),joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        GC.gc();CUDA.reclaim()
        # Hard 12 GPU-hour allocation includes compilation and validation overhead.
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+c["wall_seconds"]+300)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-300);@assert c["wall_seconds"]>60
        _,meta=run_case(c,joinpath(out,"parent");arch=GPU())
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_extension"=>false))
    else
        error("Unknown mode")
    end
end
main()
