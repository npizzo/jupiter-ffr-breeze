function config()
    c=TOML.parsefile("configs/omega3-anelastic-long-v1.toml");@assert c["omega"]==3.;c
end
function small_config()
    merge(config(),Dict("nx"=>16,"ny"=>16,"nz"=>24,"stop_time"=>.1,"diagnostic_interval"=>.025,"movie_interval"=>.025,"field_interval"=>.05,"checkpoint_interval"=>.05,"wall_seconds"=>240.,"max_dt"=>.0025))
end
function difference(a,b)
    fa=Oceananigans.prognostic_fields(a);fb=Oceananigans.prognostic_fields(b)
    @assert keys(fa)==keys(fb)
    maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
end
function readback(out)
    rows=[parse.(Float64,split(l,',')) for l in readlines(joinpath(out,"diagnostics.csv"))[2:end]]
    @assert length(rows)>1 && all(all(isfinite,r) for r in rows)
    fs=filter(f->startswith(f,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        @assert all(isfinite,a["vorticity_slices"])&&all(isfinite,a["tracer_slices"])
    end
    fs=readdir(joinpath(out,"fields3d");join=true);@assert !isempty(fs)
    a=load(last(sort(fs)));@assert all(isfinite,a["theta"])&&minimum(a["temperature"])>0
end
function restart_check(c,out,arch)
    split=merge(c,Dict("stop_time"=>.05));first,_=run_case(split,joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    resumed,_=run_case(merge(c,Dict("restart_time"=>.05)),joinpath(out,"restart");arch,pickup=cp)
    direct,_=run_case(c,joinpath(out,"direct");arch)
    for f in ["first","restart","direct"];readback(joinpath(out,f));end
    err=difference(resumed.model,direct.model);@assert err<1e-7
    err
end
function checks(c,out)
    m,p=build_model(c;arch=CPU())
    # At a common thermodynamic reference, the entropy flux is precisely the
    # original dry constitutive law; at nonzero perturbation use EOS density.
    T=[reference(zcell(k,p)).T*(1+.01cospi(2*(i-.5)/p.nx)) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    P=repeat(Array(interior(p.pref)),p.nx,p.ny,1);rho=P./T
    set!(p.nu,.001);Oceananigans.BoundaryConditions.fill_halo_regions!(p.nu)
    F=(T=T,ρᵈ=rho);cp=CompressibleControl.parameters(c);cp=merge(cp,(;nu=p.nu))
    qp=merge(p,(;rho));qs=[sgs_heating(i,j,k,F,qp) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    old=[CompressibleControl.sgs_heating(i,j,k,F,cp) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    @assert maximum(abs,qs.-old)<1e-10 && abs(mean(qs))<1e-11
    Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    qt=mean(2p.kc.*(T[:,:,end].-1)./p.dz)
    @assert abs(mean(Q)-(p.flux-qt))<1e-12
    # Full production heat source, with zero velocity and SGS viscosity, maps
    # into rho0*theta tendency as Q theta/(cp T), using Breeze's native kernel.
    set!(m;u=0.,v=0.,w=0.)
    Oceananigans.TimeSteppers.update_state!(m)
    Fm=Oceananigans.fields(m)
    G=Array(interior(m.timestepper.Gⁿ.ρθ))
    expected=[heat_source(i,j,k,m.grid,m.clock,Fm,p)*Fm.ρθ[i,j,k]/p.rho[i,j,k]/(2.5Fm.T[i,j,k]) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    herr=maximum(abs,G.-expected);@assert herr<1e-8 "Native heating tendency mismatch"
    # Hydrostatic reference must remain at rest with all external heat switched off.
    restcfg=merge(c,Dict("heat"=>0.));rest,rp=build_model(restcfg;arch=CPU(),rest=true)
    sim=Simulation(rest;Δt=.001,stop_time=.01);run!(sim)
    ra=arrays(rest);resterr=max(maximum(abs,ra.uf),maximum(abs,ra.vf),maximum(abs,ra.wf))
    @assert resterr<1e-9
    write_toml(joinpath(out,"OPERATOR_PASS.toml"),Dict("status"=>"PASS","entropy_constitutive_max_error"=>maximum(abs,qs.-old),"native_heating_tendency_max_error"=>herr,"rest_speed_max"=>resterr))
end
function main()
    mode,out=ARGS[1:2];ispath(out)&&error("Refuse overwrite");mkpath(out)
    if mode=="preflight"
        checks(small_config(),out)
        err=restart_check(small_config(),out,CPU())
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_max_abs_error"=>err,"scope"=>"Actual driver, native heating, entropy constitutive equivalence, rest equilibrium, initial mapping, diagnostic/movie/full-field readback, serialized restart"))
    elseif mode=="run"
        CUDA.functional()||error("GPU required")
        err=restart_check(small_config(),joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_max_abs_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        GC.gc();CUDA.reclaim()
        origin=TOML.parsefile("runs/omega3-anelastic-long-v1/ORIGIN.toml")
        @assert origin["status"]=="PASS" && origin["restart_time"]==5.
        @assert bytes2hex(open(sha256,origin["review"]))==origin["review_sha256"]
        cp=origin["checkpoint"];@assert bytes2hex(open(sha256,cp))==origin["checkpoint_sha256"]
        _,meta=run_case(merge(config(),Dict("restart_time"=>5.)),joinpath(out,"parent");pickup=cp)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_long_release"=>false))
    else
        error("Unknown mode")
    end
end
module CompressibleControl
include("../omega3_long_v1/model.jl")
include("../omega3_long_v1/initialization.jl")
end
include("model.jl")
main()
