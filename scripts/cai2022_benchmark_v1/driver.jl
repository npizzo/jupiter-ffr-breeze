include("model.jl")
include("run.jl")
function config(case)
    @assert case in ("B0","B7")
    c=TOML.parsefile(joinpath(@__DIR__,case*".toml"))
    @assert c["nx"]==c["ny"]==512 && c["nz"]==101
    @assert c["bottom_flux"]==.001 && c["latitude"]==90 && c["beta"]==c["jet_speed"]==c["initial_rms"]==0
    c
end
small(c)=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>16,"stop_time"=>.04,"initial_dt"=>.001,"max_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"field_interval"=>.02,"checkpoint_interval"=>.02,"wall_seconds"=>600.,"budget_gate_start"=>.01))
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
    fs=filter(f->startswith(f,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        for key in ("vorticity_slices","tracer_slices","temperature_slices","vertical_velocity_slices")
            @assert all(isfinite,a[key])
        end
    end
    a=load(joinpath(out,last(sort(filter(f->startswith(f,"full_"),readdir(out))))))
    @assert size(a["u_face"])==size(a["v_face"])==size(a["temperature"])
    @assert maximum(abs,a["w_face"][:,:,[1,end]])==0
end
function exercise(c,out,arch)
    _,meta=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,bmeta=run_case(merge(c,Dict("restart_time"=>.02)),joinpath(out,"restart");arch,pickup=cp)
    a,ameta=run_case(c,joinpath(out,"direct");arch)
    @assert ameta["status"]==bmeta["status"]==meta["status"]=="COMPLETE"
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
    @assert err<1e-7
    for part in ("first","restart","direct");readback(joinpath(out,part));end
    err
end
function physics_checks(c,out,arch)
    checks=Dict{String,Any}()
    # Analytic conductive equilibrium: T=5-4z, rho=T, p=T², g=8.
    p=parameters(c);T=[reference(zcell(k,p)).T for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    cond=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    @assert maximum(abs,cond)<1e-12
    @assert isapprox(p.kc*(GRAV/2.5)/p.flux,.8;atol=1e-14)
    checks["conductive_equilibrium_max_source"]=maximum(abs,cond)
    # CPU manufacturer: verify constitutive stress and work against finite differences of native tendencies.
    if arch isa CPU
        m,mp=build_model(merge(c,Dict("perturbation"=>0.,"heat"=>0.));arch)
        set!(m;u=(x,y,z)->.1sinpi(y/2),v=0.,w=0.)
        Oceananigans.TimeSteppers.update_state!(m);F=Oceananigans.fields(m)
        nu=Array(interior(m.closure_fields.νₑ))
        expect=[c["deardorff_coefficient"]^2*sqrt(mp.dx*mp.dy)*mp.dz*sqrt(2ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,m.grid,F.u,F.v,F.w)) for i=1:mp.nx,j=1:mp.ny,k=1:mp.nz]
        @assert maximum(abs,nu.-expect)<1e-12
        sources=[stress_work(i,j,k,m.grid,m.clock,F,mp) for i=1:mp.nx,j=1:mp.ny,k=1:mp.nz]
        @assert abs(sum(sources)*mp.dx*mp.dy*mp.dz)<1e-12
        # For u(y), momentum advection and x pressure force vanish: native tendency is viscosity alone.
        gu=Array(interior(m.timestepper.Gⁿ.ρu));uf=Array(interior(m.velocities.u))
        diss=-sum(uf.*gu)*mp.dx*mp.dy*mp.dz
        @assert diss>0
        # Validate actual native symmetric stress + trace correction for divergent u(x).
        set!(m;u=(x,y,z)->.1sinpi(x/2),v=0.,w=0.)
        Oceananigans.TimeSteppers.update_state!(m);F=Oceananigans.fields(m)
        # Compare additional trace force to -1/3 of symmetric longitudinal force.
        err=maximum(abs(fx(i,j,k,m.grid,m.clock,F,mp)+(2/3)*(F.ρᵈ[i,j,k]*mp.nu[i,j,k]*(F.u[i+1,j,k]-F.u[i,j,k])/mp.dx-F.ρᵈ[i-1,j,k]*mp.nu[i-1,j,k]*(F.u[i,j,k]-F.u[i-1,j,k])/mp.dx)/mp.dx) for i=1:mp.nx,j=1:mp.ny,k=1:mp.nz)
        @assert err<1e-12
        checks["cai_viscosity_formula_error"]=maximum(abs,nu.-expect)
        checks["closed_wall_stress_work_integral"]=sum(sources)*mp.dx*mp.dy*mp.dz
        checks["shear_kinetic_energy_dissipation_rate"]=diss
        checks["compressive_trace_correction_error"]=err
        m=nothing;GC.gc()
    end
    # Moving, adiabatically closed integration: energy retained while viscosity acts.
    q=merge(c,Dict("heat"=>0.,"smooth_velocity"=>true))
    sim,meta=run_case(q,joinpath(out,"closed_moving");arch)
    @assert meta["status"]=="COMPLETE";readback(joinpath(out,"closed_moving"))
    lines=readlines(joinpath(out,"closed_moving/diagnostics.csv"));names=split(lines[1],',');idx=findfirst(==("energy_proxy"),names)
    E=[parse(Float64,split(l,',')[idx]) for l in lines[2:end]]
    err=maximum(abs,E.-E[1])/abs(E[1]);@assert err<1e-10
    checks["closed_total_energy_relative_drift"]=err
    sim=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    checks
end
function main()
    mode,case,out=ARGS[1:3];@assert mode in ("preflight","gpu_validate","run")
    @assert !ispath(out);mkpath(out);c=config(case);arch=mode=="preflight" ? CPU() : GPU()
    mode!="preflight" && @assert CUDA.functional()
    if mode=="preflight"
        result=physics_checks(small(c),joinpath(out,"physics"),arch)
        for name in ("B0","B7")
            result[name]=Dict("restart_error"=>exercise(small(config(name)),joinpath(out,name),arch))
        end
        result["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    elseif mode=="gpu_validate"
        result=physics_checks(small(c),joinpath(out,"physics"),arch)
        result["restart_error"]=exercise(small(c),joinpath(out,"restart"),arch)
        GC.gc();CUDA.reclaim()
        # Production-size heated pilot tests memory, acoustic substeps, closures, I/O and real throughput.
        pilot=merge(c,Dict("stop_time"=>2.,"diagnostic_interval"=>.5,"movie_interval"=>1.,"field_interval"=>2.,"checkpoint_interval"=>2.,"wall_seconds"=>900.,"budget_gate_start"=>.5))
        sim,meta=run_case(pilot,joinpath(out,"target_grid");arch)
        @assert meta["status"]=="COMPLETE";readback(joinpath(out,"target_grid"))
        result["status"]="PASS";result["gpu"]=CUDA.name(CUDA.device())
        result["target_grid_time"]=meta["final_time"];result["target_grid_wall_seconds"]=meta["wall_seconds_actual"]
        result["note"]="Startup timing only; convection may be slower. No automatic production release."
        write_toml(joinpath(out,"GPU_PASS.toml"),result)
    else
        # Production requires exact GPU receipt and explicit cost acceptance by release.py.
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+c["wall_seconds"]+300)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-300)
        @assert c["wall_seconds"]>60
        _,meta=run_case(c,joinpath(out,"parent");arch)
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_extension"=>false))
    end
    println("CAI_DRIVER_DONE ",mode," ",case);flush(stdout)
end
main()
