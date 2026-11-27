include("model.jl")
include("run.jl")
config()=TOML.parsefile(joinpath(@__DIR__,"config.toml"))
function small(c)
    τ=1/(c["rossby"]*c["f"])
    merge(c,Dict("nx"=>32,"ny"=>32,"nz"=>32,"turnovers"=>.04/τ,"initial_dt"=>.001,"diagnostic_turnovers"=>.01/τ,"movie_turnovers"=>.01/τ,"full_turnovers"=>.02/τ,"checkpoint_turnovers"=>.02/τ,"wall_seconds"=>300.,"test"=>true))
end
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    fs=filter(x->startswith(x,"frame_"),readdir(out));@assert !isempty(fs)
    for f in fs
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        for k in ("vorticity_slices","tracer_slices","vertical_velocity_slices","buoyancy_slices");@assert all(isfinite,a[k]);end
    end
    fs=filter(x->startswith(x,"full_"),readdir(out));@assert !isempty(fs)
    a=load(joinpath(out,last(sort(fs))));@assert maximum(abs,a["w_face"][:,:,[1,end]])==0
end
function exercise(c,ratio,out,arch)
    τ=1/(c["rossby"]*c["f"])
    a,ma=run_case(c,ratio,joinpath(out,"direct");arch)
    _,mf=run_case(merge(c,Dict("turnovers"=>.02/τ)),ratio,joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first","checkpoints");join=true))
    b,mb=run_case(c,ratio,joinpath(out,"restart");arch,pickup=cp)
    @assert ma["status"]==mf["status"]==mb["status"]=="COMPLETE"
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k]))-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-9
    # The cumulative cap heat is checkpoint state, not a reset segment budget.
    function last_cap_heat(part)
        lines=readlines(joinpath(out,part,"diagnostics.csv"))
        idx=findfirst(==("cumulative_cap_heat"),split(lines[1],','))
        parse(Float64,split(lines[end],',')[idx])
    end
    @assert isapprox(last_cap_heat("direct"),last_cap_heat("restart");atol=1e-13,rtol=1e-9)
    for part in ("direct","first","restart");readback(joinpath(out,part));end
    err
end
function initial_checks(c,out)
    result=Dict{String,Any}();common=nothing
    for (label,ratio) in (("weak",.5),("strong",2.))
        m,p,info=build(c,ratio;arch=CPU());a=fields_cpu(m)
        if common===nothing
            common=(copy(a.uf),copy(a.vf))
            info["velocity_match_relative_error"]=0.
        else
            # Native velocities are reconstructed from rho*u / rho, so allow
            # roundoff from the two slightly different reference densities.
            err=max(maximum(abs,a.uf-common[1]),maximum(abs,a.vf-common[2]))/p.seed.U
            info["velocity_match_relative_error"]=err
            println("VELOCITY_MATCH relative_error=$err");flush(stdout)
            @assert err<1e-12
        end
        result[label]=info
        # A short actual-model step quantifies imbalance using the projected velocity.
        sim=Simulation(m;Δt=.001,stop_time=.01);attach_cap_control!(sim,p);run!(sim);aa=fields_cpu(m)
        info["initial_adjustment_w_over_U"]=maximum(abs,aa.w)/p.seed.U
        @assert info["initial_adjustment_w_over_U"]<.01
        @assert constraint(m,p)<c["constraint_tolerance"]
        GC.gc()
    end
    write_toml(joinpath(out,"INITIAL_PAIR_PASS.toml"),result)
    result
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out);c=config()
    @assert c["nx"]==c["ny"]==256 && c["nz"]==80 && c["turnovers"]==10
    if mode=="preflight"
        # Full production initialization/balance plus complete small-grid driver paths.
        initial_checks(c,out)
        # Constant horizontal source has exactly zero anomalous component.
        sc=small(c);m,p,_=build(sc,.5;arch=CPU())
        interior(p.control.mean).=interior(p.control.target).+1e-4
        k=argmin(abs.(p.z.-.8))
        rate=mean_cap_source(1,1,k,m.grid,m.clock,nothing,p.control)
        @assert rate<0 && rate==mean_cap_source(2,3,k,m.grid,m.clock,nothing,p.control)
        @assert isapprox(rate/p.rho[k],-1e-4*Array(interior(p.control.mask))[1,1,k]/c["cap_relaxation_time"];rtol=1e-8)
        result=Dict{String,Any}("status"=>"PASS")
        for (label,ratio) in (("weak",.5),("strong",2.))
            sc=small(c);result[label*"_restart_error"]=exercise(sc,ratio,joinpath(out,label),CPU())
            m,p,_=build(sc,ratio;arch=CPU(),at_rest=true)
            sim=Simulation(m;Δt=.001,stop_time=.01);run!(sim);a=fields_cpu(m)
            err=max(maximum(abs,a.uf),maximum(abs,a.vf),maximum(abs,a.wf));@assert err<1e-9
            result[label*"_rest_speed"]=err
        end
        column_result=Dict{String,Any}()
        for (label,ratio) in (("weak",.5),("strong",2.))
            cc=merge(c,Dict("nx"=>8,"ny"=>8,"column"=>true,"turnovers"=>c["column_turnovers"],"wall_seconds"=>600.))
            sim,meta=run_case(cc,ratio,joinpath(out,"column_"*label);arch=CPU())
            @assert meta["status"]=="COMPLETE"
            a=fields_cpu(sim.model);n2=c["gravity"].*diff(log.(a.θ);dims=3)/(c["lz"]/c["nz"])
            z=[(k-.5)*c["lz"]/c["nz"] for k=1:c["nz"]-1];idx=findall(x->.7<x<.9,z)
            mincap=minimum(vec(mean(n2;dims=(1,2)))[idx]);N2=meta["initialization"]["N2"]
            err=maximum(abs,mean(a.θ;dims=(1,2))./reshape(meta["initialization"]["reference_theta"],1,1,:).-1)
            @assert mincap>.8N2
            @assert maximum(abs,a.w)<1e-8
            column_result[label]=Dict("status"=>"PASS","cap_min_N2_fraction"=>mincap/N2,"max_theta_mean_departure"=>err,"time"=>Float64(sim.model.clock.time),"maximum_w"=>maximum(abs,a.w),"mass_divergence_absolute"=>constraint_components(sim.model,(;dx=cc["lx"]/cc["nx"],dz=cc["lz"]/cc["nz"])).absolute)
            GC.gc()
        end
        write_toml(joinpath(out,"COLUMN_PASS.toml"),column_result)
        result["columns"]=column_result
        write_toml(joinpath(out,"DRIVER_PASS.toml"),result)
    elseif mode=="run"
        @assert CUDA.functional();ratio=parse(Float64,ARGS[3]);@assert ratio in (.5,2.)
        err=exercise(small(c),ratio,joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"gpu"=>CUDA.name(CUDA.device())))
        GC.gc();CUDA.reclaim()
        jobend=parse(Float64,get(ENV,"SLURM_JOB_END_TIME",string(time()+c["wall_seconds"]+300)))
        c["wall_seconds"]=min(c["wall_seconds"],jobend-time()-300)
        _,meta=run_case(c,ratio,joinpath(out,"parent");arch=GPU())
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"automatic_extension"=>false,"scope"=>"10-turnover convection/cap pilot; assess resolution, cap retention and cost before longer run"))
    else
        error("Unknown mode")
    end
end
main()
