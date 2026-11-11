include("run.jl")
include("check_ferrel.jl")
include("common.jl")
length(ARGS)==2 || error("usage: preflight.jl ROOT quick|512|1024")
root,phase=ARGS
phase in ("quick","512","1024") || error("Unknown preflight phase")
verify_gate(joinpath(root,"UNIT_PASS.toml"))
verify_gate(joinpath(root,"cpu_smoke","CPU_SMOKE_PASS.toml"))
phase!="quick" && verify_gate(joinpath(root,"QUICK_PASS.toml"))
phase=="1024" && verify_gate(joinpath(root,"PASS512.toml"))

function experiment(name,case,changes;pickup=false)
    cfg=merge(TOML.parsefile(config_path(case)),Dict(changes))
    config=joinpath(root,name*".toml")
    isfile(config) && error("Refusing repeated preflight: $name")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_ferrel_output(out)
    GC.gc();CUDA.reclaim()
    return out
end

if phase=="quick"
    small=Dict("nx"=>96,"ny"=>96,"nz"=>48,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,
        "controller_seconds"=>60.0,"checkpoint_seconds"=>600.0,"turbulence_wavelength"=>3200000.0,
        "native_budget_interval"=>1,"native_budget_first_steps"=>8)
    full=experiment("restart_full","b162",merge(small,Dict("stop_seconds"=>1200.0)))
    part=experiment("restart_first","b162",merge(small,Dict("stop_seconds"=>600.0)))
    resumed=experiment("restart_second","b162",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
    @testset "Independent restart continuation with native budgets" begin
        a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","native_rhou","native_rhov","native_rhow","native_rhou_rhs")
            @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
        end
        metadata=TOML.parsefile(joinpath(full,"metadata.toml"))
        derived=check_derived_streamfunction_restart(a,b,metadata)
        open(joinpath(root,"gpu_restart_streamfunction_check.toml"),"w") do io;TOML.print(io,derived);end
    end
    for case in ("axis0","axis162")
        experiment("quick_$(case)",case,merge(small,Dict("stop_seconds"=>1200.0)))
    end
    write_gate(joinpath(root,"QUICK_PASS.toml");grid=[96,96,48],simulated_seconds=1200,
        note="Diagnostic/restart gate; 3200-km coarse seed differs from 1200-km physics runs")
else
    n=parse(Int,phase);initial=nothing;estimates=Float64[]
    cases=n==512 ? ("b0","b162","axis0","axis162") : ("b162",)
    for case in cases
        out=experiment("pilot$(n)_$(case)",case,Dict("nx"=>n,"ny"=>n,"nz"=>96,"stop_seconds"=>21600.0,
            "output_seconds"=>3600.0,"diagnostic_seconds"=>600.0,"checkpoint_seconds"=>10800.0))
        if case in ("b0","b162")
            current=load(joinpath(out,"frame_00000.jld2"))
            if initial!==nothing
                for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","native_rhou","native_rhov","native_rhow")
                    @test isapprox(current[key],initial[key];rtol=1e-12,atol=1e-13)
                end
            end
            global initial=current
        end
        m=TOML.parsefile(joinpath(out,"metadata.toml"))
        hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*48
        @test isfinite(hours) && hours>0
        n==512 && (@test hours<4.0)
        push!(estimates,hours)
        run(`tools/field-export-env/bin/python scripts/ferrel2026_v3/transfer.py $out`)
    end
    write_gate(joinpath(root,"PASS$(n).toml");grid=[n,n,96],estimated_2day_wall_hours=estimates,
        pilot_simulated_hours=6,throughput_caveat="Early flow only; not mature convection, convergence or a steady Ferrel cell",
        gpu_name=CUDA.name(CUDA.device()))
end
