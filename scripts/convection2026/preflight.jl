include("run.jl")
include("check_output.jl")
include("common.jl")
length(ARGS)==2 || error("usage: preflight.jl ROOT GRID")
root=ARGS[1];n=parse(Int,ARGS[2]);n in (512,1024) || error("Invalid pilot grid")
verify_gate(joinpath(root,"UNIT_PASS.toml"))
n==1024 && verify_gate(joinpath(root,"PASS512.toml"))
base=TOML.parsefile(config_path(162))
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    isfile(config) && error("Refusing repeated pilot: $name")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_shear_output(out)
    for path in filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true))
        jldopen(path) do d
            for key in ("wcenter_slices","theta_slices")
                @test size(d[key])==(cfg["nx"],cfg["ny"],3)
                @test all(isfinite,d[key])
            end
            for key in ("zonal_u","zonal_v","zonal_w","zonal_theta","zonal_tracer","zonal_rho_uv","zonal_rho_uw","zonal_rho_vtheta","zonal_rho_wtheta","zonal_rho_vtracer","zonal_rho_wtracer","zonal_theta_source")
                @test size(d[key])==(cfg["ny"],cfg["nz"])
                @test all(isfinite,d[key])
            end
        end
    end
    GC.gc();CUDA.reclaim()
    return out
end
if n==512
    small=Dict("nx"=>64,"ny"=>64,"nz"=>96,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0,"turbulence_wavelength"=>2400000.0)
    full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
    part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
    resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
    @testset "96-level checkpoint continuation" begin
        a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","zonal_u","zonal_rho_uv","zonal_theta_source")
            @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
        end
    end
end
initial=nothing;estimates=Float64[]
for heat in (n==512 ? HEATS : (162,))
    c=TOML.parsefile(config_path(heat))
    out=experiment("pilot$(n)_f$(heat)",merge(c,Dict("nx"=>n,"ny"=>n,"nz"=>96,"stop_seconds"=>21600.0,"output_seconds"=>3600.0,"diagnostic_seconds"=>600.0,"checkpoint_seconds"=>10800.0)))
    current=load(joinpath(out,"frame_00000.jld2"))
    if initial!==nothing
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
            @test current[key]==initial[key]
        end
    end
    global initial=current
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*48
    @test isfinite(hours) && hours>0
    # Early-time throughput only. Later convection may shorten the timestep.
    n==512 && (@test hours<4.0)
    push!(estimates,hours)
    run(`tools/field-export-env/bin/python scripts/convection2026/transfer.py $out`)
end
write_gate(joinpath(root,"PASS$(n).toml");grid=[n,n,96],estimated_2day_wall_hours=estimates,
    pilot_simulated_hours=6,throughput_caveat="Six-hour initial-state pilot; not mature-convection throughput or resolution convergence",gpu_name=CUDA.name(CUDA.device()))
