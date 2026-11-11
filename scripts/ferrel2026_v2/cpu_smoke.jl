# Run only on Unity. This small CPU check is a required predecessor of the v2 GPU gates.
# Import the frozen driver and checks without changing their source or gate files.
include("run.jl")
include("check_ferrel.jl")
include("common.jl")

length(ARGS)==2 || error("usage: smoke.jl NEW_ROOT MAIN_PREFLIGHT_ROOT")
root,main_preflight=abspath.(ARGS)
isdir(root) && !isempty(readdir(root)) && error("CPU smoke root must be fresh")
mkpath(root)
verify_gate(joinpath(main_preflight,"UNIT_PASS.toml"))

small=Dict{String,Any}("nx"=>64,"ny"=>64,"nz"=>24,
    "turbulence_wavelength"=>6000000.0,"initial_dt"=>1.0,"max_dt"=>5.0,
    "output_seconds"=>60.0,"diagnostic_seconds"=>10.0,"controller_seconds"=>10.0,
    "checkpoint_seconds"=>60.0,"native_budget_interval"=>1,"native_budget_first_steps"=>8)

function experiment(name,case,stop_seconds;pickup=false)
    cfg=merge(TOML.parsefile(config_path(case)),small,Dict("stop_seconds"=>Float64(stop_seconds)))
    config=joinpath(root,name*".toml")
    isfile(config) && error("Refusing repeated CPU smoke experiment")
    open(config,"w") do io;TOML.print(io,cfg);end
    output=joinpath(root,name)
    args=[config,output,"CPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args)
    check_ferrel_output(output)
    metadata=TOML.parsefile(joinpath(output,"metadata.toml"))
    @test metadata["architecture"]=="CPU"
    @test metadata["native_budget_checks"]>=10
    @test metadata["final_dt"]<=5.0
    @test metadata["initial_zonal_pressure_gradient_post_halo_fill"]<1e-10
    println("Pressure initialization native mean before/after halo fill: ",
        metadata["initial_zonal_pressure_gradient_pre_halo_fill"]," / ",
        metadata["initial_zonal_pressure_gradient_post_halo_fill"])
    GC.gc()
    return output
end

started=time()
full=experiment("continuous","b162",120)
first_segment=experiment("first_segment","b162",60)
resumed=experiment("resumed","b162",120;pickup=latest_checkpoint(first_segment))

@testset "CPU native closure and restart smoke" begin
    a=load(joinpath(full,"frame_00002.jld2"))
    b=load(joinpath(resumed,"frame_00002.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar",
                "native_rhou","native_rhov","native_rhow","mass_streamfunction","native_rhou_rhs")
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
    for output in (full,first_segment,resumed)
        m=TOML.parsefile(joinpath(output,"metadata.toml"))
        # Every completed step is sampled, including the initial step after restart.
        rows=readlines(joinpath(output,"native_budget","checks.csv"))
        iterations=[parse(Int,first(split(row,','))) for row in rows[2:end]]
        @test all(diff(iterations).==1)
        @test length(iterations)==m["native_budget_checks"]
    end
end

# These zero-noise controls exercise jet balance and symmetry through the same APIs.
for case in ("axis0","axis162")
    experiment(case,case,120)
end

write_gate(joinpath(root,"CPU_SMOKE_PASS.toml");
    architecture="CPU",grid=[64,64,24],max_dt_seconds=5.0,simulated_seconds=120,
    wall_seconds=time()-started,
    smoke_script_sha256=bytes2hex(sha256(read(@__FILE__))),
    frozen_unit_gate=joinpath(main_preflight,"UNIT_PASS.toml"),
    scope="Independent API, native SSPRK3 closure, restart and axisymmetric smoke check; not a physical-resolution gate and is required before the v2 GPU gates")
println("Independent CPU smoke PASS: ",joinpath(root,"CPU_SMOKE_PASS.toml"))
