include("common.jl")
include("run.jl")
include("check_native.jl")
include("restart_comparison.jl")
module FrozenReference
include("../cascade35/run.jl")
end
root=ARGS[1];validation_arch=length(ARGS)==2 ? ARGS[2] : "GPU"
validation_arch in ["CPU","GPU"] || error("Validation architecture")
mkpath(root)
validation_arch=="GPU" && verify_gate(joinpath(root,"cpu","PASS.toml"))
base=TOML.parsefile(joinpath(CONFIG_ROOT,"diagnostic.toml"));wide=TOML.parsefile(joinpath(CONFIG_ROOT,"wide.toml"))
original=TOML.parsefile("configs/cascade35_heated.toml")
@testset "Only authorized physical differences" begin
    for (k,v) in original
        k in ["name","stop_seconds"] && continue
        @test base[k]==v
        k in ["ny","ly"] && continue
        @test wide[k]==v
    end
    @test wide["ny"]==2base["ny"] && wide["ly"]==2base["ly"]
    @test base["nx"]==wide["nx"]==256
    @test base["lx"]/base["nx"]==wide["lx"]/wide["nx"]
    @test base["ly"]/base["ny"]==wide["ly"]/wide["ny"]
    @test base["stop_seconds"]==210*86400.0 && wide["stop_seconds"]==200*86400.0
    @test read(joinpath(@__DIR__,"initial_flow.jl"))==read("scripts/cascade35/initial_flow.jl")
    @test read(joinpath(@__DIR__,"thermal_budget_controller.jl"))==read("scripts/cascade35/thermal_budget_controller.jl")
end
function experiment(name,cfg;pickup=false,reference=false)
    c=copy(cfg);c["name"]=name
    config=writeconfig(root,name,c);out=joinpath(root,name)
    args=[config,out,validation_arch];pickup===false || push!(args,pickup)
    if reference
        FrozenReference.run_eddy(args);check_shear_output(out)
    else
        run_eddy(args);check_native_output(out)
    end
    GC.gc();validation_arch=="GPU" && CUDA.reclaim();out
end
function compare_fields(a,b)
    for key in ["u","v","w","theta","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","zonal_u","zonal_v","zonal_w","eddy_uv","eddy_uw","native_rhov","native_rhow"]
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
    @test streamfunction_restart_report(a,b,base["height"])["passed"]
end
lastframe(out)=load(last(sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))))
small=merge(wide,Dict("nx"=>64,"ny"=>128,"nz"=>32,"turbulence_wavelength"=>2400000.0,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0))
full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
@testset "Restart fields and integrated native terms" begin
    a,b,c=lastframe(full),lastframe(part),lastframe(resumed)
    compare_fields(a,c)
    @test isapprox(a["native_integrated_parts"],b["native_integrated_parts"]+c["native_integrated_parts"];rtol=1e-8,atol=1e-9)
end
if validation_arch=="CPU"
    gate(joinpath(root,"PASS.toml");architecture="CPU",restart_comparison="PASS",cumulative_native_closure="PASS")
    println("CPU numerical validation PASSED")
    exit()
end
# A mature actual day-200 checkpoint tests that new observers leave the physical evolution unchanged.
old=TOML.parsefile(joinpath(ORIGINAL_DAY200,"metadata.toml"));@test old["status"]=="COMPLETE" && old["simulated_seconds"]==200*86400.0
checkpoint=latest_checkpoint(ORIGINAL_DAY200)
probe=merge(base,Dict("stop_seconds"=>200*86400.0+1200,"output_seconds"=>600.0,"diagnostic_seconds"=>600.0,"checkpoint_seconds"=>600.0))
baseline=experiment("mature_original",probe;pickup=checkpoint,reference=true)
instrumented=experiment("mature_instrumented",probe;pickup=checkpoint)
@testset "Diagnostic observer leaves mature evolution unchanged" begin
    compare_fields(lastframe(baseline),lastframe(instrumented))
end
# Full-size runtime/memory and output check in the wider channel.
pilot=experiment("wide_full_grid",merge(wide,Dict("stop_seconds"=>21600.0,"output_seconds"=>3600.0,"checkpoint_seconds"=>21600.0)))
m=TOML.parsefile(joinpath(pilot,"metadata.toml"))
hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*240
@test hours<8
@test isapprox(m["initial_horizontal_rms_lower_m_s"],old["initial_horizontal_rms_lower_m_s"];rtol=1e-12)
run(`tools/field-export-env/bin/python $(joinpath(@__DIR__,"transfer.py")) $pilot`)
gate(joinpath(root,"PASS.toml");estimated_10day_wide_hours=hours,mature_checkpoint_sha256=bytes2hex(open(sha256,checkpoint)),mature_physics_comparison="PASS",restart_comparison="PASS",every_step_native_closure="PASS",full_grid=[256,512,64])
println("FERREL WIDTH PREFLIGHT PASSED; estimated wide 10-day hours=",hours)
