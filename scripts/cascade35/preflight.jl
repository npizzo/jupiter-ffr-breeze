include("run.jl")
include("check_output.jl")
include("restart_comparison.jl")
root=only(ARGS);mkpath(root)
base=TOML.parsefile("configs/cascade35_heated.toml")
reference=TOML.parsefile("configs/cascade100_f_heated.toml")
@test Set(keys(base))==Set(keys(reference))
for key in keys(reference)
    key in ("name","rotation_model","beta","latitude","turbulence_rossby") && continue
    @test base[key]==reference[key]
end
@test base["beta"] ≈ 2base["rotation_rate"]*cosd(base["latitude"])/base["planet_radius"]
@test base["rotation_model"]=="beta-plane"
oldpass=TOML.parsefile("data/cascade100-preflight-v1/PASS.toml")
for file in ("initial_flow.jl","thermal_budget_controller.jl","check_output.jl","check_large.jl","transfer.py")
    @test bytes2hex(sha256(read(joinpath(@__DIR__,file))))==oldpass["source_hashes"][file]
end
# Strip only the three output additions and verify the unchanged dynamical source.
source=read(joinpath(@__DIR__,"run.jl"),String)
source=replace(source,"include(\"zonal_output.jl\")\n"=>"")
source=replace(source,"    zonal = zonal_output(m, uc, vc, wcenter, th, dye, Ly/Ny, dz)\n"=>"")
source=replace(source,"jldsave(file; zonal..., ubar,vbar"=>"jldsave(file; ubar,vbar")
@test bytes2hex(sha256(source))==oldpass["source_hashes"]["run.jl"]
@test base["latitude"]==35.0
@test base["turbulence_rossby"]*sind(base["latitude"]) ≈ reference["turbulence_rossby"]*sind(reference["latitude"])
reference_dir="data/cascade100-v1/cascade100_f_heated"
@test TOML.parsefile(joinpath(reference_dir,"segment_10","metadata.toml"))["status"]=="COMPLETE"
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_shear_output(out);check_zonal_frames(out);GC.gc();CUDA.reclaim()
    return out
end
small=Dict("nx"=>64,"ny"=>64,"nz"=>64,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0,"turbulence_wavelength"=>2400000.0)
full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
@testset "Beta-plane restart" begin
    a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","zonal_u","zonal_v","zonal_w","eddy_uv","eddy_uw","native_rhov","native_rhow")
        @testset "$key" begin
            @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
        end
    end
    report=streamfunction_restart_report(a,b,base["height"])
    @testset "Depth-integrated mass streamfunction" begin
        @test report["native_flux_pass"]
        @test report["depth_normalized_pass"]
        @test report["propagated_error_pass"]
    end
    open(joinpath(root,"restart_streamfunction.toml"),"w") do io;TOML.print(io,report);end
end
out=experiment("day1_beta",Dict("stop_seconds"=>86400.0,"checkpoint_seconds"=>86400.0))
@testset "Actual saved initial fields match completed f-plane reference" begin
    a=load(joinpath(out,"frame_00000.jld2"))
    b=load(joinpath(reference_dir,"segment_1","frame_00000.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy")
        @test isapprox(a[key],b[key];rtol=1e-12,atol=1e-13)
    end
end
m=TOML.parsefile(joinpath(out,"metadata.toml"))
hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*240
@test hours<4
run(`tools/field-export-env/bin/python scripts/cascade35/transfer.py $out`)
hashes=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if isfile(f))
open(joinpath(root,"PASS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","initial_fields_match_reference_tolerance"=>"rtol=1e-12, atol=1e-13; pressure estimate uses local f","estimated_10day_wall_hours"=>hours,"source_hashes"=>hashes,"configuration_sha256"=>bytes2hex(sha256(read("configs/cascade35_heated.toml"))),"environment_hashes"=>Dict(f=>bytes2hex(sha256(read(f))) for f in ("Project.toml","Manifest.toml"))))
end
