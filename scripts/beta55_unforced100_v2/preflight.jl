include("run.jl")
include("check_output.jl")
root=only(ARGS);mkpath(root)
base=TOML.parsefile("configs/beta55_unforced100.toml")
reference=TOML.parsefile("configs/cascadebeta100_heated.toml")
@test Set(keys(base)) == union(Set(keys(reference)), Set(["forcing_mode"]))
changed=Set(["name","heat_flux","cooling_model","balance_heat","forcing_multiplier"])
for key in keys(reference)
    key in changed && continue
    @test base[key] == reference[key]
end
@test base["forcing_mode"] == "none"
@test base["heat_flux"] == 0 && base["cooling_model"] == "none" && !base["balance_heat"]
@test base["beta"] ≈ 2base["rotation_rate"]*cosd(base["latitude"])/base["planet_radius"]
@test base["rotation_model"] == "beta-plane" && base["latitude"] == 55.0
@test read(joinpath(@__DIR__,"initial_flow.jl")) == read("scripts/cascadebeta100/initial_flow.jl")
reference_dir="data/cascadebeta100-v1/cascadebeta100_heated"
@test TOML.parsefile(joinpath(reference_dir,"segment_10","metadata.toml"))["status"] == "COMPLETE"
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_shear_output(out);GC.gc();CUDA.reclaim()
    return out
end
small=Dict("nx"=>64,"ny"=>64,"nz"=>64,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0,"turbulence_wavelength"=>2400000.0)
full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
@testset "Beta-plane restart" begin
    a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","theta","theta_xz","w_xz")
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
end
out=experiment("day1_beta",Dict("stop_seconds"=>86400.0,"checkpoint_seconds"=>86400.0))
@testset "Actual saved initial fields match completed heated beta-plane reference" begin
    a=load(joinpath(out,"frame_00000.jld2"))
    b=load(joinpath(reference_dir,"segment_1","frame_00000.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","theta","theta_xz","w_xz")
        @test a[key]==b[key]
    end
end
m=TOML.parsefile(joinpath(out,"metadata.toml"))
hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*240
@test hours<4
run(`tools/field-export-env/bin/python scripts/beta55_unforced100_v2/transfer.py $out`)
hashes=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if isfile(f))
open(joinpath(root,"PASS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","initial_fields_equal_reference"=>true,"estimated_10day_wall_hours"=>hours,"source_hashes"=>hashes,"configuration_sha256"=>bytes2hex(sha256(read("configs/beta55_unforced100.toml"))),"environment_hashes"=>Dict(f=>bytes2hex(sha256(read(f))) for f in ("Project.toml","Manifest.toml"))))
end
