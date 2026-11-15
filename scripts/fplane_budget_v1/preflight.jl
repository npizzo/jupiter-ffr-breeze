using Test,TOML,SHA
include("run.jl")
root=only(ARGS);out=joinpath(root,"preflight");mkpath(out)
TOML.parsefile("data/fplane-budget-cpu-pilot-v2/PASS.toml")["status"]=="PASS"||error("CPU observer gate")
# Repeat the unchanged-state observer test on the GPU backend.
for audit in [false,true]
 cfg=TOML.parsefile("configs/fplane-point4-v2/early_replay_heated.toml")
 cfg["nx"]=64;cfg["ny"]=64;cfg["turbulence_wavelength"]=2400000.;cfg["stop_seconds"]=120.;cfg["output_seconds"]=60.;cfg["diagnostic_seconds"]=60.;cfg["full_field_seconds"]=60.;cfg["max_dt"]=10.;cfg["initial_dt"]=10.;cfg["native_budget"]=audit;cfg["reduce_exports"]=false;cfg["name"]=audit ? "audit" : "baseline"
 path=joinpath(out,cfg["name"]*".toml");open(path,"w") do io;TOML.print(io,cfg);end
 run_eddy([path,joinpath(out,cfg["name"]),"GPU"]);GC.gc();CUDA.reclaim()
end
name=last(sort(readdir(joinpath(out,"baseline","fields3d"))));a=load(joinpath(out,"baseline","fields3d",name));b=load(joinpath(out,"audit","fields3d",name))
for k in ["u_face","v_face","w_face","theta","dye_density"];@test a[k]==b[k];end
a=nothing;b=nothing;GC.gc();CUDA.reclaim()
resources=Dict{String,Any}()
for nx in [512,1024]
 cfg=TOML.parsefile("configs/fplane-budget-v1/mature$(nx).toml")
 cfg["stop_seconds"]=90*86400.0+1200.;cfg["output_seconds"]=600.
 path=joinpath(out,"mature$(nx).toml");open(path,"w") do io;TOML.print(io,cfg);end
 dir=joinpath(out,"mature$(nx)");run_eddy([path,dir,"GPU"])
 meta=TOML.parsefile(joinpath(dir,"metadata.toml"));@test meta["status"]=="COMPLETE"
 @test isempty(readdir(joinpath(dir,"fields3d")))
 @test length(filter(f->endswith(f,".json"),readdir(joinpath(dir,"reduced"))))==3
 # Independent scalar and energy remap accounting against the saved coarse state.
 frame=load(joinpath(dir,@sprintf("frame_%05d.jld2",Int(90*86400/600))))
 resources[string(nx)]=Dict("initial_ke"=>frame["kinetic_energy"],"initial_theta_profile"=>frame["theta_profile"],"initial_tracer_column"=>frame["tracer_column"],"estimated_18h_runtime_hours"=>meta["post_warmup_wall_seconds"]/meta["post_warmup_simulated_seconds"]*18)
 @test resources[string(nx)]["estimated_18h_runtime_hours"]<14
 GC.gc();CUDA.reclaim()
end
@test abs(resources["1024"]["initial_ke"]/resources["512"]["initial_ke"]-1)<.1
@test isapprox(resources["1024"]["initial_theta_profile"],resources["512"]["initial_theta_profile"];rtol=1e-12)
@test isapprox(resources["1024"]["initial_tracer_column"],resources["512"]["initial_tracer_column"];rtol=1e-12)
files=["run.jl","remap.jl","initial_flow.jl","thermal_budget_controller.jl","point4_state.jl","native_vorticity.jl","diagnostics.py","reduce.py","cuda_launch.jl","preflight.jl","campaign.jl"]
pass=Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"resources"=>resources,"bitwise_observer_test"=>true,"source_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in files),"config_hashes"=>Dict("configs/fplane-budget-v1/mature$(n).toml"=>bytes2hex(sha256(read("configs/fplane-budget-v1/mature$(n).toml"))) for n in [512,1024]),"bundle_sha256"=>bytes2hex(open(sha256,"data/fplane-point4-v2/preflight/day90_heated_bundle.jld2")))
open(joinpath(out,"GPU_PASS.toml"),"w") do io;TOML.print(io,pass);end
println("NATIVE_BUDGET_FULL_GRID_GPU_PASS");flush(stdout)
