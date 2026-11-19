using TOML, SHA, CUDA
started_job=parse(Float64,ENV["FFR_JOB_START_EPOCH"])
pass=TOML.parsefile("data/cuda122-probe-v1/PASS.toml")
@assert pass["status"]=="PASS"
for (file,h) in pass["hashes"]
    bytes2hex(sha256(read(file)))==h || error("CUDA runtime changed: $file")
end
CUDA.functional() || error("GPU required")
@assert v"12.2" <= CUDA.runtime_version() < v"12.3"
gpu=CUDA.name(CUDA.device())
(occursin("A100",gpu)||occursin("H100",gpu)) || error("A100/H100 required, got $gpu")
for (file,h) in TOML.parsefile("scripts/strong256_v1/frozen-source.toml")
    bytes2hex(sha256(read(file)))==h || error("Frozen model changed: $file")
end
include("../shear_layer_v1/run.jl")
root,case,oldroot,prior_seconds_string=ARGS
prior_seconds=parse(Int,prior_seconds_string)
@assert 0 < prior_seconds <= 600
@assert case in ("stable_none","stable_cyclonic","stable_anticyclonic")
ispath(root) && error("Refuse duplicate output")
cfg=TOML.parsefile("configs/strong256-v1/$case.toml")
old=joinpath(oldroot,"day0_1")
oldmeta=TOML.parsefile(joinpath(old,"metadata.toml"))
@assert oldmeta["status"]=="COMPLETE" && oldmeta["simulated_seconds"]==86400
for (k,v) in cfg
    k in ("stop_seconds","wall_time_limit") && continue
    oldmeta[k]==v || error("Configuration changed: $k")
end
@assert oldmeta["script_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"../shear_layer_v1/run.jl"))))
checkpoint=only(filter(f->endswith(f,".jld2"),readdir(joinpath(old,"checkpoints");join=true)))
mkpath(root)
function save_toml(path,data)
    open(path,"w") do io;TOML.print(io,data);end
end
save_toml(joinpath(root,"experiment.toml"),Dict("case"=>case,"gpu"=>gpu,
    "oldroot"=>abspath(oldroot),"checkpoint"=>abspath(checkpoint),
    "checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint)),
    "prior_allocated_gpu_seconds"=>prior_seconds,"new_allocation_limit_seconds"=>21000,
    "combined_allocation_limit_seconds"=>prior_seconds+21000,
    "physical_changes"=>String[],"numerical_configuration_changes"=>String[],
    "interpretation"=>"Direct continuation from saved day 1 to day 50; no intermediate runtime gate"))
# Reuse the finished first-day output for the existing two-segment movie assembler.
symlink(abspath(old),joinpath(root,"day0_1"))
cfg["wall_time_limit"]=20700-(time()-started_job)
prodconfig=joinpath(root,"continuation.toml");save_toml(prodconfig,cfg)
println("DIRECT_DAY1_TO50_START gpu=",gpu);flush(stdout)
run_eddy([prodconfig,joinpath(root,"day1_50"),"GPU",checkpoint])
final=TOML.parsefile(joinpath(root,"day1_50","metadata.toml"))
@assert final["status"]=="COMPLETE" && final["segment_start_seconds"]==86400 && final["simulated_seconds"]==50*86400
for key in ("vorticity_slices","tracer_slices","ubar","vbar")
    a=load(joinpath(old,"frame_00008.jld2"),key)
    b=load(joinpath(root,"day1_50","frame_00008.jld2"),key)
    maximum(abs,a.-b)<=1e-10+1e-6*maximum(abs,a) || error("Restart field changed: $key")
end
save_toml(joinpath(root,"COMPLETE.toml"),Dict("status"=>"COMPLETE","days"=>50,
    "total_new_wall_seconds"=>time()-started_job,"prior_allocated_gpu_seconds"=>prior_seconds,
    "oldroot"=>abspath(oldroot),"restart_fields_check"=>"PASS","segments"=>["day0_1","day1_50"]))
println("STRONG256_DIRECT_FIFTY_DAYS_COMPLETE")
