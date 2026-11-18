# Continue the frozen six-case experiment; never change its running source/configs.
using TOML,SHA,CUDA,JLD2
pass=TOML.parsefile("data/cuda122-probe-v1/PASS.toml");pass["status"]=="PASS"||error("CUDA gate")
for (f,h) in pass["hashes"];bytes2hex(sha256(read(f)))==h||error("Runtime changed");end
CUDA.functional()||error("GPU required");@assert v"12.2"<=CUDA.runtime_version()<v"12.3"
include("../shear_layer_v1/run.jl")
oldroot,newroot,case=ARGS
case in ["$(a)_$(b)" for a in ("stable","neutral") for b in ("cyclonic","none","anticyclonic")]||error("Unknown case")
p=TOML.parsefile(joinpath(oldroot,"preflight","PASS.toml"));p["status"]=="PASS"||error("Missing gate")
for (file,h) in p["source_hashes"]
 bytes2hex(sha256(read(file)))==h||error("Frozen source changed: $file")
end
cfg=TOML.parsefile("configs/shear-layer-v1/$case.toml")
previous=joinpath(oldroot,case);m=TOML.parsefile(joinpath(previous,"metadata.toml"))
m["status"]=="COMPLETE"||error("Ten-day segment incomplete")
@assert isapprox(m["simulated_seconds"],10*86400.;atol=1e-6)
for (k,v) in cfg
 m[k]==v||error("Checkpoint configuration differs: $k")
end
m["script_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"../shear_layer_v1/run.jl"))))||error("Driver changed")
checkpoints=filter(f->endswith(f,".jld2"),readdir(joinpath(previous,"checkpoints");join=true))
length(checkpoints)==1||error("Expected one retained final checkpoint")
checkpoint=only(checkpoints)
mkpath(joinpath(newroot,"configs"))
config=joinpath(newroot,"configs",case*".toml");isfile(config)&&error("Duplicate continuation config")
cfg["stop_seconds"]=50*86400.;cfg["wall_time_limit"]=8400.
open(config,"w") do io;TOML.print(io,cfg);end
open(joinpath(newroot,"configs",case*"-restart.toml"),"w") do io
 TOML.print(io,Dict("checkpoint"=>abspath(checkpoint),"checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint)),
 "source_metadata_sha256"=>bytes2hex(sha256(read(joinpath(previous,"metadata.toml")))) ,
 "start_days"=>10,"target_days"=>50,"changed_parameters"=>["stop_seconds","wall_time_limit"]))
end
out=joinpath(newroot,case)
run_eddy([config,out,"GPU",checkpoint])
m2=TOML.parsefile(joinpath(out,"metadata.toml"))
@assert isapprox(m2["segment_start_seconds"],10*86400.;atol=1e-6)
@assert isapprox(m2["simulated_seconds"],50*86400.;atol=1e-6)
println("SHEAR_LAYER_FIFTY_DAY_CASE_COMPLETE ",case);flush(stdout)
