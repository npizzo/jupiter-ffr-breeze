# Continue the frozen 35 N model; never edit or overwrite the day 0-100 campaign.
using TOML, SHA
source_dir=joinpath(@__DIR__,"..","cascade35")
length(ARGS)==4 || error("usage: campaign.jl SEGMENT EXTENSION_ROOT ORIGINAL_ROOT PREFLIGHT")
segment=parse(Int,ARGS[1]);segment in 11:30 || error("Authorized extension is segments 11:30")
root,original,preflight=ARGS[2:4]
pass=TOML.parsefile(joinpath(preflight,"PASS.toml"));pass["status"]=="PASS" || error("Preflight incomplete")
for (file,hash) in pass["source_hashes"]
    bytes2hex(sha256(read(joinpath(source_dir,file))))==hash || error("Original source changed: $file")
end
for (file,hash) in pass["environment_hashes"]
    bytes2hex(sha256(read(file)))==hash || error("Environment changed: $file")
end
frozen=TOML.parsefile(joinpath(@__DIR__,"FROZEN.toml"))
for (file,hash) in frozen["files"]
    bytes2hex(sha256(read(joinpath(@__DIR__,file))))==hash || error("Extension source changed: $file")
end
include(joinpath(source_dir,"run.jl"))
include(joinpath(source_dir,"check_output.jl"))
name="cascade35_heated";config_source="configs/"*name*".toml"
bytes2hex(sha256(read(config_source)))==pass["configuration_sha256"] || error("Configuration changed")
cfg=TOML.parsefile(config_source);cfg["stop_seconds"]=segment*10*86400.0
previous_root=segment==11 ? original : root
previous=joinpath(previous_root,name,"segment_"*string(segment-1))
old=TOML.parsefile(joinpath(previous,"metadata.toml"))
old["status"]=="COMPLETE" || error("Previous segment incomplete")
old["simulated_seconds"]==(segment-1)*10*86400.0 || error("Wrong restart time")
for (k,v) in cfg
    k=="stop_seconds" && continue
    old[k]==v || error("Restart configuration changed: $k")
end
for (key,file) in [("script_sha256","run.jl"),("initial_channel_script_sha256","initial_flow.jl"),("controller_sha256","thermal_budget_controller.jl")]
    old[key]==pass["source_hashes"][file] || error("Previous source mismatch: $file")
end
checkpoint=latest_checkpoint(previous)
out=joinpath(root,name,"segment_"*string(segment));mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs",name*"_segment_"*string(segment)*".toml")
open(config,"w") do io;TOML.print(io,cfg);end
run_eddy([config,out,"GPU",checkpoint])
check_shear_output(out);check_zonal_frames(out)
current=TOML.parsefile(joinpath(out,"metadata.toml"))
@test current["segment_start_seconds"]==(segment-1)*10*86400.0
@test current["simulated_seconds"]==segment*10*86400.0
open(joinpath(out,"continuation.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","source_checkpoint"=>abspath(checkpoint),"source_checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint)),"start_day"=>(segment-1)*10,"end_day"=>segment*10,"extension_driver_sha256"=>bytes2hex(sha256(read(@__FILE__)))))
end
