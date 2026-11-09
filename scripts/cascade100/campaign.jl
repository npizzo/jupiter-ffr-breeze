include("run.jl")
include("check_output.jl")
length(ARGS)==3 || error("usage: campaign.jl INDEX ROOT PREFLIGHT_ROOT")
index=parse(Int,ARGS[1]); index in 1:20 || error("index must be 1:20")
root,preflight=ARGS[2:3]
pass=TOML.parsefile(joinpath(preflight,"PASS.toml"))
pass["status"]=="PASS" || error("Preflight incomplete")
for (file,hash) in pass["source_hashes"]
    bytes2hex(sha256(read(joinpath(@__DIR__,file))))==hash || error("Source changed after preflight: $file")
end
names=("cascade100_f_heated","cascade100_f_control")
name=names[mod1(index,2)];segment=(index-1)÷2+1
config_source=joinpath("configs",name*".toml")
bytes2hex(sha256(read(config_source)))==pass["configuration_hashes"][name] || error("Configuration changed after preflight")
cfg=TOML.parsefile(config_source);cfg["stop_seconds"]=segment*10*86400.0
out=joinpath(root,name,"segment_"*string(segment));mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs",name*"_segment_"*string(segment)*".toml")
open(config,"w") do io;TOML.print(io,cfg);end
args=[config,out,"GPU"]
if segment>1
    previous=joinpath(root,name,"segment_"*string(segment-1))
    old=TOML.parsefile(joinpath(previous,"metadata.toml"))
    old["status"]=="COMPLETE" || error("Previous segment incomplete")
    for (k,v) in cfg
        k=="stop_seconds" && continue
        old[k]==v || error("Restart configuration changed: $k")
    end
    for (key,file) in [("script_sha256","run.jl"),("initial_channel_script_sha256","initial_flow.jl"),("controller_sha256","thermal_budget_controller.jl")]
        old[key]==bytes2hex(sha256(read(joinpath(@__DIR__,file)))) || error("Restart source changed")
    end
    push!(args,latest_checkpoint(previous))
end
run_eddy(args);check_shear_output(out)
if segment==1 && index==2
    a=load(joinpath(root,first(names),"segment_1","frame_00000.jld2"))
    b=load(joinpath(out,"frame_00000.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
        @test a[key]==b[key]
    end
end
