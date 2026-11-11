include("run.jl")
include("check_output.jl")
include("common.jl")
length(ARGS)==4 || error("usage: campaign.jl HEAT SEGMENT ROOT PREFLIGHT_ROOT")
heat=parse(Int,ARGS[1]);segment=parse(Int,ARGS[2]);heat in HEATS || error("Invalid heating")
segment in 1:5 || error("First batch is limited to ten days")
root,preflight=ARGS[3:4];verify_gate(joinpath(preflight,"PASS512.toml"))
cfg=TOML.parsefile(config_path(heat));cfg["stop_seconds"]=segment*2*86400.0
name=cfg["name"];out=joinpath(root,name,"segment_$(segment)");mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs","$(name)_segment_$(segment).toml")
isfile(config) && error("Refusing duplicate segment configuration")
open(config,"w") do io;TOML.print(io,cfg);end
args=[config,out,"GPU"]
if segment>1
    previous=joinpath(root,name,"segment_$(segment-1)")
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
