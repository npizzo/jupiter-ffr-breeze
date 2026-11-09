include("run.jl")
include("check_output.jl")
length(ARGS)==3 || error("usage: campaign.jl INDEX ROOT PREFLIGHT_ROOT")
index=parse(Int,ARGS[1]); index in 1:8 || error("index must be 1:8")
root,preflight=ARGS[2:3]
TOML.parsefile(joinpath(preflight,"PASS.toml"))["status"]=="PASS" || error("Preflight incomplete")
names=("eddy20_f_heated","eddy20_f_control","eddy20_beta_heated","eddy20_beta_control")
name=names[mod1(index,4)];segment=(index-1)÷4+1
cfg=TOML.parsefile(joinpath("configs",name*".toml"));cfg["stop_seconds"]=segment*10*86400.0
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
    old["script_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"run.jl")))) || error("Driver changed between segments")
    old["initial_channel_script_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"initial_flow.jl")))) || error("Initializer changed between segments")
    old["controller_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"thermal_budget_controller.jl")))) || error("Controller changed between segments")
    push!(args,latest_checkpoint(previous))
end
run_eddy(args);check_shear_output(out)
if segment==1 && index>1
    a=load(joinpath(root,first(names),"segment_1","frame_00000.jld2"))
    b=load(joinpath(out,"frame_00000.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices")
        @test a[key]==b[key]
    end
end
