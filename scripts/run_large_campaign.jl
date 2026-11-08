include("run_large.jl")
include("check_large.jl")
length(ARGS)==3 || error("usage: run_large_campaign.jl INDEX ROOT PREFLIGHT_ROOT")
index=parse(Int,ARGS[1]); index in 1:8 || error("index must be 1:8")
root,preflight=ARGS[2:3]
TOML.parsefile(joinpath(preflight,"PASS.toml"))["status"]=="PASS" || error("Preflight did not pass")
case=isodd(index) ? "f" : "beta"; segment=(index+1)÷2
name="large40_channel_"*case
cfg=TOML.parsefile(joinpath("configs",name*".toml"))
cfg["stop_seconds"]=segment*10*86400.0
out=joinpath(root,name,"segment_"*string(segment))
mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs",name*"_segment_"*string(segment)*".toml")
open(config,"w") do io; TOML.print(io,cfg); end
args=[config,out,"GPU"]
if segment>1
    previous=joinpath(root,name,"segment_"*string(segment-1))
    old=TOML.parsefile(joinpath(previous,"metadata.toml"))
    old["status"]=="COMPLETE" || error("Previous segment incomplete")
    for (k,v) in cfg
        k=="stop_seconds" && continue
        old[k]==v || error("Restart configuration changed: $k")
    end
    push!(args,latest_checkpoint(previous))
end
run_large(args)
check_large_output(out)
if case=="beta" && segment==1
    a=load(joinpath(root,"large40_channel_f","segment_1","frame_00000.jld2"))
    b=load(joinpath(out,"frame_00000.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","ubar","vbar")
        @test a[key]==b[key]
    end
end
