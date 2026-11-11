include("run.jl")
include("check_output.jl")
length(ARGS)==3 || error("usage: campaign.jl SEGMENT ROOT PREFLIGHT")
segment=parse(Int,ARGS[1]);segment in 1:10 || error("segment must be 1:10")
root,preflight=ARGS[2:3]
pass=TOML.parsefile(joinpath(preflight,"PASS.toml"));pass["status"]=="PASS" || error("Preflight incomplete")
for (file,hash) in pass["source_hashes"]
    bytes2hex(sha256(read(joinpath(@__DIR__,file))))==hash || error("Source changed: $file")
end
for (file,hash) in pass["environment_hashes"]
    bytes2hex(sha256(read(file)))==hash || error("Environment changed: $file")
end
name="cascade35_heated";config_source="configs/"*name*".toml"
bytes2hex(sha256(read(config_source)))==pass["configuration_sha256"] || error("Configuration changed")
cfg=TOML.parsefile(config_source);cfg["stop_seconds"]=segment*10*86400.0
out=joinpath(root,name,"segment_"*string(segment));mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs",name*"_segment_"*string(segment)*".toml")
open(config,"w") do io;TOML.print(io,cfg);end
args=[config,out,"GPU"]
if segment>1
    previous=joinpath(root,name,"segment_"*string(segment-1))
    old=TOML.parsefile(joinpath(previous,"metadata.toml"));old["status"]=="COMPLETE" || error("Previous segment incomplete")
    for (k,v) in cfg
        k=="stop_seconds" && continue
        old[k]==v || error("Restart configuration changed: $k")
    end
    push!(args,latest_checkpoint(previous))
end
run_eddy(args);check_shear_output(out);check_zonal_frames(out)
if segment==1
    a=load(joinpath(out,"frame_00000.jld2"));b=load("data/cascade100-v1/cascade100_f_heated/segment_1/frame_00000.jld2")
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy")
        @test isapprox(a[key],b[key];rtol=1e-12,atol=1e-13)
    end
end
