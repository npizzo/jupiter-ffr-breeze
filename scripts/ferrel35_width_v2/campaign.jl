include("common.jl")
include("run.jl")
include("check_native.jl")
length(ARGS)==4 || error("usage campaign.jl diagnostic|wide SEGMENT ROOT PREFLIGHT")
mode=ARGS[1];segment=parse(Int,ARGS[2]);root=ARGS[3];preflight=ARGS[4]
pass=verify_gate(joinpath(preflight,"PASS.toml"))
mode in ["diagnostic","wide"] || error("Unknown campaign")
cfg=TOML.parsefile(joinpath(CONFIG_ROOT,mode*".toml"))
checkpoint=false
if mode=="diagnostic"
    segment==1 || error("Only day200-210 diagnostic branch authorized")
    previous=ORIGINAL_DAY200;checkpoint=latest_checkpoint(previous)
    bytes2hex(open(sha256,checkpoint))==pass["mature_checkpoint_sha256"] || error("Day200 checkpoint changed")
else
    1<=segment<=20 || error("Only 200 wide-channel days authorized")
    verify_gate(joinpath(root,"diagnostic","DIAGNOSTIC_PASS.toml"))
    cfg["stop_seconds"]=segment*10*86400.0
    if segment>1
        previous=joinpath(root,"wide","segment_$(segment-1)")
        old=TOML.parsefile(joinpath(previous,"metadata.toml"))
        old["status"]=="COMPLETE" && old["simulated_seconds"]==(segment-1)*10*86400.0 || error("Previous wide segment incomplete")
        checkpoint=latest_checkpoint(previous)
    end
end
out=joinpath(root,mode,"segment_$segment")
config=writeconfig(joinpath(root,"configs"),mode*"_segment_$segment",cfg)
args=[config,out,"GPU"];checkpoint===false || push!(args,checkpoint)
run_eddy(args);check_native_output(out)
if mode=="diagnostic"
    run(`tools/gallery-env/bin/python $(joinpath(@__DIR__,"review_budget.py")) $out`)
    gate(joinpath(root,"diagnostic","DIAGNOSTIC_PASS.toml");scope="Numerical diagnostic gate, not certification of Ferrel cells",completed_earth_day=210)
end
println("Completed authorized ",mode," segment ",segment)
