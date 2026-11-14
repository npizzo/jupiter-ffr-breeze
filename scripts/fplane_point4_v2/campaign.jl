using TOML,SHA
include("run.jl")
root,case=ARGS
pass=TOML.parsefile(joinpath(root,"preflight","PASS.toml"))
pass["status"]=="PASS" || error("Preflight did not pass")
for (f,h) in pass["source_hashes"]
 bytes2hex(sha256(read(joinpath(@__DIR__,f))))==h || error("Frozen source changed: $f")
end
cases=TOML.parsefile(joinpath(@__DIR__,"cases.toml"));entry=cases[case]
config=String(entry["config"]);bytes2hex(sha256(read(config)))==pass["config_hashes"][config] || error("Configuration changed after validation")
checkpoint=String(entry["checkpoint"])
cfg=TOML.parsefile(config)
if haskey(cfg,"restart_bundle")
 bytes2hex(open(sha256,cfg["restart_bundle"]))==pass["bundle_sha256"] || error("Restart bundle changed")
end
args=[config,joinpath(root,case),"GPU"]
isempty(checkpoint) || push!(args,checkpoint)
run_eddy(args)
println("POINT4_CASE_COMPLETE ",case);flush(stdout)
