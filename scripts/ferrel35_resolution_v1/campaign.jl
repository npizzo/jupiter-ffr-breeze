using TOML,SHA
include("run.jl")
root,case=ARGS;case in ["coarse","fine"] || error("Unknown resolution")
gate=TOML.parsefile(joinpath(root,"preflight","PASS.toml"));gate["status"]=="PASS" || error("Preflight not passed")
for (f,h) in gate["source_hashes"]
 bytes2hex(sha256(read(joinpath(@__DIR__,f))))==h || error("Frozen source changed: $f")
end
bundle=joinpath(root,"preflight","day300_bundle.jld2")
bytes2hex(open(sha256,bundle))==gate["bundle_sha256"] || error("Restart bundle changed")
cfg=TOML.parsefile("configs/ferrel35-resolution-v1/"*case*".toml");cfg["restart_bundle"]=bundle
config=joinpath(root,case*".toml");open(config,"w") do io;TOML.print(io,cfg);end
run_eddy([config,joinpath(root,case),"GPU"])
println("RESOLUTION_CASE_COMPLETE ",case);flush(stdout)
