using TOML,SHA
include("run.jl")
root,case=ARGS;p=TOML.parsefile(joinpath(root,"preflight","PASS.toml"));p["status"]=="PASS"||error("Unvalidated native budget")
for (f,h) in p["source_hashes"];bytes2hex(sha256(read(joinpath(@__DIR__,f))))==h||error("Changed frozen source $f");end
cfg="configs/fplane-budget-v1/$(case).toml";bytes2hex(sha256(read(cfg)))==p["config_hashes"][cfg]||error("Changed config")
bytes2hex(open(sha256,"data/fplane-point4-v2/preflight/day90_heated_bundle.jld2"))==p["bundle_sha256"]||error("Changed restart")
run_eddy([cfg,joinpath(root,case),"GPU"])
println("DENSE_BUDGET_CASE_COMPLETE ",case);flush(stdout)
