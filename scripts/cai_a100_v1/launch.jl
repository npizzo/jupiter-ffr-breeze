using CUDA,TOML,SHA
CUDA.functional() || error("GPU required; no CPU fallback")
@assert v"12.2" <= CUDA.runtime_version() < v"12.3"
gpu=CUDA.name(CUDA.device())
occursin("A100",gpu) || error("This relaunch requires an A100; got $gpu")
config,out=ARGS
frozen=TOML.parsefile(joinpath(@__DIR__,"frozen-source.toml"))
for (file,hash) in frozen
    bytes2hex(sha256(read(file)))==hash || error("Frozen file changed: $file")
end
receipt=TOML.parsefile(joinpath(@__DIR__,"implementation-pass.toml"))
@assert receipt["status"]=="IMPLEMENTATION_PASS_COST_CHECK_FAILED"
@assert receipt["model_sha256"]==frozen["scripts/cai_benchmark_v1/model.jl"]
include("../cai_benchmark_v1/model.jl")
cfg=TOML.parsefile(config)
@assert cfg["stop_time"]==5000 && cfg["wall_seconds"]==13800
@assert [cfg["nx"],cfg["ny"],cfg["nz"]]==[128,128,64]
# Preserve the tested model/config verbatim. Timing is recorded during this run;
# it does not cancel production. The existing wall-time stop retains a checkpoint.
println("A100_PRODUCTION_START gpu=",gpu," config=",config," reused_validation=",receipt["validation_job"]);flush(stdout)
sim,meta=run_case(cfg,out)
write_toml(joinpath(out,"relaunch-provenance.toml"),Dict("gpu"=>gpu,"validation_job"=>receipt["validation_job"],"model_sha256"=>receipt["model_sha256"],"config_sha256"=>frozen[config],"physical_changes"=>String[],"numerical_configuration_changes"=>String[],"cost_forecast_is_stop_condition"=>false,"allocated_wall_seconds"=>14400,"final_status"=>meta["status"]))
