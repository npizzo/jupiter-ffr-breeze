using CUDA,TOML,SHA
CUDA.functional()||error("GPU required")
@assert v"12.2"<=CUDA.runtime_version()<v"12.3"
config,out,validation=ARGS
pass=TOML.parsefile(joinpath(validation,"PASS.toml"))
@assert pass["status"]=="PASS"
@assert pass["model_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl"))))
include("model.jl")
cfg=TOML.parsefile(config)
run_case(cfg,out)
