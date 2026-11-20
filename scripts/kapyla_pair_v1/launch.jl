include("model.jl")
CUDA.functional()||error("GPU required")
occursin("A100",CUDA.name(CUDA.device()))||error("A100 required")
cfgfile,out,passfile=ARGS
receipt=TOML.parsefile(passfile)
@assert receipt["status"]=="PASS"
@assert receipt["model_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl"))))
cfg=TOML.parsefile(cfgfile)
@assert [cfg["nx"],cfg["ny"],cfg["nz"]]==[256,256,128]
@assert cfg["wall_seconds"]<=27900
run_case(cfg,out)
