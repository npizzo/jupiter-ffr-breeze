include("model.jl")
CUDA.functional()||error("GPU required")
occursin("A100",CUDA.name(CUDA.device()))||error("A100 required")
cfgfile,out,passfile=ARGS
receipt=TOML.parsefile(passfile);@assert receipt["status"]=="PASS"
for (f,h) in receipt["source_hashes"]
    @assert h==bytes2hex(sha256(read(joinpath(@__DIR__,f))))
end
@assert receipt["config_hashes"][basename(cfgfile)]==bytes2hex(sha256(read(cfgfile)))
cfg=TOML.parsefile(cfgfile)
@assert [cfg["nx"],cfg["ny"],cfg["nz"]]==[256,256,80]
@assert cfg["wall_seconds"]<=27900
run_case(cfg,out)
