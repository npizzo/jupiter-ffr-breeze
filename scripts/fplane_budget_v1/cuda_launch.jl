using TOML,SHA,CUDA
pass=TOML.parsefile("data/cuda122-probe-v1/PASS.toml");pass["status"]=="PASS"||error("CUDA gate")
for (f,h) in pass["hashes"];bytes2hex(sha256(read(f)))==h||error("Runtime changed");end
CUDA.functional()||error("GPU required");@assert v"12.2"<=CUDA.runtime_version()<v"12.3"
println("Verified runtime ",CUDA.runtime_version()," GPU ",CUDA.name(CUDA.device()));flush(stdout)
mode=popfirst!(ARGS);mode in ["pilot","preflight","campaign"]||error("mode")
include(joinpath(@__DIR__,mode*".jl"))
