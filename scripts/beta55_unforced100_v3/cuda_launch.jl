# Verify the runtime repair, then enter an unchanged campaign driver.
using TOML,SHA,CUDA
pass=TOML.parsefile("data/cuda122-probe-v1/PASS.toml")
pass["status"]=="PASS" || error("CUDA/Breeze GPU probe has not passed")
for (f,h) in pass["hashes"]
 bytes2hex(sha256(read(f)))==h || error("Runtime environment changed: $f")
end
CUDA.functional() || error("GPU unavailable; refusing CPU fallback")
@assert v"12.2" <= CUDA.runtime_version() < v"12.3"
println("Verified CUDA runtime=",CUDA.runtime_version()," GPU=",CUDA.name(CUDA.device()))
script=popfirst!(ARGS)
script in ["scripts/beta55_unforced100_v3/preflight.jl","scripts/beta55_unforced100_v3/campaign.jl"] || error("Unapproved unforced55 v2 entry point")
include(abspath(script))
