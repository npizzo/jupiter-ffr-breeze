# Check the isolated runtime, then execute an unchanged frozen v3 entry point.
include("environment.jl")
length(ARGS)>=3 || error("usage: launch.jl PREFLIGHT CPU|GPU SCRIPT [SCRIPT_ARGUMENTS...]")
preflight,architecture,script=ARGS[1:3]
architecture in ("CPU","GPU") || error("Unknown architecture")
verify_repaired_environment(abspath(preflight))
if architecture=="GPU"
    @eval using CUDA
    CUDA.functional(true) || error("CUDA unavailable in repaired environment; refusing CPU fallback")
    runtime=CUDA.runtime_version()
    (runtime.major,runtime.minor)==(12,2) || error("Runtime selection changed")
    println("Runtime gate: GPU=",CUDA.name(CUDA.device())," runtime=",runtime," driver=",CUDA.driver_version())
end
script=abspath(script)
dirname(script)==joinpath(REPAIR_REPO,"scripts","ferrel2026_v3") || error("Wrapper accepts only frozen v3 entry points")
script_args=copy(ARGS[4:end]);empty!(ARGS);append!(ARGS,script_args)
include(script)
