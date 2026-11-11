# Reuse the immutable repaired-runtime verifier; do not broaden its old launcher.
include("../ferrel2026_gpu/environment.jl")
haskey(ENV,"SLURM_JOB_ID") || error("Seed-test launch requires a Unity allocation")
length(ARGS)>=3 || error("usage: launch.jl V3_PREFLIGHT_ROOT CPU|GPU SCRIPT [ARGUMENTS...]")
prior,architecture,script=ARGS[1:3]
architecture in ("CPU","GPU") || error("Unknown architecture")
normpath(Base.active_project())==joinpath(REPAIR_ENV,"Project.toml") || error("Wrong Julia project")
verify_repaired_environment(abspath(prior))
script=abspath(script)
dirname(script)==(@__DIR__) || error("Only seedtest_v1 entry points are accepted")
expected=architecture=="CPU" ? "validate_initialization.jl" : "campaign.jl"
basename(script)==expected || error("Unexpected entry point for architecture")
if architecture=="GPU"
    @eval using CUDA
    CUDA.functional(true) || error("CUDA unavailable; refusing CPU fallback")
    runtime=CUDA.runtime_version()
    (runtime.major,runtime.minor)==(12,2) || error("Repaired CUDA runtime changed")
    println("Seed-test runtime gate: GPU=",CUDA.name(CUDA.device())," runtime=",runtime," driver=",CUDA.driver_version())
end
script_args=copy(ARGS[4:end]);empty!(ARGS);append!(ARGS,script_args)
include(script)
