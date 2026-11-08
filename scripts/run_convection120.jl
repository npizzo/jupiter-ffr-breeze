include("run_controlled_campaign.jl")
include("check_initial_turbulence.jl")

const convection120_cases = (
    "convection120_f030_heated_h256_z64",
    "convection120_f030_cooled_h256_z64",
    "convection120_f030_heated_h256_z128",
    "convection120_f030_cooled_h256_z128",
)
length(ARGS) >= 1 || error("Usage: run_convection120.jl CASE_INDEX [OUTPUT_ROOT] [CPU|GPU]")
index = parse(Int, ARGS[1])
index in eachindex(convection120_cases) || error("CASE_INDEX must be 1–4")
root = length(ARGS) >= 2 ? ARGS[2] : "data/convection120-v1"
arch = length(ARGS) >= 3 ? ARGS[3] : "GPU"
name = convection120_cases[index]
out = joinpath(root, name)
check_initial_turbulence()
main(["configs/" * name * ".toml", out, arch])
check_output(out)
check_controlled_diagnostics(out)
check_turbulent_output(out)
