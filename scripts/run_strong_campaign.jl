include("run_controlled_campaign.jl")
include("check_initial_turbulence.jl")
root=length(ARGS)>0 ? ARGS[1] : "data/strong48-v1"
arch=length(ARGS)>1 ? ARGS[2] : "GPU"
check_initial_turbulence()
cases=("strong48_f030_heated",
       "strong48_f030_cooled",
       "strong48_f100_heated",
       "strong48_f100_cooled",
       "strong48_f100_heated_h256",
       "strong48_f100_cooled_h256")
for name in cases
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_controlled_diagnostics(out)
    check_turbulent_output(out)
    GC.gc()
end
