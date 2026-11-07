include("run_controlled_campaign.jl")
include("check_initial_turbulence.jl")

root=length(ARGS)>0 ? ARGS[1] : "data/turbulence48-v1"
arch=length(ARGS)>1 ? ARGS[2] : "GPU"
check_initial_turbulence()
for name in ("turbulence_check","turbulence48_heated","turbulence48_cooled","turbulence48_unforced")
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_controlled_diagnostics(out)
    check_turbulent_output(out)
    GC.gc()
end
