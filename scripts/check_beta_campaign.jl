using Test, TOML, JLD2
include("run.jl")
root=length(ARGS)>0 ? ARGS[1] : "data/regional-check-v1"
main(["configs/regional_check_periodic.toml",joinpath(root,"legacy_periodic"),"GPU"])
include("run_regional.jl")
include("check_regional.jl")
check_initial_turbulence()
for name in ("regional_check_periodic","regional_check_channel_f","regional_check_channel_beta0","regional_check_channel_beta","regional_check_rest")
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,"GPU"])
    check_regional_output(out)
end
for (a,b) in (("legacy_periodic","regional_check_periodic"),("regional_check_channel_f","regional_check_channel_beta0"))
    fa=load(joinpath(root,a,"frame_00002.jld2")); fb=load(joinpath(root,b,"frame_00002.jld2"))
    for key in ("u","v","w","theta_profile","w_xz")
        @test isapprox(fa[key],fb[key];rtol=1e-10,atol=1e-12)
    end
end
# Verify the pinned library evaluates the intended Coriolis gradient on this grid.
grid=RectilinearGrid(size=(8,8,8),x=(0,1200000),y=(-600000,600000),z=(0,120000),topology=(Periodic,Bounded,Bounded))
beta=2.8216700143325537e-12; f=2.880957739764384e-4
c=BetaPlane(f₀=f,β=beta)
@test Oceananigans.Coriolis.fᶠᶠᵃ(1,1,1,grid,c) ≈ f-beta*600000
@test Oceananigans.Coriolis.fᶠᶠᵃ(1,9,1,grid,c) ≈ f+beta*600000
println("BETA CAMPAIGN PREFLIGHT PASSED: regression, beta-zero equivalence, walls, budgets, rest and Coriolis gradient")
