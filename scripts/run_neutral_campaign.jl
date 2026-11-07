include("run_controlled_campaign.jl")

root=length(ARGS)>0 ? ARGS[1] : "data/neutral48-v1"
arch=length(ARGS)>1 ? ARGS[2] : "GPU"
cases=("neutral48_heated","neutral48_cooled","neutral48_unforced","neutral48_refined")
for name in cases
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_controlled_diagnostics(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test m["convective_n2"] == 0
    @test m["beta"] == 0
    @test m["jet_speed"] == 0
    initial=load(joinpath(out,"frame_00000.jld2"))
    z=initial["z"]
    @test maximum(abs.(initial["theta_profile"][z .< m["cap_height"]] .- m["theta0"])) < 1e-9
    println("Neutral initial-profile check passed: ",out)
    GC.gc()
end
