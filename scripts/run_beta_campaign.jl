include("run_regional.jl")
include("check_regional.jl")
cases = ("regional20_periodic_f", "regional20_channel_f", "regional20_channel_beta", "regional48_channel_f_h512", "regional48_channel_beta_h512")
index=parse(Int,ARGS[1]); index in eachindex(cases) || error("Invalid case")
root=length(ARGS)>1 ? ARGS[2] : "data/regional20-v1"
name=cases[index]; out=joinpath(root,name)
main(["configs/"*name*".toml",out,"GPU"])
check_regional_output(out)

if index in (3,5)
    other=joinpath(root,cases[index-1],"frame_00000.jld2")
    if isfile(other)
        a=load(other);b=load(joinpath(out,"frame_00000.jld2"))
        for key in ("u","v","w","theta_profile","horizontal_ke_profile","vertical_ke_profile","ubar","vbar")
            @test a[key] == b[key]
        end
        println("Matched initial fields verified against corresponding f-plane channel")
    end
end
