using JLD2, TOML, Test
root=ARGS[1]
m=TOML.parsefile(joinpath(root,"metadata.toml"))
files=sort(filter(f->endswith(f,".jld2"),readdir(root;join=true)))
@test m["status"] == "COMPLETE"
@test m["simulated_seconds"] >= m["stop_seconds"]
@test length(files) >= 2
initial=load(first(files)); final=load(last(files))
@test all(isfinite,final["theta"])
@test all(>(0),final["pressure"])
@test all(>(0),final["density"])
@test final["pressure"][end] < final["pressure"][1]
@test final["t"] >= initial["t"]
if m["name"] == "rest"
    @test maximum(abs,final["w_xz"]) < 1e-5
elseif m["name"] == "heating"
    dz=m["height"]/m["nz"]
    change=sum(final["density"] .* (final["theta_profile"]-initial["theta_profile"])) * dz
    exner=(m["base_pressure"]/1e5)^((8.314462618/m["molar_mass"])/m["heat_capacity"])
    expected=m["heat_flux"]/(m["heat_capacity"]*exner)*(final["t"]-initial["t"])
    @test change > 0
    @test isapprox(change,expected;rtol=0.01)
else
    @test final["w_rms"] > 0
end
println("Output checks passed: ",root)
