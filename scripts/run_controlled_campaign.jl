include("run.jl")
include("check_output.jl")
using Test

function check_controlled_diagnostics(root)
    m=TOML.parsefile(joinpath(root,"metadata.toml"))
    lines=readlines(joinpath(root,"diagnostics.csv"))
    names=split(first(lines), ',')
    data=permutedims(reduce(hcat,[parse.(Float64,split(line,',')) for line in lines[2:end]]))
    col(name)=data[:,findfirst(==(name),names)]
    @test all(isfinite,data)
    @test all(diff(col("time_seconds")) .> 0)
    @test all(col("kinetic_energy_J_m2") .>= 0)
    @test all(col("sponge_power_W_m2") .<= 1e-12)
    if m["cooling_model"] == "fixed_power"
        @test all(col("upper_enthalpy_W_m2") .<= 1e-12)
        @test all(abs.(col("upper_enthalpy_W_m2") .+ m["upper_cooling_flux"]) .< 1e-10)
        @test all(abs.(col("external_enthalpy_net_W_m2") .- (m["heat_flux"]-m["upper_cooling_flux"])) .< 1e-10)
        dt=last(col("time_seconds"))-first(col("time_seconds"))
        expected=(first(col("bottom_theta_flux_kg_K_m2_s"))+first(col("upper_theta_rate_kg_K_m2_s")))*dt
        change=last(col("theta_column_kg_K_m2"))-first(col("theta_column_kg_K_m2"))
        @test isapprox(change,expected;rtol=1e-3,atol=1e-3)
        if startswith(m["name"],"controlled_check")
            @test change < 0 # Equal enthalpy powers need not give equal rho-theta powers.
            @test isapprox(change,expected;rtol=1e-4,atol=1e-6)
        end
    end
    println("Source and inventory checks passed: ",root)
end

root=length(ARGS)>0 ? ARGS[1] : "data/controlled48-v1"
arch=length(ARGS)>1 ? ARGS[2] : "GPU"
for name in ("rest","heating","controlled_check_sink","controlled_check_balance")
    out=joinpath(root,"check_"*name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_controlled_diagnostics(out)
    GC.gc()
end
cases=("controlled48_on_p020","controlled48_off_p020",
       "controlled48_on_p002","controlled48_off_p002",
       "controlled48_unforced","controlled48_refined")
for name in cases
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_controlled_diagnostics(out)
    GC.gc()
end
