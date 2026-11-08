using Test, TOML, JLD2
include("check_output.jl")
include("check_initial_turbulence.jl")
function check_regional_budget(root)
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

function check_regional_output(out)
    check_output(out); check_regional_budget(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    frame=load(joinpath(out,"frame_00000.jld2"))
    final=load(last(sort(filter(f->endswith(f,".jld2"),readdir(out;join=true)))))
    @test all(isfinite,final["ubar"]) && all(isfinite,final["Rxy"])
    @test final["wall_normal_max"] < 1e-10
    @test m["f_south"] ≈ m["f"]-m["beta"]*m["ly"]/2
    @test m["f_north"] ≈ m["f"]+m["beta"]*m["ly"]/2
    if m["initial_velocity"] == "turbulence"
        check_turbulent_output(out)
    else
        @test maximum(abs,final["w_xz"]) < 1e-4
    end
    println("Regional boundaries, initial state and depth-averaged outputs passed: ",out)
end
