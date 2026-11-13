using Test, TOML, JLD2, Statistics
include("check_large.jl")
function check_shear_output(out)
    check_large_output(out;budget=false)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    rows=readlines(joinpath(out,"diagnostics.csv"));names=split(rows[1],',')
    a=permutedims(reduce(hcat,[parse.(Float64,split(row,',')) for row in rows[2:end]]))
    col(key)=a[:,findfirst(==(key),names)]
    @test maximum(abs.(col("tracer_column_kg_m2")./first(col("tracer_column_kg_m2")).-1)) < 1e-8
    @test minimum(col("tracer_min")) > -1e-9
    @test maximum(col("tracer_max")) < 1+1e-9
    @test m["installed_external_forcing"] == "none"
    @test !m["balance_heat"]
    for key in ("bottom_theta_flux_kg_K_m2_s", "upper_theta_rate_kg_K_m2_s", "external_enthalpy_net_W_m2", "sponge_power_W_m2", "cap_enthalpy_W_m2", "budget_alpha_s_inv", "integrated_theta_source_kg_K_m2", "integrated_external_enthalpy_J_m2")
        @test all(iszero, col(key))
    end
    @test maximum(abs.(col("theta_column_kg_K_m2")./first(col("theta_column_kg_K_m2")).-1)) < 1e-8
    if m["balance_heat"]
        @test maximum(abs,col("external_enthalpy_net_W_m2")) < 1e-7
        @test last(col("max_preupdate_net_W_m2")) < 0.005*162
        dt=last(col("time_seconds"))-first(col("time_seconds"))
        imbalance=last(col("integrated_external_enthalpy_J_m2"))-first(col("integrated_external_enthalpy_J_m2"))
        @test abs(imbalance)/max(dt,1) < 0.001*162
        predicted=last(col("integrated_theta_source_kg_K_m2"))-first(col("integrated_theta_source_kg_K_m2"))
        measured=last(col("theta_column_kg_K_m2"))-first(col("theta_column_kg_K_m2"))
        @test abs(measured-predicted) < max(0.01,0.001*max(abs(predicted),dt*abs(first(col("bottom_theta_flux_kg_K_m2_s")))))
    end
    if m["initial_velocity"]=="rest"
        @test maximum(col("w_rms_center_m_s")) < 1e-3
    end
    for path in sort(filter(p->startswith(basename(p),"frame_") && endswith(p,".jld2"),readdir(out;join=true)))
        jldopen(path) do d
            @test size(d["tracer_slices"]) == (m["nx"],m["ny"],3)
            @test all(isfinite,d["tracer_slices"])
        end
    end
    println("Unforced sources, tracer and potential-temperature conservation and three-depth fields passed: ",out)
end
