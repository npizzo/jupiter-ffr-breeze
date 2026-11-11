using Test, TOML, JLD2, Statistics

function check_large_output(out; budget=true)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test m["status"] == "COMPLETE"
    rows=readlines(joinpath(out,"diagnostics.csv")); names=split(rows[1],',')
    a=permutedims(reduce(hcat,[parse.(Float64,split(row,',')) for row in rows[2:end]]))
    col(key)=a[:,findfirst(==(key),names)]
    @test all(isfinite,a)
    @test all(diff(col("time_seconds")) .> 0)
    @test minimum(col("cap_core_min_N2_s2")) >= m["cap_min_fraction"]*m["cap_n2"]
    @test all(col("kinetic_energy_J_m2") .>= 0)
    @test all(col("sponge_power_W_m2") .<= 1e-12)
    @test maximum(abs.(col("external_enthalpy_net_W_m2") .- m["heat_flux"] .- col("upper_enthalpy_W_m2"))) < 1e-9
    if budget && size(a,1)>2
        times=col("time_seconds"); rates=col("bottom_theta_flux_kg_K_m2_s")+col("upper_theta_rate_kg_K_m2_s")
        predicted=sum(diff(times).*(rates[1:end-1]+rates[2:end])./2)
        measured=last(col("theta_column_kg_K_m2"))-first(col("theta_column_kg_K_m2"))
        # Normalize by total absolute source throughput, not a nearly cancelling net flux.
        absolute=col("bottom_theta_flux_kg_K_m2_s")+abs.(col("upper_theta_rate_kg_K_m2_s"))
        throughput=sum(diff(times).*(absolute[1:end-1]+absolute[2:end])./2)
        residual=abs(measured-predicted)/max(throughput,1e-8)
        @test residual < 0.01
        println("Integrated rho-theta source residual / absolute throughput: ",residual)
    end
    files=sort(filter(p->startswith(basename(p),"frame_") && endswith(p,".jld2"),readdir(out;join=true)))
    times=Float64[]
    for path in files
        jldopen(path) do d
            push!(times,d["t"])
            @test size(d["vorticity_slices"]) == (m["nx"],m["ny"],3)
            @test d["video_heights"] == m["video_heights_actual_m"]
            @test all(isfinite,d["vorticity_slices"])
            @test Float32.(d["zeta"]) == d["vorticity_slices"][:,:,1]
            @test d["wall_normal_max"] < 1e-10
        end
    end
    @test all(diff(times) .≈ m["output_seconds"])
    @test last(times) ≈ m["stop_seconds"]
    @test m["simulated_seconds"] ≈ m["stop_seconds"]
    @test length(filter(f->endswith(f,".jld2"),readdir(joinpath(out,"checkpoints")))) == 1
    println("Stable-cap, budgets, wall, checkpoint and three-depth outputs passed: ",out)
end

function latest_checkpoint(out)
    files=filter(f->endswith(f,".jld2"),readdir(joinpath(out,"checkpoints");join=true))
    length(files)==1 || error("Expected one retained checkpoint: $out")
    return only(files)
end
