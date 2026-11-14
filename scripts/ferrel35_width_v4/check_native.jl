include("check_output.jl")
function check_native_output(out)
    check_shear_output(out);check_zonal_frames(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test m["native_budget_checks"]>0
    lines=readlines(joinpath(out,"native_budget","checks.csv"))
    rows=[parse.(Float64,split(l,',')) for l in lines[2:end]]
    @test length(rows)==m["native_budget_checks"]
    @test rows[1][2]==m["segment_start_seconds"]
    @test rows[end][3]==m["stop_seconds"]
    @test all(rows[i][3]==rows[i+1][2] for i in 1:length(rows)-1)
    @test isapprox(sum(r[4] for r in rows),m["stop_seconds"]-m["segment_start_seconds"];atol=1e-7,rtol=1e-12)
    @test m["native_budget_cumulative_relative_residual"]<1e-5
    frames=sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))
    lastd=load(last(frames))
    @test lastd["native_budget_origin_time"]==m["segment_start_seconds"]
    @test lastd["native_budget_end_time"]==m["stop_seconds"]
    @test lastd["native_budget_steps"]==m["native_budget_checks"]
    for p in frames
        d=load(p)
        for k in ["native_rhou_rhs","native_integrated_parts","native_integrated_residual","centered_native_advection_remainder"]
            @test all(isfinite,d[k])
        end
        @test size(d["native_integrated_parts"])==(m["ny"],m["nz"],4)
        @test d["native_integrated_relative_residual"]<1e-5
        @test d["mass_divergence_max"]<1e-12
    end
    println("Every-step and cumulative native budget checks passed: ",out)
end
