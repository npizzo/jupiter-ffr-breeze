include("check_output.jl")
function check_native_output(out)
    check_shear_output(out);check_zonal_frames(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test m["native_budget_checks"]>0
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
