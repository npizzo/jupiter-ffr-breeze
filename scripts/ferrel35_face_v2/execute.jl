using TOML, JLD2, Test, SHA, CUDA
include("run.jl")
module Reference
include("../cascade35/run.jl")
end
root=ARGS[1];mkpath(root)
checkpoint=only(filter(p->endswith(p,".jld2"),readdir("data/cascade35-300-v1/cascade35_heated/segment_30/checkpoints";join=true)))
base=TOML.parsefile("configs/ferrel35-face-v2/diagnostic.toml")
original=TOML.parsefile("configs/cascade35_heated.toml")
for (k,v) in original
    k in ["name","stop_seconds"] && continue
    @test base[k]==v
end
function trial(label,cfg;reference=false,pickup=checkpoint)
    out=joinpath(root,label);conf=joinpath(root,label*".toml")
    open(conf,"w") do io;TOML.print(io,cfg);end
    args=[conf,out,"GPU"];pickup===nothing || push!(args,pickup)
    reference ? Reference.run_eddy(args) : run_eddy(args)
    GC.gc();CUDA.reclaim()
    files=sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))
    return load(last(files))
end
# A short real mature-state comparison, not a new random state.
probe=merge(base,Dict("stop_seconds"=>300*86400.0+600,"output_seconds"=>60.0,"diagnostic_seconds"=>60.0,"checkpoint_seconds"=>600.0))
a=trial("reference",probe;reference=true)
b=trial("instrumented",probe)
@testset "Diagnostic observer does not alter day-300 evolution" begin
    for key in ["u","v","w","theta","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","zonal_u","zonal_v","zonal_w","eddy_uv","eddy_uw","native_rhov","native_rhow"]
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
    @test b["face_steps"]>0
    @test b["face_origin_time"]==300*86400.0
    @test b["face_end_time"]==300*86400.0+600
    @test maximum(abs,b["face_heat_residual"])<1e-8
end
# Zero-covariance deterministic face-product check (same collocated constants).
f=(fill(6.0,4,2),fill(8.0,3,3),fill(2.0,4,2),fill(2.0,3,3),fill(3.0,4,2),fill(4.0,3,3),fill(6.0,4,2),fill(8.0,3,3))
@test all(iszero,flux_partition(f,1.0))
open(joinpath(root,"PREFLIGHT_PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","checkpoint"=>checkpoint,"checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint)),"physical_field_comparisons"=>20,"start_day"=>300,"production_stop_day"=>301,"gpu"=>CUDA.name(CUDA.device()),"cuda"=>string(CUDA.runtime_version())))
end
println("FACE_PREFLIGHT_PASS; starting isolated day300–301 diagnostic restart");flush(stdout)
a=nothing;b=nothing;GC.gc();CUDA.reclaim()
trial("production",base)
println("FACE_PRODUCTION_COMPLETE");flush(stdout)
