# Reuse completed matched initializations and checkpoint runs; do not rerun them.
using Test
include("run.jl")
root=only(ARGS)
@testset "Restart maximum pointwise error audit" begin
 a=load(joinpath(root,"stable_anticyclonic","frame_00002.jld2"))
 b=load(joinpath(root,"restart_second","frame_00002.jld2"))
 for k in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
  err=maximum(abs,a[k].-b[k]);scale=maximum(abs,a[k])
  println("RESTART_ERROR ",k," max_abs=",err," max_relative=",err/max(scale,1e-30))
  @test err <= 1e-10+1e-9*scale
 end
end

function experiment(name,cfg;pickup=false)
 p=joinpath(root,name*".toml");open(p,"w") do io;TOML.print(io,cfg);end
 o=joinpath(root,name);args=[p,o,"GPU"];pickup===false||push!(args,pickup)
 run_eddy(args);GC.gc();CUDA.reclaim();return o
end
function check_output(o)
 m=TOML.parsefile(joinpath(o,"metadata.toml"));@test m["status"]=="COMPLETE"
 files=sort(filter(x->startswith(basename(x),"frame_")&&endswith(x,".jld2"),readdir(o;join=true)))
 a=load(first(files));b=load(last(files))
 @test abs(b["tracer_column"]/a["tracer_column"]-1)<1e-8
 @test b["tracer_min"] > -1e-8 && b["tracer_max"] < 1+1e-8
 @test abs(m["heat_flux"]+b["relaxation_enthalpy"])<1e-6
 @test b["wall_normal_max"]<1e-10
 @test all(isfinite,b["vorticity_slices"])
 @test size(b["vorticity_slices"])==(128,128,3)
 @test all(isfinite,b["layer_data"]["signed_skew_prime_profile"])
 if m["cap_n2"]>0
  @test b["cap_min"]>=m["cap_min_fraction"]*m["cap_n2"]
 end
 return a,b,m
end
cases=["$(layer)_$(shear)" for layer in ("stable","neutral") for shear in ("cyclonic","none","anticyclonic")]
@testset "Existing six-case output recheck" begin
 for case in cases
  a,b,m=check_output(joinpath(root,case))
  for label in ("lower","upper","full")
   @test abs(a["layer_data"]["layers"][label]["signed_skew_prime"])<1e-10
  end
 end
end

cfg=TOML.parsefile("configs/shear-layer-v1/stable_anticyclonic.toml")
merge!(cfg,Dict("stop_seconds"=>21600.,"output_seconds"=>3600.,"diagnostic_seconds"=>1800.,"checkpoint_seconds"=>21600.,"wall_time_limit"=>1200.))
o=experiment("throughput",cfg);a,b,m=check_output(o)
forecast=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*864000.
@test forecast < 1800. # leave substantial margin within the 3000 s internal stop
hashes=Dict(p=>bytes2hex(sha256(read(p))) for folder in ("scripts/shear_layer_v1","configs/shear-layer-v1") for p in readdir(folder;join=true) if isfile(p))
open(joinpath(root,"PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","forecast_10day_seconds"=>forecast,"source_hashes"=>hashes,"actual_grid"=>[128,128,64],"gpu"=>m["gpu_name"],"cap_neutral_control"=>true,"restart_checked"=>true))
end
println("SHEAR_LAYER_PREFLIGHT_PASS forecast_seconds=",forecast);flush(stdout)
