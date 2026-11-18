using Test
include("run.jl")
root=only(ARGS);mkpath(root)
cases=["$(layer)_$(shear)" for layer in ("stable","neutral") for shear in ("cyclonic","none","anticyclonic")]
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
initials=Dict(); configurations=Dict(); reference=nothing
@testset "Matched factorial initialization and output" begin
 for case in cases
  cfg=TOML.parsefile("configs/shear-layer-v1/$case.toml")
  merge!(cfg,Dict("stop_seconds"=>1200.,"output_seconds"=>600.,"diagnostic_seconds"=>600.,"controller_seconds"=>600.,"checkpoint_seconds"=>600.,"full_field_seconds"=>600.,"wall_time_limit"=>600.))
  configurations[case]=cfg
  o=experiment(case,cfg);a,b,m=check_output(o);initials[case]=a
  for label in ("lower","upper","full")
   @test abs(a["layer_data"]["layers"][label]["signed_skew_prime"])<1e-10
  end
  @test m["initial_mass_divergence_relative"]<1e-12
  @test m["initial_random_geostrophic_residual_f0"]<1e-12
  # Remove the imposed zonal flow: all six have identical random velocities.
  up=a["u"].-mean(a["u"];dims=1)
  if reference===nothing
   global reference=(copy(up),copy(a["v"]),copy(a["tracer_slices"]))
  else
   @test isapprox(up,reference[1];atol=1e-12,rtol=1e-12)
   @test isapprox(a["v"],reference[2];atol=1e-12,rtol=1e-12)
   @test isapprox(a["tracer_slices"],reference[3];atol=1e-7,rtol=1e-7)
  end
  peak=a["layer_data"]["layers"]["full"]["zonal_vorticity_over_f"]
  if cfg["shear_rossby"]!=0
   @test sign(peak[64])==sign(cfg["shear_rossby"])
   @test isapprox(peak[64],cfg["shear_rossby"];rtol=.002)
  end
 end
 @test initials["stable_none"]["theta_profile"] != initials["neutral_none"]["theta_profile"]
end
@testset "Checkpoint equivalence" begin
 cfg=copy(configurations["stable_anticyclonic"]);cfg["stop_seconds"]=600.
 o=experiment("restart_first",cfg)
 checkpoints=sort(filter(p->endswith(p,".jld2"),readdir(joinpath(o,"checkpoints");join=true)))
 cfg["stop_seconds"]=1200.
 o2=experiment("restart_second",cfg;pickup=last(checkpoints));a,b,m=check_output(o2)
 expected=load(joinpath(root,"stable_anticyclonic","frame_00002.jld2"))
 for k in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
  @test isapprox(b[k],expected[k];atol=1e-10,rtol=1e-9)
 end
end
# Throughput includes output, controller and all new diagnostic costs on the actual grid.
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
