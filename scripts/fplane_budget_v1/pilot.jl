using Test,TOML
include("run.jl")
root=only(ARGS);mkpath(root)
for audit in [false,true]
 cfg=TOML.parsefile("configs/fplane-point4-v2/early_replay_heated.toml")
 cfg["nx"]=64;cfg["ny"]=64;cfg["turbulence_wavelength"]=2400000.;cfg["stop_seconds"]=120.0;cfg["output_seconds"]=60.;cfg["diagnostic_seconds"]=60.;cfg["full_field_seconds"]=60.;cfg["max_dt"]=10.;cfg["initial_dt"]=10.;cfg["native_budget"]=audit;cfg["reduce_exports"]=false;cfg["name"]=audit ? "audit" : "baseline"
 path=joinpath(root,cfg["name"]*".toml");open(path,"w") do io;TOML.print(io,cfg);end
 run_eddy([path,joinpath(root,cfg["name"]),"GPU"]);GC.gc();CUDA.reclaim()
end
files=sort(readdir(joinpath(root,"baseline","fields3d");join=true));lastbase=load(last(files));lastaudit=load(joinpath(root,"audit","fields3d",basename(last(files))))
for k in ["u_face","v_face","w_face","theta","dye_density"]
 @test lastbase[k]==lastaudit[k]
end
open(joinpath(root,"PASS.toml"),"w") do io;TOML.print(io,Dict("status"=>"PASS","bitwise_physical_state_unchanged"=>true,"gpu"=>CUDA.name(CUDA.device())));end
println("NATIVE_VORTICITY_SMALL_GPU_PASS");flush(stdout)
