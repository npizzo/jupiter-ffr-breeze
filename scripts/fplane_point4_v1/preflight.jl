using TOML,SHA,Test
include("run.jl")
root=only(ARGS);out=joinpath(root,"preflight");mkpath(out)
cases=TOML.parsefile(joinpath(@__DIR__,"cases.toml"))
# Every distinct pathway: old checkpoints, reversal, high flux, negative f, and 512².
names=["late_heated","late_control","late_restart_reversed","rest_q540","rest_q162_flipf","rest_q162_fine"]
for name in names
 cfg=TOML.parsefile(cases[name]["config"])
 start=startswith(name,"late_restart") ? 100*86400.0 : startswith(name,"late_") ? 90*86400.0 : 0.0
 cfg["stop_seconds"]=start+600.0;cfg["output_seconds"]=600.0;cfg["diagnostic_seconds"]=600.0;cfg["full_field_seconds"]=600.0
 config=joinpath(out,name*".toml");open(config,"w") do io;TOML.print(io,cfg);end
 args=[config,joinpath(out,name),"GPU"]
 cp=cases[name]["checkpoint"];isempty(cp) || push!(args,cp)
 run_eddy(args)
 m=TOML.parsefile(joinpath(out,name,"metadata.toml"));@test m["status"]=="COMPLETE";@test m["simulated_seconds"]==start+600
 files=filter(f->endswith(f,".jld2"),readdir(joinpath(out,name,"fields3d");join=true));@test length(files)==2
 for file in files
  b=load(file);@test all(all(isfinite,b[k]) for k in ["u_face","v_face","w_face","theta","dye_density"])
  @test maximum(abs,b["v_face"][:,[1,end],:])<1e-11
  @test maximum(abs,b["w_face"][:,:,[1,end]])<1e-11
 end
 println("PREFLIGHT_CASE_PASS ",name);flush(stdout)
 GC.gc();CUDA.reclaim()
end
hashes=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in readdir(@__DIR__) if isfile(joinpath(@__DIR__,f)))
# GPU simulation gate; the CPU analytic diagnostics gate is added by the batch wrapper.
open(joinpath(out,"GPU_PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","compute_location"=>"Unity","gpu"=>CUDA.name(CUDA.device()),"source_hashes"=>hashes,"config_hashes"=>Dict(c["config"]=>bytes2hex(sha256(read(c["config"]))) for c in values(cases)),"cases"=>names))
end
