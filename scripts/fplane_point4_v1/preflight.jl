using TOML,SHA,Test
include("run.jl")
root=only(ARGS);out=joinpath(root,"preflight");mkpath(out)
cases=TOML.parsefile(joinpath(@__DIR__,"cases.toml"))
# Every distinct pathway: old checkpoints, reversal, high flux, negative f, and 512².
names=["late_heated","late_control","late_restart_reversed","rest_q540","rest_q162_flipf","rest_q162_fine","late_heated_dt15","late_heated_fine"]
# Reuse the divergence-commuting conservative prolongation already validated in
# the preceding beta-plane audit, and revalidate it on this actual f-plane state.
exportcfg=TOML.parsefile(cases["late_heated"]["config"])
exportcfg["export_bundle"]=joinpath(out,"day90_heated_bundle.jld2")
exportpath=joinpath(out,"export.toml");open(exportpath,"w") do io;TOML.print(io,exportcfg);end
run_eddy([exportpath,joinpath(out,"export"),"GPU",cases["late_heated"]["checkpoint"]]);GC.gc();CUDA.reclaim()
for name in names
 cfg=TOML.parsefile(cases[name]["config"])
 start=startswith(name,"late_restart") ? 100*86400.0 : startswith(name,"late_") ? 90*86400.0 : 0.0
 duration=name=="late_heated_fine" ? 3600.0 : 600.0
 cfg["stop_seconds"]=start+duration;cfg["output_seconds"]=600.0;cfg["diagnostic_seconds"]=600.0;cfg["full_field_seconds"]=600.0
 config=joinpath(out,name*".toml");open(config,"w") do io;TOML.print(io,cfg);end
 args=[config,joinpath(out,name),"GPU"]
 cp=cases[name]["checkpoint"];isempty(cp) || push!(args,cp)
 run_eddy(args)
 m=TOML.parsefile(joinpath(out,name,"metadata.toml"));@test m["status"]=="COMPLETE";@test m["simulated_seconds"]==start+duration
 files=filter(f->endswith(f,".jld2"),readdir(joinpath(out,name,"fields3d");join=true));@test length(files)==Int(duration/600)+1
 for file in files
  b=load(file);@test all(all(isfinite,b[k]) for k in ["u_face","v_face","w_face","theta","dye_density"])
  @test maximum(abs,b["v_face"][:,[1,end],:])<1e-11
  @test maximum(abs,b["w_face"][:,:,[1,end]])<1e-11
 end
 println("PREFLIGHT_CASE_PASS ",name);flush(stdout)
 GC.gc();CUDA.reclaim()
end
coarse=load(joinpath(out,"late_heated_dt15","frame_12960.jld2"));fine=load(joinpath(out,"late_heated_fine","frame_12960.jld2"))
ke_fraction=fine["kinetic_energy"]/coarse["kinetic_energy"]-1
@test abs(ke_fraction)<0.1
@test isapprox(fine["theta_profile"],coarse["theta_profile"];rtol=1e-12)
@test isapprox(fine["tracer_column"],coarse["tracer_column"];rtol=1e-12)
fine_meta=TOML.parsefile(joinpath(out,"late_heated_fine","metadata.toml"))
fine_ten_day_hours=fine_meta["post_warmup_wall_seconds"]/fine_meta["post_warmup_simulated_seconds"]*240
@test fine_ten_day_hours<7
hashes=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in readdir(@__DIR__) if isfile(joinpath(@__DIR__,f)))
# GPU simulation gate; the CPU analytic diagnostics gate is added by the batch wrapper.
open(joinpath(out,"GPU_PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","compute_location"=>"Unity","gpu"=>CUDA.name(CUDA.device()),"source_hashes"=>hashes,"bundle_sha256"=>bytes2hex(open(sha256,joinpath(out,"day90_heated_bundle.jld2"))),"remap_initial_ke_fraction"=>ke_fraction,"fine_ten_day_estimated_gpu_hours"=>fine_ten_day_hours,"config_hashes"=>Dict(c["config"]=>bytes2hex(sha256(read(c["config"]))) for c in values(cases)),"cases"=>names))
end
