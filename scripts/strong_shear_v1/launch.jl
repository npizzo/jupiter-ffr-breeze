using TOML, SHA, CUDA
started_job=time()
pass=TOML.parsefile("data/cuda122-probe-v1/PASS.toml")
@assert pass["status"]=="PASS"
for (file,h) in pass["hashes"]
    bytes2hex(sha256(read(file)))==h || error("CUDA runtime changed: $file")
end
CUDA.functional() || error("GPU required")
@assert v"12.2" <= CUDA.runtime_version() < v"12.3"
for (file,h) in TOML.parsefile(joinpath(@__DIR__,"frozen-source.toml"))
    bytes2hex(sha256(read(String(file))))==String(h) || error("Frozen source changed")
end
include("../shear_layer_v1/run.jl")

root,case=ARGS
@assert case in ("stable_cyclonic","stable_anticyclonic")
ispath(root) && error("Refuse duplicate experiment output")
original=TOML.parsefile("configs/strong-seed-v1/stable_none.toml")
cfg=TOML.parsefile("configs/strong-shear-v1/$case.toml")
@assert Set(keys(original))==Set(keys(cfg))
allowed=Set(["name","shear_rossby"])
changes=sort([k for k in keys(cfg) if cfg[k]!=original[k]])
Set(changes)==allowed || error("Unannounced parameter change")
@assert cfg["turbulence_rossby"]==0.1 && cfg["stop_seconds"]==50*86400
@assert cfg["shear_rossby"]==(case=="stable_cyclonic" ? 0.05 : -0.05) && cfg["beta"]==0
mkpath(root)
function save_toml(path,data)
    open(path,"w") do io; TOML.print(io,data); end
end
save_toml(joinpath(root,"experiment.toml"),Dict(
    "physical_changes"=>["Add initial sinusoidal shear, peak zeta/f=$(cfg["shear_rossby"]), to fixed-Ro large-box control"],
    "changed_config_keys"=>changes,"case_count"=>1,"target_days"=>50,
    "gpu_hours_hard_limit"=>3,"initial_gate_days"=>1,
    "reference"=>"data/strong-seed-v1-1983355",
    "interpretation"=>"Coarse original-scale control; 128 squared versus original 256 squared and current odd-mode seed retained"))
pilot=copy(cfg);pilot["stop_seconds"]=86400.;pilot["wall_time_limit"]=900.
pilotconfig=joinpath(root,"pilot.toml");save_toml(pilotconfig,pilot)
println("START_ONE_DAY_GATE ",CUDA.name(CUDA.device()));flush(stdout)
run_eddy([pilotconfig,joinpath(root,"day0_1"),"GPU"])
m=TOML.parsefile(joinpath(root,"day0_1","metadata.toml"))
@assert m["status"]=="COMPLETE" && isapprox(m["simulated_seconds"],86400.;atol=1e-6)
@assert 5.4 < m["initial_horizontal_rms_lower_m_s"] < 5.6
@assert m["segment_start_seconds"]==0 && m["initial_mass_divergence_relative"]<1e-12
forecast=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*49*86400
remaining=10200-(time()-started_job)
isfinite(forecast) && 1.25forecast < remaining || error("Runtime gate failed; no continuation")
checkpoints=filter(f->endswith(f,".jld2"),readdir(joinpath(root,"day0_1","checkpoints");join=true))
checkpoint=only(checkpoints)
save_toml(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","gate_days"=>1,
    "forecast_remaining_wall_seconds"=>forecast,"remaining_budget_seconds"=>remaining,
    "initial_rms_m_s"=>m["initial_horizontal_rms_lower_m_s"],
    "checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint))))
println("ONE_DAY_GATE_PASS remaining_forecast_seconds=",forecast);flush(stdout)
GC.gc();CUDA.reclaim()
cfg["wall_time_limit"]=remaining
prodconfig=joinpath(root,"continuation.toml");save_toml(prodconfig,cfg)
run_eddy([prodconfig,joinpath(root,"day1_50"),"GPU",checkpoint])
final=TOML.parsefile(joinpath(root,"day1_50","metadata.toml"))
@assert isapprox(final["segment_start_seconds"],86400.;atol=1e-6)
@assert isapprox(final["simulated_seconds"],50*86400.;atol=1e-6)
save_toml(joinpath(root,"COMPLETE.toml"),Dict("status"=>"COMPLETE","days"=>50,
    "total_wall_seconds"=>time()-started_job,"segments"=>["day0_1","day1_50"]))
println("STRONG_SHEAR_FIFTY_DAYS_COMPLETE")
