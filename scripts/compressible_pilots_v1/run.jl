using SHA, TOML
family,out=ARGS
specs=Dict(
 "kapyla"=>("kapyla_pair_v1","configs/kapyla-pair-v1/A1_target.toml","91b682d967a6127be5a7850c9b5fc22ebd98eca679ba02db6e2962404cdfec60"),
 "chan"=>("chan2013_v1","configs/chan2013-v1/case4.toml","eeb2d07fe1556a76405ca2293e25f6bd53f7d4d476d155983d7025903b9457be"))
dir,configfile,expected=specs[family]
modelpath=joinpath(@__DIR__,"..",dir,"model.jl")
@assert bytes2hex(sha256(read(modelpath)))==expected
include(modelpath)
cfg=TOML.parsefile(configfile)
merge!(cfg,Dict("name"=>family*"_coarse_feasibility","nx"=>64,"ny"=>64,"stop_time"=>500.,"wall_seconds"=>1200.,"diagnostic_interval"=>1.,"movie_interval"=>5.,"checkpoint_interval"=>25.,"field_interval"=>500.,"save_3d"=>false))
# Full box and vertical grid; no accelerated heating, stronger seed, or changed rotation.
sim,meta=run_case(cfg,out;arch=CPU())
a=arrays(sim.model);p=parameters(cfg);z=[zcell(k,p) for k=1:p.nz]
cp=family=="kapyla" ? 1. : 2.5
tbar=mean(a.T;dims=(1,2));rhobar=vec(mean(a.rho;dims=(1,2)))
hflux=vec(mean(cp .* a.rho .* a.w .* (a.T .- tbar);dims=(1,2)))
wrms=vec(sqrt.(mean(a.w .^ 2;dims=(1,2))))
write_toml(joinpath(out,"pilot_profiles.toml"),Dict("time"=>Float64(sim.model.clock.time),"z"=>z,"rho_mean"=>rhobar,"temperature_mean"=>vec(tbar),"w_rms"=>wrms,"enthalpy_flux_fluctuation"=>hflux,"mass_flux"=>vec(mean(a.rho .* a.w;dims=(1,2))),"bottom_flux"=>p.flux,"model_sha256"=>expected,"review_required"=>true,"interpretation"=>"Feasibility only; full thermal equilibration and LSV sign selection not established"))
println("PILOT_FINISHED_REVIEW_REQUIRED family=",family," time=",sim.model.clock.time)
