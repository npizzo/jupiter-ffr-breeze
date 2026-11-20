include("model.jl")
include("checks.jl")
CUDA.functional()||error("GPU required")
root=ARGS[1];mkpath(root)
cfg=TOML.parsefile("configs/chan2013-v1/case4.toml")
err=operator_checks(cfg)
small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>40,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>300.,"save_3d"=>false))
rest=merge(small,Dict("heat"=>0.,"perturbation"=>0.,"stop_time"=>1.,"diagnostic_interval"=>1.,"movie_interval"=>1.,"checkpoint_interval"=>1.))
s,mt=run_case(rest,joinpath(root,"rest"));restspeed=maximum(abs,arrays(s.model).w);@assert restspeed<1e-8
s1,mt1=run_case(small,joinpath(root,"first"))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"first/checkpoints");join=true)))
s2,mt2=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"resumed");pickup=cp)
s3,mt3=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"continuous"))
a=arrays(s2.model);b=arrays(s3.model);errs=Dict{String,Float64}()
for key in (:u,:v,:w,:rho,:T)
    e=maximum(abs,getproperty(a,key) .- getproperty(b,key));scale=maximum(abs,getproperty(b,key))
    @assert e<=1e-9+1e-6scale
    errs[string(key)]=e
end
println("GPU_REST_RESTART_PASS ",errs);flush(stdout)
# Both rotation extremes must survive heating and satisfy the explicit energy check.
for om in (0.,3.)
    pilot=merge(small,Dict("omega"=>om,"nx"=>32,"ny"=>32,"nz"=>80,"stop_time"=>20.,"diagnostic_interval"=>.25,"movie_interval"=>5.,"checkpoint_interval"=>20.,"wall_seconds"=>300.))
    s,mt=run_case(pilot,joinpath(root,"heated_omega_$(om)"));@assert mt["status"]=="COMPLETE"
end
s=nothing;s1=nothing;s2=nothing;s3=nothing;a=nothing;b=nothing;GC.gc();CUDA.reclaim()
target=merge(cfg,Dict("stop_time"=>1.,"diagnostic_interval"=>.25,"movie_interval"=>1.,"checkpoint_interval"=>1.,"wall_seconds"=>600.,"save_3d"=>false))
s,mt=run_case(target,joinpath(root,"target_grid"));@assert mt["status"]=="COMPLETE"
rows=[parse.(Float64,split(l,',')) for l in readlines(joinpath(root,"target_grid/diagnostics.csv"))[2:end]]
rate=(rows[end][4]-rows[end-1][4])/(rows[end][1]-rows[end-1][1])
write_toml(joinpath(root,"COST.toml"),Dict("seconds_per_unit"=>rate,"target_forecast_hours_with_margin"=>(1.5rate*cfg["stop_time"]+300)/3600,"first_allocation_hours_per_case"=>8.,"gpu"=>CUDA.name(CUDA.device()),"warning"=>"Early-time estimate; turbulent CFL and SGS diffusion can slow later integration"))
28800/(1.5rate)>=1000||error("Insufficient progress within eight-hour resource cap")
hashes=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ("model.jl","checks.jl","validate.jl","launch.jl"))
configs=Dict(f=>bytes2hex(sha256(read(joinpath("configs/chan2013-v1",f)))) for f in ("case0.toml","case3.toml","case4.toml"))
write_toml(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","source_hashes"=>hashes,"config_hashes"=>configs,"restart_errors"=>errs,"rest_max_w"=>restspeed,"stress_error"=>err,"grid"=>[256,256,80],"seconds_per_unit"=>rate,"scope"=>"Implementation/early-energy/timing gate only; not turbulent equilibration or replication certificate"))
println("CHAN_VALIDATION_PASS");flush(stdout)
