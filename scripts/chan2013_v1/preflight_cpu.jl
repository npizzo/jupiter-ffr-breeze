include("model.jl")
include("checks.jl")
mkpath(ARGS[1]);cfg=TOML.parsefile("configs/chan2013-v1/case4.toml")
operator_checks(cfg)
small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>40,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>300.,"save_3d"=>false))
rest=merge(small,Dict("heat"=>0.,"perturbation"=>0.))
s,mt=run_case(rest,joinpath(ARGS[1],"rest");arch=CPU());@assert maximum(abs,arrays(s.model).w)<1e-8
s,mt=run_case(small,joinpath(ARGS[1],"heated");arch=CPU());@assert mt["status"]=="COMPLETE"
write_toml(joinpath(ARGS[1],"PASS.toml"),Dict("status"=>"PASS","scope"=>"CPU operator and short integration checks"))
println("CHAN_CPU_PREFLIGHT_PASS")
