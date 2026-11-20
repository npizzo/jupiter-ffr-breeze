include("model.jl")
root=ARGS[1];mkpath(root)
c=TOML.parsefile("configs/anelastic-parent-v1/parent.toml")
c=merge(c,Dict("nx"=>8,"ny"=>8,"nz"=>32,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>120.))
s,mt=run_case(c,joinpath(root,"rest");arch=CPU(),rest=true)
a=arrays(s.model);@assert maximum(abs,a.w)<1e-10
@assert maximum(abs,a.theta.-1)<1e-12
uniform=merge(c,Dict("initial_delta_theta"=>0.,"noise_amplitude"=>0.,"cap_delta_theta"=>0.))
s,mt=run_case(uniform,joinpath(root,"flux_check");arch=CPU())
rows=[parse.(Float64,split(l,',')) for l in readlines(joinpath(root,"flux_check/diagnostics.csv"))[2:end]]
@assert maximum(abs,r[7] for r in rows)<1e-5
s,mt=run_case(c,joinpath(root,"heated");arch=CPU())
@assert mt["status"]=="TARGET_REACHED_REVIEW_REQUIRED"
write_toml(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","model_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl")))),"scope"=>"Rest, boundary heat flux, potential-temperature budget, mass-flux projection and dry integration"))
println("ANELASTIC_CPU_PREFLIGHT_PASS")
