include("run_shear.jl")
include("check_shear.jl")
root=ARGS[1];stage=ARGS[2]
stage in ("thermal","probe") || error("Unknown stage")
stage=="thermal" && isdir(root) && !isempty(readdir(root)) && error("Nonempty root")
mkpath(root)
base=TOML.parsefile("configs/shear20_heated.toml")
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_shear(args);check_shear_output(out);GC.gc();CUDA.reclaim()
    return out
end
if stage=="thermal"
    small=Dict("nx"=>32,"ny"=>32,"nz"=>32,"lx"=>600000.0,"ly"=>600000.0,
               "output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>1800.0)
    rest=experiment("balanced_rest",merge(small,Dict("initial_velocity"=>"rest","perturbation"=>0.0,"heat_flux"=>0.0,"stop_seconds"=>86400.0,"output_seconds"=>3600.0)))
    full=experiment("restart_full",merge(small,Dict("stop_seconds"=>3600.0)))
    part=experiment("restart_first",merge(small,Dict("stop_seconds"=>1800.0)))
    resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>3600.0));pickup=latest_checkpoint(part))
    uniform=experiment("passivity_uniform",merge(small,Dict("stop_seconds"=>3600.0,"tracer_uniform"=>true)))
    @testset "Feedback and tracer restart; tracer passivity" begin
        a=load(joinpath(full,"frame_00006.jld2"));b=load(joinpath(resumed,"frame_00006.jld2"));c=load(joinpath(uniform,"frame_00006.jld2"))
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
            @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
        end
        for key in ("u","v","w","theta_profile")
            @test isapprox(a[key],c[key];rtol=1e-10,atol=1e-10)
        end
        @test maximum(abs.(c["tracer_slices"].-0.5)) < 1e-6
    end
    # Long test at the same total thermal forcing and cap depth, with a cheaper domain.
    pilot=experiment("budget40_small",Dict("nx"=>64,"ny"=>64,"lx"=>600000.0,"ly"=>600000.0,"stop_seconds"=>40*86400.0))
    open(joinpath(root,"THERMAL_PASS.toml"),"w") do io;TOML.print(io,Dict("status"=>"PASS"));end
else
    TOML.parsefile(joinpath(root,"THERMAL_PASS.toml"))["status"]=="PASS" || error("Thermal tests not passed")
    records=Dict{String,Any}()
    for n in (1024,1536)
        out=experiment("grid_probe_"*string(n),Dict("nx"=>n,"ny"=>n,"stop_seconds"=>7200.0,"output_seconds"=>3600.0,"checkpoint_seconds"=>7200.0))
        m=TOML.parsefile(joinpath(out,"metadata.toml"))
        hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*10*24
        @test hours<18 # production has a 24-hour limit
        records[string(n)*"_estimated_10day_wall_hours"]=hours
    end
    records["status"]="PASS"
    open(joinpath(root,"PASS.toml"),"w") do io;TOML.print(io,records);end
end
