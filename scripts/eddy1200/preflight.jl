include("run.jl")
include("check_output.jl")
root=ARGS[1];stage=ARGS[2]
stage in ("thermal","probe") || error("Unknown validation stage")
TOML.parsefile(joinpath(root,"UNIT_PASS.toml"))["status"]=="PASS" || error("Initializer checks incomplete")
base=TOML.parsefile("configs/eddy20_beta_heated.toml")
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_shear_output(out);GC.gc();CUDA.reclaim()
    return out
end
if stage=="thermal"
    small=Dict("nx"=>128,"ny"=>128,"nz"=>32,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0)
    full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
    part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
    resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
    @testset "Large eddy restart" begin
        a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
            @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
        end
    end
    initial=Ref{Any}(nothing)
    for plane in ("f","beta"), heating in ("heated","control")
        cfg=TOML.parsefile("configs/eddy20_$(plane)_$(heating).toml")
        changes=merge(cfg,small,Dict("stop_seconds"=>86400.0,"output_seconds"=>3600.0,"diagnostic_seconds"=>600.0,"controller_seconds"=>600.0,"checkpoint_seconds"=>86400.0))
        out=experiment("day1_$(plane)_$(heating)",changes)
        current=load(joinpath(out,"frame_00000.jld2"))
        if initial[]!==nothing
            for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices")
                @test current[key]==initial[][key]
            end
        end
        initial[]=current
    end
    open(joinpath(root,"THERMAL_PASS.toml"),"w") do io;TOML.print(io,Dict("status"=>"PASS"));end
else
    TOML.parsefile(joinpath(root,"THERMAL_PASS.toml"))["status"]=="PASS" || error("Short runs not passed")
    out=experiment("grid1024_probe",Dict("stop_seconds"=>7200.0,"output_seconds"=>3600.0,"checkpoint_seconds"=>7200.0))
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*240
    @test hours<18
    open(joinpath(root,"PASS.toml"),"w") do io
        TOML.print(io,Dict("status"=>"PASS","estimated_10day_wall_hours"=>hours))
    end
end
