include("run.jl")
include("check_output.jl")
root=only(ARGS);mkpath(root)
unit=TOML.parsefile(joinpath(root,"UNIT_PASS.toml"))
unit["status"]=="PASS" || error("Unit checks incomplete")
for (key,file) in [("initializer_sha256","initial_flow.jl"),("transfer_sha256","transfer.py")]
    unit[key]==bytes2hex(sha256(read(joinpath(@__DIR__,file)))) || error("Source changed after unit checks")
end
base=TOML.parsefile("configs/cascade100_f_heated.toml")
@test base["shear_rossby"]==0 && base["jet_speed"]==0 && base["beta"]==0
function experiment(name,changes;pickup=false)
    cfg=merge(base,Dict(changes));cfg["name"]=name
    config=joinpath(root,name*".toml")
    open(config,"w") do io;TOML.print(io,cfg);end
    out=joinpath(root,name);args=[config,out,"GPU"]
    pickup===false || push!(args,pickup)
    run_eddy(args);check_shear_output(out);GC.gc();CUDA.reclaim()
    return out
end
small=Dict("nx"=>64,"ny"=>64,"nz"=>64,"output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"controller_seconds"=>60.0,"checkpoint_seconds"=>600.0,"turbulence_wavelength"=>2400000.0)
# Wider seed only in this inexpensive restart check to satisfy >=16 cells/wavelength.
full=experiment("restart_full",merge(small,Dict("stop_seconds"=>1200.0)))
part=experiment("restart_first",merge(small,Dict("stop_seconds"=>600.0)))
resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>1200.0));pickup=latest_checkpoint(part))
@testset "Coarse campaign restart" begin
    a=load(joinpath(full,"frame_00002.jld2"));b=load(joinpath(resumed,"frame_00002.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy")
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
end
initial=nothing;estimates=Float64[]
for mode in ("heated","control")
    c=TOML.parsefile("configs/cascade100_f_$(mode).toml")
    out=experiment("day1_$(mode)",merge(c,Dict("stop_seconds"=>86400.0,"output_seconds"=>10800.0,"checkpoint_seconds"=>86400.0)))
    current=load(joinpath(out,"frame_00000.jld2"))
    if initial!==nothing
        for key in ("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar")
            @test current[key]==initial[key]
        end
    end
    global initial=current
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*240
    @test hours<4
    push!(estimates,hours)
end
# Numerical transfer tests and real-output validation run on Unity too.
run(`tools/field-export-env/bin/python scripts/cascade100/transfer.py --self-test`)
for mode in ("heated","control")
    run(`tools/field-export-env/bin/python scripts/cascade100/transfer.py $(joinpath(root,"day1_"*mode))`)
end
hashes=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if endswith(f,".jl") || endswith(f,".py"))
config_hashes=Dict("cascade100_f_"*mode=>bytes2hex(sha256(read("configs/cascade100_f_$(mode).toml"))) for mode in ("heated","control"))
open(joinpath(root,"PASS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","estimated_10day_wall_hours"=>estimates,"source_hashes"=>hashes,"configuration_hashes"=>config_hashes))
end
