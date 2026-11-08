include("run_large.jl")
include("check_large.jl")
include("initial_channel_turbulence.jl")
root=length(ARGS)>0 ? ARGS[1] : "data/large40-preflight-v1"
stage=length(ARGS)>1 ? ARGS[2] : "thermal"
stage in ("thermal","probe") || error("Unknown preflight stage")
stage=="thermal" && isdir(root) && !isempty(readdir(root)) && error("Nonempty preflight root")
mkpath(root)
base=TOML.parsefile("configs/large40_channel_beta.toml")
function experiment(name, changes; pickup=false)
    cfg=merge(base,Dict(changes)); cfg["name"]=name
    path=joinpath(root,name*".toml")
    open(path,"w") do io; TOML.print(io,cfg); end
    out=joinpath(root,name)
    args=[path,out,"GPU"]; pickup===false || push!(args,pickup)
    run_large(args)
    check_large_output(out)
    return out
end
if stage=="thermal"
small=Dict("nx"=>32,"ny"=>32,"nz"=>32,"lx"=>600000.0,"ly"=>600000.0,
           "output_seconds"=>600.0,"diagnostic_seconds"=>60.0,"checkpoint_seconds"=>1800.0)
@testset "Optimized initializer preserves fields and constraints" begin
    cfg=merge(base,small); rho=exp.(-((1:32).-0.5)./32)
    f=2cfg["rotation_rate"]*sind(cfg["latitude"])
    old=initial_channel_turbulence(cfg,rho,f); new=initial_large_channel_turbulence(cfg,rho,f)
    for k in 1:3; @test isapprox(old[k],new[k];rtol=1e-12,atol=1e-12); end
    @test new[4]["initial_mass_divergence_relative"] < 1e-10
    @test new[4]["initial_eddy_rossby"] ≈ 0.1
end
full=experiment("restart_reference",merge(small,Dict("stop_seconds"=>3600.0)))
part=experiment("restart_first",merge(small,Dict("stop_seconds"=>1800.0)))
resumed=experiment("restart_second",merge(small,Dict("stop_seconds"=>3600.0));pickup=latest_checkpoint(part))
@testset "Restart equivalence" begin
    a=load(joinpath(full,"frame_00006.jld2")); b=load(joinpath(resumed,"frame_00006.jld2"))
    for key in ("u","v","w","theta_profile","vorticity_slices","ubar","vbar","Rxx","Rxy","Ryy")
        @test isapprox(a[key],b[key];rtol=1e-10,atol=1e-10)
    end
end
GC.gc(); CUDA.reclaim()
# Affordable 40-day test of the thermal reservoir. This small box is not a morphology result.
pilot=experiment("cap40_small_channel",Dict("nx"=>64,"ny"=>64,"nz"=>64,"lx"=>600000.0,"ly"=>600000.0,
    "stop_seconds"=>40*86400.0,"output_seconds"=>10800.0,"diagnostic_seconds"=>600.0))
rows=readlines(joinpath(pilot,"diagnostics.csv")); names=split(rows[1],','); idx=findfirst(==("cap_core_min_N2_s2"),names)
mincap=minimum(parse(Float64,split(row,',')[idx]) for row in rows[2:end])
@test mincap >= 0.5base["cap_n2"]
open(joinpath(root,"THERMAL_PASS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","minimum_40day_cap_N2"=>mincap))
end
else
thermal=TOML.parsefile(joinpath(root,"THERMAL_PASS.toml"))
thermal["status"]=="PASS" || error("Thermal checks incomplete")
mincap=thermal["minimum_40day_cap_N2"]
# Test actual production geometry, initialization cost, GPU memory, I/O and throughput.
probe=experiment("large_grid_probe",Dict("stop_seconds"=>7200.0,"output_seconds"=>3600.0,"diagnostic_seconds"=>3600.0,"checkpoint_seconds"=>7200.0))
m=TOML.parsefile(joinpath(probe,"metadata.toml"))
projected10days=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*10*86400
@test projected10days < 80*3600 # leave 40 h margin within five-day allocation
record=Dict("status"=>"PASS","minimum_40day_cap_N2"=>mincap,"estimated_10day_wall_seconds"=>projected10days,
            "note"=>"Startup throughput estimate; later turbulence may be slower. Runtime cap guard remains active.")
open(joinpath(root,"PASS.toml"),"w") do io; TOML.print(io,record); end
println("ALL LARGE-CAMPAIGN GATES PASSED")

end
