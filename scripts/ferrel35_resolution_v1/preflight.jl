using Test,TOML,JLD2,SHA,CUDA
include("run.jl")
root=ARGS[1];mkpath(root)
base=TOML.parsefile("configs/ferrel35-resolution-v1/coarse.toml")
fine=TOML.parsefile("configs/ferrel35-resolution-v1/fine.toml")
original=TOML.parsefile("configs/ferrel35-face-v3/diagnostic.toml")
for (k,v) in original
    k in ["name","max_dt"] && continue
    @test base[k]==v
    k in ["nx","ny"] || @test fine[k]==v
end
@test base["max_dt"]==fine["max_dt"]==15.0
@test fine["nx"]==2base["nx"] && fine["ny"]==2base["ny"] && fine["nz"]==base["nz"]
@test read(joinpath(@__DIR__,"face_budget.jl"))==read("scripts/ferrel35_face_v3/face_budget.jl")
@testset "Conservative staggered refinement" begin
    rng=MersenneTwister(719)
    u=randn(rng,16,16,8);v=randn(rng,16,17,8);w=randn(rng,16,16,9)
    v[:,[1,end],:].=0;w[:,:,[1,end]].=0
    M=Dict("ρu"=>u,"ρv"=>v,"ρw"=>w)
    N=Dict(k=>refine_field(a,k) for (k,a) in M)
    dc=mass_divergence(M,1,1,1);df=mass_divergence(N,.5,.5,1)
    expected=prolong_center(prolong_center(dc,1;periodic=true),2)
    @test maximum(abs,df.-expected)<1e-12
    for (k,a) in M;@test isapprox(restrict_field(N[k],k),a;atol=1e-14,rtol=1e-14);end
    q=rand(rng,16,16,8);qq=refine_field(q,"dye")
    @test minimum(qq)>=minimum(q) && maximum(qq)<=maximum(q)
    @test isapprox(restrict_field(qq,"dye"),q;atol=1e-14,rtol=1e-14)
    @test all(iszero,refine_field(zeros(16,16,8),"ρu"))
end
checkpoint=only(filter(p->endswith(p,".jld2"),readdir("data/cascade35-300-v1/cascade35_heated/segment_30/checkpoints";join=true)))
function trial(label,cfg;pickup=nothing)
    out=joinpath(root,label);conf=joinpath(root,label*".toml")
    open(conf,"w") do io;TOML.print(io,cfg);end
    args=[conf,out,"GPU"];pickup===nothing || push!(args,pickup)
    run_eddy(args);GC.gc();CUDA.reclaim()
    return out
end
bundle=joinpath(root,"day300_bundle.jld2")
trial("export",merge(base,Dict("export_bundle"=>bundle));pickup=checkpoint)
probe=merge(base,Dict("stop_seconds"=>300*86400.0+3600,"output_seconds"=>600.0,"diagnostic_seconds"=>600.0,"checkpoint_seconds"=>3600.0))
direct=trial("coarse_direct",probe;pickup=checkpoint)
identity=trial("coarse_bundle",merge(probe,Dict("restart_bundle"=>bundle)))
lastframe(out)=load(last(sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))))
@testset "Same-grid bundle reproduces direct checkpoint evolution" begin
    a,b=lastframe(direct),lastframe(identity)
    for key in ["u","v","w","theta","theta_profile","vorticity_slices","tracer_slices","ubar","vbar","Rxx","Rxy","Ryy","zonal_u","zonal_v","zonal_w","eddy_uv","eddy_uw","native_rhov","native_rhow"]
        @test isapprox(a[key],b[key];rtol=1e-9,atol=1e-10)
    end
    @test isapprox(a["native_integrated_parts"],b["native_integrated_parts"];rtol=1e-8,atol=1e-9)
    @test isapprox(a["face_integrated_momentum"],b["face_integrated_momentum"];rtol=1e-8,atol=1e-9)
    @test isapprox(a["face_integrated_heat"],b["face_integrated_heat"];rtol=1e-8,atol=1e-9)
end
pilot=trial("fine_probe",merge(fine,Dict("restart_bundle"=>bundle,"stop_seconds"=>300*86400.0+3600,"output_seconds"=>600.0,"diagnostic_seconds"=>600.0,"checkpoint_seconds"=>3600.0)))
m=TOML.parsefile(joinpath(pilot,"metadata.toml"))
estimate_hours=m["post_warmup_wall_seconds"]/m["post_warmup_simulated_seconds"]*24
@test estimate_hours<5
startframe(out)=load(joinpath(out,"frame_43200.jld2"))
a,b=startframe(identity),startframe(pilot)
ke_change=b["kinetic_energy"]/a["kinetic_energy"]-1
@test abs(ke_change)<0.1
@test isapprox(a["tracer_column"],b["tracer_column"];rtol=1e-12)
@test isapprox(a["theta_profile"],b["theta_profile"];rtol=1e-12)
source_hashes=Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in readdir(@__DIR__;join=true) if isfile(p))
open(joinpath(root,"PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","compute_location"=>"Unity","gpu"=>CUDA.name(CUDA.device()),"fine_grid"=>[512,512,64],"coarse_grid"=>[256,256,64],"max_dt"=>15.0,"source_hashes"=>source_hashes,"estimated_fine_day_gpu_hours"=>estimate_hours,"initial_fine_vs_coarse_kinetic_energy_fraction"=>ke_change,"checkpoint_sha256"=>bytes2hex(open(sha256,checkpoint)),"bundle_sha256"=>bytes2hex(open(sha256,bundle))))
end
println("RESOLUTION_PREFLIGHT_PASS estimated fine day GPU hours=",estimate_hours," remap KE fraction=",ke_change);flush(stdout)
