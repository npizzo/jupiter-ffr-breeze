using Test,TOML,SHA
include("run.jl")
root=only(ARGS);out=joinpath(root,"preflight");mkpath(out)
TOML.parsefile("data/fplane-sign-cpu-v1/PASS.toml")["status"]=="PASS"||error("CPU gate")
for name in ["baseline","plus","mirror"]
    cfg=TOML.parsefile("configs/fplane-sign-v1/plus.toml")
    cfg["name"]=name;cfg["stop_seconds"]=90*86400.+600.;cfg["diagnostic_seconds"]=60.;cfg["output_seconds"]=60.;cfg["full_field_seconds"]=600.
    cfg["native_budget"]=name!="baseline";cfg["reduce_exports"]=false
    cfg["reflect_x"]=name=="mirror";cfg["f_sign"]=name=="mirror" ? -1 : 1
    path=joinpath(out,name*".toml");open(path,"w") do io;TOML.print(io,cfg);end
    run_eddy([path,joinpath(out,name),"GPU"]);GC.gc();CUDA.reclaim()
end
state(d,t)=load(joinpath(out,d,"fields3d",@sprintf("state_%012d.jld2",Int(t))))
a=state("baseline",90*86400+600);b=state("plus",90*86400+600);c=state("mirror",90*86400+600)
for k in ["u_face","v_face","w_face","theta","dye_density"];@test a[k]==b[k];end
errors=Dict{String,Any}()
for k in ["u_face","v_face","w_face","theta","dye_density"]
    q=b[k];nx=size(q,1);expected=k=="u_face" ? -q[mod1.(2 .- (1:nx),nx),:,:] : reverse(q;dims=1)
    # Temperature error is normalized by perturbations, not the large background.
    scale=k=="theta" ? max(maximum(abs,q.-mean(q;dims=(1,2))),1e-5) : max(maximum(abs,q),1e-8)
    err=maximum(abs,c[k].-expected)/scale;errors[k]=err
    @test err<1e-6
end
meta=TOML.parsefile(joinpath(out,"plus","metadata.toml"))
estimate=meta["post_warmup_wall_seconds"]/meta["post_warmup_simulated_seconds"]*21600/3600
@test estimate<1.0 # Each six-hour run has a hard 90-minute allocation.
files=filter(f->endswith(f,".jl")||endswith(f,".py"),readdir(@__DIR__))
pass=Dict("status"=>"PASS","grid"=>[256,256,64],"gpu"=>CUDA.name(CUDA.device()),"estimated_six_hour_wall_hours"=>estimate,"reflection_relative_errors"=>errors,"physical_state_bitwise_unchanged"=>true,
 "source_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in files),
 "config_hashes"=>Dict("configs/fplane-sign-v1/$n.toml"=>bytes2hex(sha256(read("configs/fplane-sign-v1/$n.toml"))) for n in ["plus","minus"]),
 "bundle_sha256"=>bytes2hex(open(sha256,"data/fplane-point4-v2/preflight/day90_heated_bundle.jld2")))
open(joinpath(out,"GPU_PASS.toml"),"w") do io;TOML.print(io,pass);end
println("SIGN_GPU_PREFLIGHT_PASS");flush(stdout)
