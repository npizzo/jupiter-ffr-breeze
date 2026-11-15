using Test,TOML
include("run.jl")
root=only(ARGS);mkpath(root)
# Independent mathematical tests of centering, reflection and stress parity.
n=32;x=reshape(((1:n).-0.5)./n,n,1);y=reshape(((1:n).-0.5)./n,1,n)
u=sin.(2pi.*x).*cos.(pi.*y);v=cos.(2pi.*x).*sin.(pi.*y)
z=nb_ccurl(u,v,1/n,1/n)
@test isapprox(nb_ccurl(-reverse(u;dims=1),reverse(v;dims=1),1/n,1/n),-reverse(z;dims=1);atol=1e-12)
c=nb_stress(u.^2,u.*v,v.^2,1/n,1/n)
@test isapprox(nb_stress(reverse(u.^2;dims=1),-reverse(u.*v;dims=1),reverse(v.^2;dims=1),1/n,1/n),-reverse(c;dims=1);atol=1e-10)
@test maximum(abs,nb_stress(ones(n,n),zeros(n,n),ones(n,n),1/n,1/n))==0
for audit in [false,true]
 cfg=TOML.parsefile("configs/fplane-sign-v1/plus.toml");delete!(cfg,"restart_bundle")
 cfg["nx"]=64;cfg["ny"]=64;cfg["turbulence_wavelength"]=2400000.;cfg["stop_seconds"]=120.;cfg["output_seconds"]=60.;cfg["diagnostic_seconds"]=60.;cfg["full_field_seconds"]=60.;cfg["max_dt"]=10.;cfg["initial_dt"]=10.;cfg["native_budget"]=audit;cfg["reduce_exports"]=false;cfg["name"]=audit ? "audit" : "baseline"
 path=joinpath(root,cfg["name"]*".toml");open(path,"w") do io;TOML.print(io,cfg);end
 run_eddy([path,joinpath(root,cfg["name"]),"CPU"]);GC.gc()
end
name=last(sort(readdir(joinpath(root,"baseline","fields3d"))));a=load(joinpath(root,"baseline","fields3d",name));b=load(joinpath(root,"audit","fields3d",name))
for k in ["u_face","v_face","w_face","theta","dye_density"];@test a[k]==b[k];end
p=load(last(sort(readdir(joinpath(root,"audit","native_vorticity");join=true)[endswith.(readdir(joinpath(root,"audit","native_vorticity");join=true),".jld2")])) )
@test p["checks"]>0
open(joinpath(root,"PASS.toml"),"w") do io;TOML.print(io,Dict("status"=>"PASS","bitwise_observer_unchanged"=>true,"full_domain_relative_closure"=>p["max_relative_step_residual"],"checks"=>p["checks"]));end
println("SIGN_CPU_CHECK_PASS");flush(stdout)
