using JLD2,TOML,Statistics
root,out=ARGS;mkpath(out)
labels=["lower","interface","upper","full"]
open(joinpath(out,"layer_statistics.csv"),"w") do io
 println(io,"case,time_days,layer,signed_skew_prime,signed_skew_prime_interior,rms_prime_over_abs_f,mean_zeta_over_f")
 for layer in ("stable","neutral"),shear in ("cyclonic","none","anticyclonic")
  name="$(layer)_$(shear)";dir=joinpath(root,name)
  m=TOML.parsefile(joinpath(dir,"metadata.toml"));m["status"]=="COMPLETE"||error("Incomplete case $name")
  target=joinpath(out,"layers",name);mkpath(target)
  for file in sort(filter(x->startswith(basename(x),"frame_")&&endswith(x,".jld2"),readdir(dir;join=true)))
   d=load(file);ls=d["layer_data"]["layers"]
   for label in labels
    r=ls[label]
    println(io,join((name,d["t"]/86400,label,r["signed_skew_prime"],r["signed_skew_prime_interior"],r["rms_prime_over_abs_f"],r["mean_zeta_over_f"]),","))
   end
   u=cat([ls[l]["u"] for l in labels]...;dims=3)
   v=cat([ls[l]["v"] for l in labels]...;dims=3)
   jldsave(joinpath(target,basename(file));u,v,t=d["t"])
  end
 end
end
println("LAYER_STATISTICS_PREPARED")
