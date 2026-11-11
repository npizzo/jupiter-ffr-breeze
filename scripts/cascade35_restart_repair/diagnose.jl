using JLD2,LinearAlgebra,TOML
root="data/cascade35-preflight-v1"
a=load(joinpath(root,"restart_full/frame_00002.jld2"));b=load(joinpath(root,"restart_second/frame_00002.jld2"))
rows=Dict{String,Any}()
for k in sort(collect(intersect(keys(a),keys(b))))
 x,y=a[k],b[k]
 x isa AbstractArray{<:Number} && y isa AbstractArray{<:Number} || continue
 size(x)==size(y) || continue
 delta=norm(x-y);scale=max(norm(x),norm(y));mx=maximum(abs,x-y)
 rows[k]=Dict("l2_difference"=>delta,"l2_scale"=>scale,"relative_l2"=>delta/max(scale,eps()),"max_difference"=>mx,"max_scale"=>max(maximum(abs,x),maximum(abs,y)),"old_test_pass"=>isapprox(x,y;rtol=1e-9,atol=1e-10))
 println(k," ",rows[k])
end
for name in ["restart_full","restart_first","restart_second"]
 p=joinpath(root,name,"metadata.toml");m=TOML.parsefile(p)
 println("METADATA ",name," ",Dict(k=>v for (k,v) in m if occursin("step",k)||occursin("time",k)||occursin("checkpoint",k)||k in ["simulated_seconds","segment_start_seconds","status"]))
end
for k in ["mass_divergence_max","mass_flux_reconstruction_max","top_streamfunction_max"]
 println("RESIDUAL ",k," uninterrupted=",a[k]," restarted=",b[k])
end
out="data/cascade35-restart-audit-v1";mkpath(out)
open(joinpath(out,"differences.toml"),"w") do io;TOML.print(io,rows);end
