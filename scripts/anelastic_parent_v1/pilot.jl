include("model.jl")
CUDA.functional()||error("GPU required")
root,cpureceipt=ARGS;mkpath(root)
receipt=TOML.parsefile(cpureceipt)
@assert receipt["status"]=="PASS"
@assert receipt["model_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl"))))
started=time()
c=TOML.parsefile("configs/anelastic-parent-v1/parent.toml")
small=merge(c,Dict("nx"=>8,"ny"=>8,"nz"=>32,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>120.))
s,m=run_case(small,joinpath(root,"rest");rest=true);@assert maximum(abs,arrays(s.model).w)<1e-9
s1,m1=run_case(small,joinpath(root,"restart_first"))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"restart_first/checkpoints");join=true)))
s2,m2=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"restart_resumed");pickup=cp)
s3,m3=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"restart_continuous"))
a=arrays(s2.model);b=arrays(s3.model);err=Dict{String,Float64}()
for k in (:u,:v,:w,:theta)
    e=maximum(abs,getproperty(a,k).-getproperty(b,k));@assert e<1e-8;err[string(k)]=e
end
s=nothing;s1=nothing;s2=nothing;s3=nothing;a=nothing;b=nothing;GC.gc();CUDA.reclaim()
# Actual target grid, evolved as part of the parent rather than discarded timing work.
initial=merge(c,Dict("stop_time"=>5.,"checkpoint_interval"=>5.,"wall_seconds"=>300.))
s,m=run_case(initial,joinpath(root,"parent_initial"));@assert m["status"]=="TARGET_REACHED_REVIEW_REQUIRED"
rows=[parse.(Float64,split(l,',')) for l in readlines(joinpath(root,"parent_initial/diagnostics.csv"))[2:end]]
rate=(rows[end][4]-rows[end-2][4])/(rows[end][1]-rows[end-2][1])
write_toml(joinpath(root,"COST.toml"),Dict("seconds_per_unit_initial"=>rate,"gpu"=>CUDA.name(CUDA.device()),"elapsed_validation_seconds"=>time()-started,"warning"=>"Early flow only; measure turbulent throughput and turnover count before any branch"))
remaining=min(c["wall_seconds"],3300-(time()-started))
remaining>120||error("Validation used available budget; no parent continuation")
remaining/(1.5rate)>=100||error("Projected progress under t100 within remaining budget: stop review")
write_toml(joinpath(root,"VALIDATION_PASS.toml"),Dict("status"=>"PASS","model_sha256"=>receipt["model_sha256"],"restart_errors"=>err,"gpu"=>CUDA.name(CUDA.device()),"source_scope"=>"Implementation only. No production release authorization."))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"parent_initial/checkpoints");join=true)))
s=nothing;GC.gc();CUDA.reclaim()
s,m=run_case(merge(c,Dict("wall_seconds"=>remaining)),joinpath(root,"parent");pickup=cp)
write_toml(joinpath(root,"REVIEW_REQUIRED.toml"),Dict("status"=>"REVIEW_REQUIRED","final_time"=>m["final_time"],"turnovers_depth"=>m["turnovers_depth"],"elapsed_total_seconds"=>time()-started,"next"=>"Review convection, thermal drift, resolution, Rossby number and cost. Do not launch seeded cases automatically."))
println("ANELASTIC_PARENT_REVIEW_REQUIRED");flush(stdout)
