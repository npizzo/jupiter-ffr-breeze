include("calibration.jl")
preflight = !isempty(ARGS) && ARGS[1]=="--preflight"
root,cpureceipt = preflight ? ARGS[2:3] : ARGS
arch = preflight ? CPU() : GPU()
preflight || CUDA.functional() || error("GPU required")
mkpath(root)
receipt=TOML.parsefile(cpureceipt)
@assert receipt["status"]=="PASS"
@assert receipt["model_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl"))))
started=time()
c=TOML.parsefile("configs/anelastic-sequence-v3/baseline.toml")
c=merge(c,Dict("stop_time"=>2000.,"initial_dt"=>.001,"max_dt"=>.2,"wall_seconds"=>2700.,"diagnostic_interval"=>1.,"movie_interval"=>5.,"checkpoint_interval"=>50.))
design=TOML.parsefile(joinpath(dirname(cpureceipt),"DESIGN.toml"))["weak"]
@assert design["reference_downward_diffusion_peak_over_input"]<1 "Cap thermal forcing needs redesign before calibration"
@assert .99<design["initial_cap_cooling_over_bottom_input"]<1.01

small=merge(c,Dict("nx"=>8,"ny"=>8,"nz"=>32,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>120.))
s,m=run_case(small,joinpath(root,"rest");arch,rest=true);@assert maximum(abs,arrays(s.model).w)<1e-9
s1,m1=run_case(small,joinpath(root,"restart_first");arch)
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"restart_first/checkpoints");join=true)))
s2,m2=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"restart_resumed");arch,pickup=cp)
s3,m3=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"restart_continuous");arch)
a=arrays(s2.model);b=arrays(s3.model);err=Dict{String,Float64}()
for k in (:u,:v,:w,:theta)
    e=maximum(abs,getproperty(a,k).-getproperty(b,k));@assert e<1e-8;err[string(k)]=e
end
# Exercise and read back the real diagnostics and output formats on either architecture.
for name in ("rest","restart_first","restart_resumed","restart_continuous")
    lines=readlines(joinpath(root,name,"diagnostics.csv"));ncols=length(split(lines[1],','))
    @assert length(lines)>1
    for row in lines[2:end]
        vals=parse.(Float64,split(row,','));@assert length(vals)==ncols && all(isfinite,vals)
    end
    frames=filter(x->startswith(basename(x),"frame_")&&endswith(x,".jld2"),readdir(joinpath(root,name);join=true))
    @assert !isempty(frames)
    jldopen(last(frames),"r") do file
        @assert size(file["vorticity_slices"],3)==3 && size(file["tracer_slices"],3)==3
        @assert all(isfinite,file["w_spectrum_horizontal"])
    end
end
if preflight
    write_toml(joinpath(root,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_errors"=>err,"scope"=>"Exact driver imports, spectral diagnostics, three-depth output readback, checkpoint and restart on CPU; GPU checks still required."))
    println("EXACT_DRIVER_CPU_PASS");flush(stdout)
    exit(0)
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
remaining/(1.5rate)>=1000||error("Projected progress under t1000 within remaining budget: stop review")
write_toml(joinpath(root,"VALIDATION_PASS.toml"),Dict("status"=>"PASS","model_sha256"=>receipt["model_sha256"],"restart_errors"=>err,"gpu"=>CUDA.name(CUDA.device()),"source_scope"=>"Implementation only. No production release authorization."))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"parent_initial/checkpoints");join=true)))
s=nothing;GC.gc();CUDA.reclaim()
s,m=run_case(merge(c,Dict("wall_seconds"=>remaining)),joinpath(root,"parent");pickup=cp)
write_toml(joinpath(root,"REVIEW_REQUIRED.toml"),Dict("status"=>"REVIEW_REQUIRED","final_time"=>m["final_time"],"turnovers_depth"=>m["turnovers_depth"],"elapsed_total_seconds"=>time()-started,"next"=>"Review convection, thermal drift, resolution, Rossby number and cost. Do not launch seeded cases automatically."))
println("ANELASTIC_SEQUENCE_CALIBRATION_REVIEW_REQUIRED");flush(stdout)
