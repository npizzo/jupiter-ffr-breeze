function main()
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
c=merge(c,Dict("stop_time"=>6000.,"initial_dt"=>.001,"max_dt"=>.2,"wall_seconds"=>2700.,"diagnostic_interval"=>1.,"movie_interval"=>5.,"checkpoint_interval"=>50.))
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
# This same continuation branch is exercised at small size in CPU/GPU validation.
function continuation(c, root, checkpoint; arch, target, wall_seconds, expected_start, source_frame=nothing)
    @assert isfile(checkpoint)
    if source_frame !== nothing
        # Read the actual production checkpoint at its saved time without advancing.
        mkpath(root)
        restored_model, restored_params=build(c;arch)
        loaded=Simulation(restored_model;Δt=c["initial_dt"],stop_time=target)
        saved_pickup=SavedStep(c["initial_dt"])
        add_callback!(loaded,saved_pickup,IterationInterval(1);name=:saved_step)
        set!(loaded;checkpoint)
        Oceananigans.initialize!(loaded) # auxiliary fields only; no time_step! or run!
        @assert isapprox(loaded.model.clock.time,expected_start;atol=1e-9)
        a=arrays(loaded.model);rho=a.rho
        ub=dropdims(sum(rho .* a.u;dims=3)./sum(rho;dims=3);dims=3)
        vb=dropdims(sum(rho .* a.v;dims=3)./sum(rho;dims=3);dims=3)
        jldopen(source_frame,"r") do file
            @assert Float32.(ub)==file["ubar"]
            @assert Float32.(vb)==file["vbar"]
            @assert isapprox(vec(mean(a.theta;dims=(1,2))),file["theta_profile"];rtol=1e-12,atol=1e-14)
            z=[(k-.5)*c["height"]/c["nz"] for k=1:c["nz"]]
            planes=[argmin(abs.(z.-zz)) for zz in (.2,.6,.9)]
            @assert Float32.((Array(interior(loaded.model.tracers.dye))./rho)[:,:,planes])==file["tracer_slices"]
        end
        write_toml(joinpath(root,"PICKUP_PASS.toml"),Dict("status"=>"PASS","checkpoint"=>checkpoint,"source_frame"=>source_frame,"start_time"=>expected_start,"scope"=>"Production-grid checkpoint restored; velocity depth means, theta profile and tracer planes equal saved t2000 fields before advancement"))
        loaded=nothing;restored_model=nothing;restored_params=nothing;a=nothing;rho=nothing;ub=nothing;vb=nothing;GC.gc()
        arch isa GPU && CUDA.reclaim()
    end
    run_case(merge(c,Dict("stop_time"=>target,"wall_seconds"=>wall_seconds)),joinpath(root,"parent");arch,pickup=checkpoint)
end
# Use the same continuation wrapper, diagnostic outputs and checkpoint code in preflight.
small_frame=joinpath(root,"restart_first","frame_00000020.jld2")
continued, cm=continuation(small,joinpath(root,"continuation_test"),cp;arch,target=.4,wall_seconds=120.,expected_start=.2,source_frame=small_frame)
ca=arrays(continued.model)
@assert maximum(abs,ca.u.-arrays(s3.model).u)<1e-8
@assert isfile(joinpath(root,"continuation_test/parent/diagnostics.csv"))
if preflight
    write_toml(joinpath(root,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_errors"=>err,"scope"=>"Exact continuation driver, imports, diagnostics, three-depth output/readback and continuous-versus-restart validation on CPU"))
    println("EXACT_DRIVER_CPU_PASS");flush(stdout);return nothing
end
s=nothing;s1=nothing;s2=nothing;s3=nothing;a=nothing;b=nothing;continued=nothing;ca=nothing;GC.gc();CUDA.reclaim()
spec=TOML.parsefile("configs/anelastic-sequence-v3/continuation.toml")
@assert bytes2hex(sha256(read(spec["checkpoint"])))==spec["checkpoint_sha256"]
@assert spec["start_time"]==2000. && spec["stop_time"]==6000.
remaining=min(3000.0, 3300.0 - (time()-started));@assert remaining>120
write_toml(joinpath(root,"VALIDATION_PASS.toml"),Dict("status"=>"PASS","model_sha256"=>receipt["model_sha256"],"restart_errors"=>err,"gpu"=>CUDA.name(CUDA.device()),"physics"=>"Unchanged weak-cap baseline; no other branches"))
s,m=continuation(c,root,spec["checkpoint"];arch,target=spec["stop_time"],wall_seconds=remaining,expected_start=spec["start_time"],source_frame=spec["source_frame"])
@assert m["final_time"]>spec["start_time"] "Continuation made no progress"
write_toml(joinpath(root,"REVIEW_REQUIRED.toml"),Dict("status"=>m["status"],"final_time"=>m["final_time"],"additional_turnovers_depth"=>m["turnovers_depth"],"elapsed_total_seconds"=>time()-started,"next"=>"Review thermal equilibration, cap and large-scale energy. No automatic continuation or parameter change."))
println("ANELASTIC_CONTINUATION_REVIEW_REQUIRED");flush(stdout)

end

include("calibration.jl")
main()
