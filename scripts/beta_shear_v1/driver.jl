function config(case)
    @assert case in ["control","shear","pilot"]
    c=TOML.parsefile("configs/beta-shear-v1-"*case*".toml");@assert c["omega"]==3.;c
end
small(c)=merge(c,Dict("nx"=>32,"ny"=>64,"nz"=>80,"stop_time"=>.04,"restart_time"=>.02,"initial_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"field_interval"=>.02,"wall_seconds"=>480.))
function difference(a,b)
    fa=Oceananigans.prognostic_fields(a);fb=Oceananigans.prognostic_fields(b)
    @assert keys(fa)==keys(fb)
    maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa))
end
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in rows[2:end])
    fs=sort(filter(f->startswith(f,"frame_"),readdir(out)));@assert !isempty(fs)
    a=load(joinpath(out,last(fs)));@assert size(a["vorticity_slices"],3)==3 && all(isfinite,a["tracer_slices"])
    @assert a["ly"]==8. && haskey(a,"zonal_u")
    fs=sort(readdir(joinpath(out,"fields3d");join=true));@assert !isempty(fs)
    a=load(last(fs));@assert minimum(a["density"])>0 && all(isfinite,a["temperature"])
    @assert maximum(abs,a["v_face"][:,[1,end],:])==0. && maximum(abs,a["w_face"][:,:,[1,end]])==0.
end
function exercise(c,out,arch)
    a,_=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,_=run_case(c,joinpath(out,"restart");arch,pickup=cp)
    e,_=run_case(c,joinpath(out,"direct");arch)
    err=difference(b.model,e.model);@assert err<1e-7
    for n in ["first","restart","direct"];readback(joinpath(out,n));end
    err
end
function main()
    mode,out=ARGS[1:2];ispath(out)&&error("Refuse existing output");mkpath(out)
    if mode=="preflight"
        errors=Dict{String,Any}()
        for case in ["control","shear","pilot"]
            c=config(case);operator_checks(c,out)
            m,p=build_model(merge(c,Dict("nx"=>64,"ny"=>128,"nz"=>80));arch=CPU(),perturb=false)
            errors[case*"_balance"]=balance_check(m,p);m=nothing;GC.gc()
        end
        for case in ["control","shear","pilot"]
            errors[case*"_restart_error"]=exercise(small(config(case)),joinpath(out,case),CPU())
        end
        # With jets disabled and original geometry, the physical implementation must agree.
        c=merge(small(config("control")),Dict("ly"=>4.,"ny"=>32));a,_=build_model(c;arch=CPU());b,_=Previous.build_model(c;arch=CPU());err=difference(a,b);@assert err<1e-9
        errors["original_no_shear_initial_state_error"]=err
        errors["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),errors)
    elseif mode=="run"
        case=ARGS[3];c=config(case);CUDA.functional()||error("GPU required")
        err=exercise(small(c),joinpath(out,"gpu_check"),GPU())
        write_toml(joinpath(out,"GPU_PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"restart_error"=>err,"omega"=>3.,"beta"=>c["beta"]))
        GC.gc();CUDA.reclaim()
        _,meta=run_case(c,joinpath(out,"parent"))
        write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>meta["status"],"final_time"=>meta["final_time"],"case"=>case))
    else;error("Unknown mode");end
end
include("model.jl")
include("initialization.jl")
include("checks.jl")
module Previous
include("../rotation_channel_v1/model.jl")
include("../rotation_channel_v1/initialization.jl")
end
main()
