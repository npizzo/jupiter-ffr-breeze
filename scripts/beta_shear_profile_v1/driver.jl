include("model.jl")
include("initialization.jl")
include("cap_initialization.jl")
include("guards.jl")
include("capture.jl")
include("native_audit.jl")
include("stages.jl")
include("probe.jl")
include("checks.jl")
include("run.jl")
function readback(out)
    lines=readlines(joinpath(out,"diagnostics.csv"));@assert length(lines)>2
    @assert all(all(isfinite,parse.(Float64,split(l,','))) for l in lines[2:end])
    fs=filter(f->startswith(f,"frame_"),readdir(out;join=false));@assert !isempty(fs)
    a=load(joinpath(out,last(sort(fs))));@assert size(a["vorticity_slices"],3)==3 && all(isfinite,a["tracer_slices"])
    fs=filter(f->startswith(f,"yz_"),readdir(out));a=load(joinpath(out,last(sort(fs))));@assert minimum(a["rho"])>0
    fs=filter(f->startswith(f,"full_"),readdir(out));@assert !isempty(fs);a=load(joinpath(out,last(sort(fs))))
    @assert maximum(abs,a["v_face"][:,[1,end],:])==0 && maximum(abs,a["w_face"][:,:,[1,end]])==0
    @assert length(readlines(joinpath(out,"stage_terms.csv")))>2
end
function exercise(c,out,arch)
    _,_=run_case(merge(c,Dict("stop_time"=>.02)),joinpath(out,"first");arch)
    cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
    b,_=run_case(merge(c,Dict("replay_duration"=>.02)),joinpath(out,"restart");arch,pickup=cp)
    a,_=run_case(c,joinpath(out,"direct");arch)
    fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
    err=maximum(maximum(abs,Array(interior(fa[k])).-Array(interior(fb[k]))) for k in keys(fa));@assert err<1e-7
    for x in ["first","restart","direct"];readback(joinpath(out,x));end
    err
end
function main()
    mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
    arch=mode=="preflight" ? CPU() : GPU();mode=="preflight" || CUDA.functional() || error("GPU required")
    base=merge(TOML.parsefile(joinpath(@__DIR__,"pilot.toml")),Dict("initial_rms"=>0.,"perturbation"=>0.,"initial_dt"=>.001,"max_dt"=>1/120,"wall_seconds"=>240.,"diagnostic_interval"=>.25,"movie_interval"=>1.,"checkpoint_interval"=>10.,"save_3d"=>false,"mach_stop"=>.5,"temperature_stop"=>.5,"profile_drift_stop"=>.05,"heat"=>0.,"conductivity"=>.0004*.03125/(GRAV/2.5)))
    proofs=Dict{String,Any}("guards"=>guard_checks(base),"sources"=>source_checks(base))
    for (name,opts) in [("original",Dict()),("traditional",Dict("horizontal_rotation_factor"=>0.)),("cap",Dict("jet_design"=>"cap","jet_speed"=>.015625,"cap_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"probe_amplitude"=>1e-5,"probe_low"=>.8,"probe_high"=>1.))]
        c=merge(base,opts,Dict("nx"=>8,"ny"=>32,"nz"=>80,"lx"=>.25,"stop_time"=>.04,"initial_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"save_3d"=>true,"record_stages"=>true))
        proofs[name]=exercise(c,joinpath(out,"preflight",name),arch)
    end
    proofs["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),proofs)
    mode=="preflight" && return
    # Bounded jet-only mean-profile acceptance tests, not convection/LSV experiments.
    base=merge(base,Dict("nx"=>8,"ny"=>128,"nz"=>80,"lx"=>.25,"stop_time"=>30.,"wall_seconds"=>240.,"allow_guard_stop"=>true,"cap_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"bottom_factor"=>0.,"save_3d"=>false))
    outcomes=Dict{String,Any}()
    for (name,opts) in [("cap_anticyclonic",Dict("jet_design"=>"cap","jet_speed"=>.015625)),("cap_cyclonic",Dict("jet_design"=>"cap","jet_speed"=>-.015625)),("traditional_anticyclonic",Dict("horizontal_rotation_factor"=>0.,"jet_speed"=>.25)),("traditional_cyclonic",Dict("horizontal_rotation_factor"=>0.,"jet_speed"=>-.25))]
        cfg=merge(base,opts);mm,pp=build_model(cfg;arch=CPU());screen=spatial_screen(mm,pp)
        # No dynamic test of an already-invalid balanced cap state.
        if screen["cap_PV_min"]<=0 || screen["cap_N2_min"]<=0 || screen["Tmin"]<=.5 || screen["vertical_absolute_vorticity_min"]<=0
            outcomes[name]=Dict("status"=>"REJECTED_INITIAL_STATE","screen"=>screen)
        else
            println("BEGIN_PROFILE ",name);flush(stdout)
            sim,meta=run_case(cfg,joinpath(out,name);arch);outcomes[name]=meta
        end
        write_toml(joinpath(out,"PROFILE_OUTCOMES.toml"),outcomes);GC.gc();CUDA.reclaim()
    end
end
main()
