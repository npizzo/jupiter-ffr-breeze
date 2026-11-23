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
    for (name,opts) in [("original",Dict()),("centered",Dict("momentum_scheme"=>"centered2","scalar_scheme"=>"centered2")),("cap",Dict("jet_design"=>"cap","jet_speed"=>.015625,"cap_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"probe_amplitude"=>1e-5,"probe_low"=>.8,"probe_high"=>1.))]
        c=merge(base,opts,Dict("nx"=>8,"ny"=>32,"nz"=>80,"lx"=>.25,"stop_time"=>.04,"initial_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"save_3d"=>true,"record_stages"=>true))
        proofs[name]=exercise(c,joinpath(out,"preflight",name),arch)
    end
    proofs["status"]="PASS";write_toml(joinpath(out,"DRIVER_PASS.toml"),proofs)
    mode=="preflight" && return
    if mode=="forensic"
        cp="data/beta-shear-forensic-v1-2138558/replay/checkpoints/checkpoint_iteration22093.jld2"
        c=merge(base,Dict("record_stages"=>true,"allow_guard_stop"=>true,"mach_stop"=>3.,"temperature_stop"=>0.,"theta_log_departure_stop"=>Inf,"profile_drift_stop"=>Inf,"diagnostic_interval"=>.02,"movie_interval"=>.1,"checkpoint_interval"=>.5,"replay_duration"=>1.,"wall_seconds"=>240.,"max_dt"=>.02))
        for (name,opts) in [("original",Dict("cap_factor"=>0.)),("restored_cap",Dict("cap_factor"=>1.))]
            run_case(merge(c,opts),joinpath(out,name);arch,pickup=cp)
        end
        return
    end
    # Native CPU operator comparison on identical fields; no new physical source.
    auditrows=Dict{String,Any}()
    bench=merge(base,Dict("nx"=>24,"ny"=>48,"nz"=>40,"lx"=>4.,"initial_rms"=>.1,"jet_speed"=>0.,"entropy_diffusion_factor"=>0.,"stop_time"=>10.,"allow_guard_stop"=>true,"save_3d"=>true))
    for (name,opts) in [("weno",Dict()),("centered_momentum",Dict("momentum_scheme"=>"centered2")),("centered_both",Dict("momentum_scheme"=>"centered2","scalar_scheme"=>"centered2"))]
        c=merge(bench,opts);mm,pp=build_model(c;arch=CPU());auditrows[name*"_initial"]=audit(mm,pp)
        println("BEGIN_BENCH ",name);flush(stdout);sim,meta=run_case(c,joinpath(out,name);arch)
        a=arrays(sim.model);set!(mm;ρ=a.rho,θ=(a.rho.*a.T).^.6 ./a.rho,u=a.uf,v=a.vf,w=a.wf,qᵗ=0.)
        Oceananigans.TimeSteppers.update_state!(mm);auditrows[name*"_final"]=audit(mm,pp)
        write_toml(joinpath(out,"NATIVE_ENERGY.toml"),auditrows);GC.gc();CUDA.reclaim()
    end
    # Static screening of both signs and designs; full-tilt cap candidate gets refinement.
    screens=Dict{String,Any}();capbase=merge(base,Dict("nx"=>8,"ny"=>128,"nz"=>80,"lx"=>.25,"jet_design"=>"cap","cap_factor"=>1.,"conduction_factor"=>1.,"top_factor"=>1.,"bottom_factor"=>0.,"stop_time"=>30.,"allow_guard_stop"=>true,"probe_low"=>.8,"probe_high"=>1.))
    cases=Tuple{String,Dict{String,Any}}[]
    for speed in [.25,.125,.0625,.03125,.015625], sign in [-1.,1.]
        push!(cases,("cap_$(sign*speed)",merge(capbase,Dict("jet_speed"=>sign*speed))))
    end
    for sign in [-1.,1.]
        push!(cases,("cap_fine_$(sign*.015625)",merge(capbase,Dict("jet_speed"=>sign*.015625,"ny"=>256,"nz"=>160,"smagorinsky_coefficient"=>.16*2^(2/3)))))
        push!(cases,("traditional_deep_$(sign*.25)",merge(capbase,Dict("jet_design"=>"deep","horizontal_rotation_factor"=>0.,"jet_speed"=>sign*.25))))
        push!(cases,("traditional_cap_$(sign*.015625)",merge(capbase,Dict("horizontal_rotation_factor"=>0.,"jet_speed"=>sign*.015625))))
    end
    for (name,c) in cases
        mm,pp=build_model(c;arch=CPU());screens[name]=spatial_screen(mm,pp)
        write_toml(joinpath(out,"JET_SCREENS.toml"),screens);println("SCREEN ",name," ",screens[name]);flush(stdout);GC.gc()
    end
    write_toml(joinpath(out,"JET_SCREENS.toml"),screens)
    # No automatic production nor arbitrary thermal/energy fixes. Save the decisions for review.
    write_toml(joinpath(out,"OUTCOME.toml"),Dict("status"=>"BENCHMARKS_AND_STATIC_SCREENS_COMPLETE","production_release"=>false,"next"=>"Review measured transport energy errors and static cap stability before paired dynamic screening"))
end
main()
