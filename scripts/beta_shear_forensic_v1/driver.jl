include("model.jl")
include("initialization.jl")
include("capture.jl")
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
 mode,out=ARGS[1:2];@assert !ispath(out);mkpath(out)
 arch=mode=="preflight" ? CPU() : GPU()
 mode=="preflight" || CUDA.functional() || error("GPU required")
 c=TOML.parsefile(joinpath(@__DIR__,"pilot.toml"))
 tiny=merge(c,Dict("nx"=>8,"ny"=>32,"stop_time"=>.04,"restart_time"=>.02,"initial_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"field_interval"=>.02))
 err=exercise(tiny,joinpath(out,"preflight"),arch)
 # Exercise forensic serialization separately without taking a nonphysical step.
 m,p=build_model(tiny;arch);bad=joinpath(out,"capture_check");mkpath(bad);save_bad(m,bad,"serialization path check on a valid state")
 a=load(joinpath(bad,"FIRST_BAD_STATE.jld2"));@assert minimum(a["rho"])>0 && all(isfinite,a["rhotheta"])
 write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err))
 mode=="preflight" && return
 c=merge(c,Dict("stop_time"=>160.,"restart_time"=>150.,"wall_seconds"=>600.,"diagnostic_interval"=>.05,"movie_interval"=>.5,"checkpoint_interval"=>1.,"save_3d"=>false,"mach_stop"=>3.))
 cp="data/beta-shear-v1-pilot-2134379/parent/checkpoints/checkpoint_iteration20549.jld2"
 try
  run_case(c,joinpath(out,"replay");arch,pickup=cp)
 catch err
  write_toml(joinpath(out,"REPLAY_OUTCOME.toml"),Dict("status"=>"STOPPED","reason"=>sprint(showerror,err)))
  println("REPLAY_STOPPED ",sprint(showerror,err));flush(stdout)
 end
end
main()
