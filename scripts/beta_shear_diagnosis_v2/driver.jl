# Diagnostic ablations only. No production submission or changed original files.
include("model.jl")
include("initialization.jl")
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
 base=merge(TOML.parsefile(joinpath(@__DIR__,"pilot.toml")),Dict("nx"=>8,"lx"=>.25,"initial_rms"=>0.,"perturbation"=>0.,"stop_time"=>30.,"initial_dt"=>.01,"wall_seconds"=>180.,"movie_interval"=>5.,"field_interval"=>10.,"checkpoint_interval"=>10.,"mach_stop"=>.5))
 variants=[("baseline",Dict()),("half_dt",Dict("max_dt"=>1/120)),("no_entropy_diffusion",Dict("entropy_diffusion_factor"=>0.)),("no_horizontal_rotation",Dict("horizontal_rotation_factor"=>0.))]
 arch=mode=="preflight" ? CPU() : GPU()
 mode=="preflight" || CUDA.functional() || error("GPU required")
 checks=Dict{String,Any}()
 for (name,overrides) in variants
  c=merge(base,overrides,Dict("nx"=>8,"ny"=>32,"lx"=>.25,"stop_time"=>.04,"restart_time"=>.02,"initial_dt"=>.001,"diagnostic_interval"=>.01,"movie_interval"=>.01,"checkpoint_interval"=>.02,"field_interval"=>.02))
  checks[name]=exercise(c,joinpath(out,"preflight",name),arch)
 end
 checks["status"]="PASS"
 checks["files"]=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if endswith(f,".jl")||endswith(f,".toml"))
 write_toml(joinpath(out,"DRIVER_PASS.toml"),checks)
 mode=="preflight" && return
 for (name,overrides) in variants
  println("BEGIN_CASE ",name);flush(stdout)
  try
   run_case(merge(base,overrides),joinpath(out,name);arch)
  catch e
   write_toml(joinpath(out,name*"_ERROR.toml"),Dict("error"=>sprint(showerror,e)))
   println("CASE_ERROR ",name," ",sprint(showerror,e));flush(stdout)
  end
  GC.gc()
 end
end
main()
