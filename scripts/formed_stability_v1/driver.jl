function main()
    mode,out=ARGS[1:2];isdir(out)&&error("Refuse existing output");mkpath(out)
    cfg=TOML.parsefile("configs/formed-stability-v1.toml")
    function run(name,c;checkpoint=nothing,arch="CPU")
        path=joinpath(out,name*".toml");open(path,"w") do io;TOML.print(io,c);end
        args=[path,joinpath(out,name),arch];checkpoint===nothing||push!(args,checkpoint)
        run_eddy(args)
    end
    if mode=="preflight"
        c=merge(cfg,Dict("nx"=>32,"ny"=>32,"nz"=>32,"turbulence_wavelength"=>4000000.,"stop_seconds"=>60.,"output_seconds"=>30.,"diagnostic_seconds"=>30.,"checkpoint_seconds"=>60.,"full_field_seconds"=>60.,"controller_seconds"=>30.,"max_dt"=>5.,"initial_dt"=>5.,"wall_limit"=>180.,"deformation_enabled"=>false))
        first=run("initial",c);cp=only(readdir(joinpath(out,"initial/checkpoints");join=true));first=nothing;GC.gc()
        for (name,eps) in [("control",0.),("half",.005),("full",.01)]
            q=merge(c,Dict("stop_seconds"=>120.,"deformation_enabled"=>true,"deformation_fraction"=>eps))
            ss=run(name,q;checkpoint=cp)
            @assert ss.model.clock.time>=120 .- 1e-8
            frame=load(joinpath(out,name,"frame_00004.jld2"));@assert size(frame["vorticity_slices"],3)==3 && size(frame["tracer_slices"],3)==3
            rows=readlines(joinpath(out,name,"diagnostics.csv"));@assert length(rows)>1 && all(all(isfinite,parse.(Float64,split(s,','))) for s in rows[2:end])
            ss=nothing;GC.gc()
        end
        cp2=only(readdir(joinpath(out,"control/checkpoints");join=true))
        q=merge(c,Dict("stop_seconds"=>180.))
        resumed=run("resumed",q;checkpoint=cp2);continuous=run("continuous",q;checkpoint=cp)
        err=maximum(maximum(abs,Array(interior(f)) .- Array(interior(Oceananigans.prognostic_fields(continuous.model)[k]))) for (k,f) in pairs(Oceananigans.prognostic_fields(resumed.model)))
        @assert err<1e-8
        open(joinpath(out,"DRIVER_PASS.toml"),"w") do io;TOML.print(io,Dict("status"=>"PASS","restart_error"=>err,"scope"=>"Actual driver and imports, both perturbation amplitudes, output/readback, serialized restart"));end
    else
        mode in ["control","half","full"]||error("Unknown mode");CUDA.functional()||error("GPU required")
        cp="data/fplane-point4-v2/late_heated/checkpoints/checkpoint_iteration289093.jld2"
        @assert bytes2hex(open(sha256,cp))==strip(read("configs/formed-stability-v1-checkpoint.sha256",String))
        cfg["source_state"]="data/fplane-point4-v2/late_heated/fields3d/state_000008640000.jld2"
        cfg["deformation_fraction"]=mode=="half" ? .005 : mode=="full" ? .01 : 0.
        ss=run("parent",cfg;checkpoint=cp,arch="GPU")
        @assert ss.model.clock.time>=cfg["stop_seconds"]-1e-8 "Wall limit; partial result"
        open(joinpath(out,"COMPLETE.toml"),"w") do io;TOML.print(io,Dict("status"=>"COMPLETE","final_day"=>ss.model.clock.time/86400));end
    end
end
include("run.jl")
main()
