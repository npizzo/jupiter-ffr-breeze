using TOML, SHA
const SEED_CASES=("zero512_b0","zero512_b162","modes512_b162","modes1024_b162")
seed_config_path(c)=joinpath("configs","ferrel2026_seedtest_v1",c*".toml")
seed_sha(path)=bytes2hex(sha256(read(path)))
function require_allocation()
    haskey(ENV,"SLURM_JOB_ID") || error("Scientific validation requires a Unity allocation")
end
function seed_source_hashes()
    result=Dict(basename(p)=>seed_sha(p) for p in readdir(@__DIR__;join=true)
         if any(endswith(p,s) for s in (".jl",".py",".json",".sh")))
    for name in ("unit","512","1024")
        file="ferrel2026-seedtest-v1-$(name).sbatch"
        result["slurm/"*file]=seed_sha(joinpath(@__DIR__,"..","..","slurm",file))
    end
    result
end
seed_configuration_hashes()=Dict(c=>seed_sha(seed_config_path(c)) for c in SEED_CASES)
function fresh_directory(path)
    ispath(path) && error("Refusing an existing output path: $path")
    mkpath(path)
end
function verify_seed_configs()
    expected_extra=Set(["thermal_seed_model","thermal_seed_reference_nx","thermal_seed_reference_ny","precision_diagnostic_heights"])
    for c in SEED_CASES
        cfg=TOML.parsefile(seed_config_path(c))
        basecase=endswith(c,"_b0") ? "b0" : "b162"
        baseline=TOML.parsefile(joinpath("configs","ferrel2026",basecase*".toml"))
        Set(setdiff(keys(cfg),keys(baseline)))==expected_extra || error("Unexpected config keys")
        for key in keys(baseline)
            key in ("name","nx","ny","stop_seconds") && continue
            cfg[key]==baseline[key] || error("Unreviewed physical config change: $c/$key")
        end
        n=startswith(c,"modes1024") ? 1024 : 512
        cfg["name"]=="ferrel2026_seedtest_v1_"*c || error("Wrong output name")
        cfg["nx"]==cfg["ny"]==n && cfg["nz"]==96 || error("Unexpected grid")
        cfg["stop_seconds"]==21600.0 || error("Only6h tests authorized here")
        cfg["thermal_seed_model"]==(startswith(c,"zero") ? "zero" : "horizontal_modes") || error("Wrong seed type")
        cfg["thermal_seed_reference_nx"]==cfg["thermal_seed_reference_ny"]==512 || error("Wrong variance reference")
        cfg["precision_diagnostic_heights"]==[5625.0,9375.0,20000.0,60000.0,100000.0] || error("Unexpected precision planes")
    end
end
function audit_seed_source()
    verify_seed_configs()
    read(`python3 $(joinpath(@__DIR__,"verify_inheritance.py"))`,String)
end
function write_seed_gate(path;details...)
    isfile(path) && error("Refusing gate overwrite")
    d=Dict{String,Any}("status"=>"PASS", "source_hashes"=>seed_source_hashes(),
        "configuration_hashes"=>seed_configuration_hashes(),"slurm_job_id"=>ENV["SLURM_JOB_ID"])
    for (k,v) in details;d[string(k)]=v;end
    open(path,"w") do io;TOML.print(io,d);end
end
function verify_seed_gate(path)
    d=TOML.parsefile(path)
    d["status"]=="PASS" || error("Initialization gate failed")
    d["source_hashes"]==seed_source_hashes() || error("Seed-test source changed after gate")
    d["configuration_hashes"]==seed_configuration_hashes() || error("Seed-test config changed after gate")
    d
end
