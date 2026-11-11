include("run.jl")
include("common.jl")
include("../ferrel2026_v3/common.jl")
include("check_case.jl")

require_allocation()
length(ARGS)==4 || error("usage: campaign.jl CASE OUTPUT_ROOT SEED_GATE_ROOT V3_PREFLIGHT_ROOT")
case,root,gate,prior=ARGS
case in SEED_CASES || error("Unknown bounded seed experiment")
verify_gate(joinpath(prior,"PASS512.toml"))
verify_gate(joinpath(prior,"PASS1024.toml"))
initialization=verify_seed_gate(joinpath(gate,"INITIALIZATION_PASS.toml"))
initialization["prior_pass512_sha256"]==seed_sha(joinpath(prior,"PASS512.toml")) || error("Changed reference gate")
initialization["runtime_marker_sha256"]==seed_sha(joinpath(prior,"runtime.toml")) || error("Changed repaired runtime marker")
initialization["thermal_mode_reference_sha256"]==seed_sha(joinpath(gate,"thermal_mode_reference.toml")) || error("Changed mode reference")
audit=audit_seed_source()
cfg=TOML.parsefile(seed_config_path(case))
cfg["thermal_seed_reference_file"]=abspath(joinpath(gate,"thermal_mode_reference.toml"))
cfg["thermal_seed_reference_sha256"]=initialization["thermal_mode_reference_sha256"]
out=joinpath(root,cfg["name"])
ispath(out) && error("Refusing existing case output")
configroot=joinpath(root,"configs");mkpath(configroot)
config=joinpath(configroot,case*".toml")
ispath(config) && error("Refusing existing case configuration")
open(config,"w") do io;TOML.print(io,cfg);end
run_eddy([config,out,"GPU"])
check_seed_case(out,cfg,joinpath(gate,"thermal_mode_reference.toml"))
write(joinpath(out,"source_inheritance_audit.json"),audit)
run(`tools/field-export-env/bin/python scripts/ferrel2026_v3/transfer.py $out`)
write_seed_gate(joinpath(out,"CASE_PASS.toml");case_name=case,simulated_hours=6,
    initialization_gate_sha256=seed_sha(joinpath(gate,"INITIALIZATION_PASS.toml")),
    prior_pass512_sha256=seed_sha(joinpath(prior,"PASS512.toml")),
    prior_pass1024_sha256=seed_sha(joinpath(prior,"PASS1024.toml")),
    runtime_marker_sha256=seed_sha(joinpath(prior,"runtime.toml")),
    thermal_seed_record_sha256=seed_sha(joinpath(out,"thermal_seed.toml")),
    metadata_sha256=seed_sha(joinpath(out,"metadata.toml")),
    interpretation="Six-hour horizontal thermal-seed diagnostic only; physics checks inherited explicitly, source differences audited. No convergence or developed-convection claim.")
