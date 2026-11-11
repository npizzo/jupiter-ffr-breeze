include("run.jl")
include("check_ferrel.jl")
include("common.jl")
length(ARGS)==3 || error("usage: campaign.jl CASE ROOT PREFLIGHT_ROOT")
case,root,preflight=ARGS
case in ("b0","b162","boundary162") || error("Only authorized two-day cases")
verify_gate(joinpath(preflight,"PASS512.toml"))
cfg=TOML.parsefile(config_path(case))
cfg["stop_seconds"]==172800.0 || error("First screening is bounded at two days")
name=cfg["name"];out=joinpath(root,name);mkpath(joinpath(root,"configs"))
config=joinpath(root,"configs","$(name).toml")
isfile(config) && error("Refusing duplicate case")
open(config,"w") do io;TOML.print(io,cfg);end
run_eddy([config,out,"GPU"]);check_ferrel_output(out)
