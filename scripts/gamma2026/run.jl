include("core.jl")
using .Gamma2026, TOML
length(ARGS)>=3 || error("CONFIG OUTPUT GPU_GATE [resume]")
c=TOML.parsefile(ARGS[1])
Gamma2026.run(c,ARGS[2];gate=ARGS[3],resume=length(ARGS)>3 && ARGS[4]=="resume")
