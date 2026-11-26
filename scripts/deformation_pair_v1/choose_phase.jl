# Bounded initialization-only search. Never integrate or evaluate vortex preference.
include("model.jl")
c=TOML.parsefile(joinpath(@__DIR__,"config.toml"))
out=ARGS[1];@assert !ispath(out);mkpath(out)
for seed=261128:261147
    cc=merge(c,Dict("seed"=>seed))
    try
        m,p,info=build(cc,2.;arch=CPU())
        write_toml(joinpath(out,"PHASE_PASS.toml"),Dict("seed"=>seed,"status"=>"PASS","initial_min_N2_over_background"=>info["minimum_initial_N2"]/info["N2"],"criterion"=>"first deterministic phase realization with min N2 > 10% background; no evolution or sign-outcome selection"))
        println("PHASE_PASS seed=$seed margin=$(info["minimum_initial_N2"]/info["N2"])");flush(stdout)
        exit()
    catch e
        msg=sprint(showerror,e)
        occursin("INITIAL_STABILITY_MARGIN",msg)||rethrow()
        println("PHASE_REJECT seed=$seed $msg");flush(stdout)
    end
    GC.gc()
end
error("No initial phase realization passed the predefined static-stability margin")
