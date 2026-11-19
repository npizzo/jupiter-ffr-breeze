include("model.jl")
cfg=TOML.parsefile("configs/cai-benchmark-v1/B7_anticyclone.toml")
small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>16,"stop_time"=>.1,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.1,"wall_seconds"=>120.))
m,p=build_model(small;arch=CPU(),perturb=false)
println("CPU_MODEL_BUILD_PASS fields=",keys(Oceananigans.fields(m)));flush(stdout)
T=[5-4((k-.5)*p.dz) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
@assert maximum(abs,[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz])<1e-12
sim,meta=run_case(small,ARGS[1];arch=CPU())
@assert meta["status"]=="COMPLETE"
println("CPU_REFERENCE_PASS")
