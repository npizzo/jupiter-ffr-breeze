include("model.jl")
CUDA.functional()||error("GPU required")
root=ARGS[1];mkpath(root)
cfg=TOML.parsefile("configs/cai-benchmark-v1/B7_anticyclone.toml")
small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>16,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"wall_seconds"=>300.))
m,p=build_model(small;arch=CPU(),perturb=false)
@assert :T in keys(Oceananigans.fields(m)) && :ρᵈ in keys(Oceananigans.fields(m))
# Exact finite-volume source balance and manufactured profile.
T=[5-4((k-.5)*p.dz) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
@assert maximum(abs,Q)<1e-12
T .+=reshape([.03cos(2pi*(i-.5)/p.nx) for i=1:p.nx],p.nx,1,1)
Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
qt=mean(2p.conductivity .* (T[:,:,end] .- 1) ./ p.dz)
@assert abs(mean(Q)*1-(p.flux-qt))<1e-12
for k in 2:p.nz-1,j in 1:p.ny,i in 1:p.nx
 expected=-p.conductivity*4sinpi(1/p.nx)^2/p.dx^2*.03cos(2pi*(i-.5)/p.nx)
 @assert abs(Q[i,j,k]-expected)<1e-12
end
@assert m.coriolis.rotation.fy≈0 atol=1e-12
@assert m.coriolis.rotation.fz==2cfg["omega"]
a=state_arrays(m)
@assert maximum(abs,a.u)==0 && maximum(abs,a.v)==0 && maximum(abs,a.w)==0
@assert maximum(abs,a.T .- reshape([5-4((k-.5)*p.dz) for k=1:p.nz],1,1,:))<.01
# Verify the momentum-name adapter against an analytic cross product at both latitudes.
M=(ρu=fill(1.,16,16,16),ρv=fill(2.,16,16,16),ρw=fill(3.,16,16,16))
for latitude in (90.,22.5)
 c=MomentumCoriolis(ConstantCartesianCoriolis(;rotation_rate=.64,latitude))
 got=(Oceananigans.Coriolis.x_f_cross_U(8,8,8,m.grid,c,M),Oceananigans.Coriolis.y_f_cross_U(8,8,8,m.grid,c,M),Oceananigans.Coriolis.z_f_cross_U(8,8,8,m.grid,c,M))
 want=(3c.rotation.fy-2c.rotation.fz,c.rotation.fz,-c.rotation.fy)
 @assert maximum(abs,collect(got).-collect(want))<1e-12
 @assert abs(got[1]+2got[2]+3got[3])<1e-12
end
println("CONDUCTION_INITIALIZATION_CORIOLIS_PASS");flush(stdout)
# Same GPU calculation run continuously and via an actual serialized checkpoint.
s1,meta1=run_case(small,joinpath(root,"segment1"))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"segment1/checkpoints");join=true)))
s2,meta2=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"resumed");pickup=cp)
s3,meta3=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"continuous"))
a=state_arrays(s2.model);b=state_arrays(s3.model)
errors=Dict{String,Float64}()
for k in (:u,:v,:w,:rho,:T)
 err=maximum(abs,getproperty(a,k) .- getproperty(b,k));scale=maximum(abs,getproperty(b,k))
 @assert err<=1e-9+1e-6scale
 errors[string(k)]=err
end
println("SERIALIZED_RESTART_PASS ",errors);flush(stdout)
s1=nothing;s2=nothing;s3=nothing;m=nothing;GC.gc();CUDA.reclaim()
# Actual target-grid smoke test and timing; no coarse-grid runtime extrapolation.
target=merge(cfg,Dict("stop_time"=>2.,"diagnostic_interval"=>.5,"movie_interval"=>1.,"checkpoint_interval"=>2.,"wall_seconds"=>400.))
s,meta=run_case(target,joinpath(root,"target_grid"))
@assert meta["status"]=="COMPLETE"
# Cost check uses actual production-grid stepping after compilation. It does
# not prove late-time CFL or physical resolution; production retains its wall cap.
rows=[parse.(Float64,split(l,",")) for l in readlines(joinpath(root,"target_grid/diagnostics.csv"))[2:end]]
rate=(rows[end][4]-rows[end-1][4])/(rows[end][1]-rows[end-1][1])
forecast=1.25rate*cfg["stop_time"]+300
write_toml(joinpath(root,"COST.toml"),Dict("seconds_per_simulation_unit"=>rate,"forecast_seconds_with_25percent_margin_plus_startup"=>forecast,"allocation_limit_seconds"=>14400,"gpu"=>CUDA.name(CUDA.device())))
forecast<14400 || error("Target-grid forecast exceeds four-hour per-case budget; retain validation, do not release production")
write_toml(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","gpu"=>CUDA.name(CUDA.device()),"restart_errors"=>errors,"conduction_operator"=>"manufactured solution and global boundary-flux balance PASS","target_grid"=>[128,128,64],"target_grid_final_time"=>meta["final_time"],"target_grid_wall_seconds"=>meta["wall_elapsed"],"model_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl")))),"scope"=>"implementation checks only; no physical-fidelity claim"))
println("CAI_BENCHMARK_VALIDATION_PASS");flush(stdout)
