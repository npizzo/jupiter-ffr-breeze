include("model.jl")
using Breeze: AtmosphereModels
using Oceananigans.BoundaryConditions: fill_halo_regions!
const AM=AtmosphereModels
CUDA.functional()||error("GPU required")
root=ARGS[1];mkpath(root)
cfg=TOML.parsefile("configs/kapyla-pair-v1/A1_target.toml")
small=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>32,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.2,"checkpoint_interval"=>.2,"wall_seconds"=>300.,"save_3d"=>false))
m,p=build_model(small;arch=CPU(),perturb=false)
# Paper-defined nondimensional inputs must be identities, not estimated output.
@assert isapprox(p.nu/p.chi,cfg["pr"];rtol=1e-12)
@assert isapprox((p.f/p.nu)^2,cfg["ta"];rtol=1e-12)
@assert isapprox(.1/(RGAS*reference(.5,cfg["xi"]).T)/(p.nu*p.chi),cfg["ra"];rtol=1e-12)
T=[reference(zcell(k,p),cfg["xi"]).T for i=1:p.nx,j=1:p.ny,k=1:p.nz]
Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
@assert maximum(abs,Q)<1e-12
T  .+= reshape([.03cospi(2*(i-.5)/p.nx) for i=1:p.nx],p.nx,1,1)
Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
qt=mean(2p.kc .* (T[:,:,end] .- p.tt) ./ p.dz)
@assert abs(mean(Q)*2-(p.flux-qt))<1e-12
# Manufactured compressive flow checks the missing trace stress explicitly.
set!(m;u=(x,y,z)->sinpi(x/2),v=0,w=0)
Oceananigans.TimeSteppers.update_state!(m)
F=Oceananigans.fields(m);rho=m.dynamics.dry_density
err=0.
for k=3:p.nz-2,j=1:p.ny,i=1:p.nx
    native=-AM.∂ⱼ_𝒯₁ⱼ(i,j,k,m.grid,rho,m.closure,m.closure_fields,m.clock,F,nothing)
    got=native+fx(i,j,k,m.grid,m.clock,F,p)
    want=4/3*p.nu*rho[i,j,k]*(-4sinpi(1/p.nx)^2/p.dx^2)*F.u[i,j,k]
    global err=max(err,abs(got-want))
end
@assert err<1e-12
println("PARAMETERS_CONDUCTION_DEVIATORIC_STRESS_PASS error=",err);flush(stdout)
# Static hydrostatic atmosphere with sources disabled must remain at rest.
rest=merge(small,Dict("heat"=>0.,"perturbation"=>0.,"stop_time"=>1.,"diagnostic_interval"=>1.,"movie_interval"=>1.,"checkpoint_interval"=>1.))
s,mt=run_case(rest,joinpath(root,"rest"));a=arrays(s.model)
restspeed=maximum(abs,a.w)
@assert restspeed<1e-8
println("HYDROSTATIC_PASS max_w=",restspeed);flush(stdout)
s1,mt1=run_case(small,joinpath(root,"first"))
cp=only(filter(x->endswith(x,".jld2"),readdir(joinpath(root,"first/checkpoints");join=true)))
s2,mt2=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"resumed");pickup=cp)
s3,mt3=run_case(merge(small,Dict("stop_time"=>.4)),joinpath(root,"continuous"))
a=arrays(s2.model);b=arrays(s3.model);errs=Dict{String,Float64}()
for key in (:u,:v,:w,:rho,:T)
    e=maximum(abs,getproperty(a,key) .- getproperty(b,key));scale=maximum(abs,getproperty(b,key))
    @assert e<=1e-9+1e-6scale
    errs[string(key)]=e
end
println("SERIALIZED_RESTART_PASS ",errs);flush(stdout)
# Small-grid heated transient catches source-driven acoustic failures.
pilot=merge(small,Dict("nx"=>32,"ny"=>32,"nz"=>64,"stop_time"=>20.,"diagnostic_interval"=>2.,"movie_interval"=>10.,"checkpoint_interval"=>20.,"wall_seconds"=>300.))
s,mt=run_case(pilot,joinpath(root,"heated_pilot"));@assert mt["status"]=="COMPLETE"
s=nothing;s1=nothing;s2=nothing;s3=nothing;m=nothing;a=nothing;b=nothing;GC.gc();CUDA.reclaim()
# Time the actual 256² x 128 grid, keeping this allocation bounded.
target=merge(cfg,Dict("stop_time"=>1.,"diagnostic_interval"=>.25,"movie_interval"=>1.,"checkpoint_interval"=>1.,"wall_seconds"=>600.,"save_3d"=>false))
s,mt=run_case(target,joinpath(root,"target_grid"));@assert mt["status"]=="COMPLETE"
rows=[parse.(Float64,split(l,',')) for l in readlines(joinpath(root,"target_grid/diagnostics.csv"))[2:end]]
rate=(rows[end][4]-rows[end-1][4])/(rows[end][1]-rows[end-1][1])
forecast=1.5rate*cfg["stop_time"]+300
write_toml(joinpath(root,"COST.toml"),Dict("seconds_per_simulation_unit"=>rate,"target_forecast_hours_with_margin"=>forecast/3600,"first_allocation_hours"=>8.,"gpu"=>CUDA.name(CUDA.device()),"warning"=>"Early-time timing only; turbulent CFL can reduce progress"))
# Eight hours must reach at least t=500; otherwise do not release two long jobs.
28800/(1.5rate)>=500||error("Insufficient progress within resource cap")
write_toml(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","model_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl")))),"restart_errors"=>errs,"rest_max_w"=>restspeed,"stress_error"=>err,"grid"=>[256,256,128],"seconds_per_simulation_unit"=>rate,"gpu"=>CUDA.name(CUDA.device()),"scope"=>"Implementation and initial timing only, not a physical-reproduction or vorticity-budget-closure certificate"))
println("KAPYLA_PAIR_VALIDATION_PASS");flush(stdout)
