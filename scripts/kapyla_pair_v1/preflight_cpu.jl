include("model.jl")
using Breeze: AtmosphereModels
using Oceananigans.BoundaryConditions: fill_halo_regions!
const AM=AtmosphereModels
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
rest=merge(small,Dict("heat"=>0.,"perturbation"=>0.,"stop_time"=>.2,"diagnostic_interval"=>.1))
s,mt=run_case(rest,joinpath(root,"rest");arch=CPU())
a=arrays(s.model);@assert maximum(abs,a.w)<1e-8
s,mt=run_case(small,joinpath(root,"heated");arch=CPU())
println("CPU_PREFLIGHT_PASS");flush(stdout)
