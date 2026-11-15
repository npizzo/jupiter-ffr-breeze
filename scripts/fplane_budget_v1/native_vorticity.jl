# Native explicit momentum terms, reduced vertically on GPU before host transfer.
# Linear curl matches the existing cell-centered diagnostic Z exactly.
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans: TendencyCallsite
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.BoundaryConditions: fill_halo_regions!
@inline nb_xadv(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.x_momentum_flux_divergence(i,j,k,g,a,M,v,dyn)
@inline nb_yadv(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.y_momentum_flux_divergence(i,j,k,g,a,M,v,dyn)
@inline nb_xcor(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.x_f_cross_U(i,j,k,g,f,M)
@inline nb_ycor(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.y_f_cross_U(i,j,k,g,f,M)
@inline nb_xdif(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,ρ,c,cf,t,F,nothing)
@inline nb_ydif(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=-AM.∂ⱼ_𝒯₂ⱼ(i,j,k,g,ρ,c,cf,t,F,nothing)
@inline nb_force(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,force,dyn)=force(i,j,k,g,t,F)
@inline nb_xp(i,j,k,g,ρ,p)=-ρ[i,j,k]*∂xᶠᶜᶜ(i,j,k,g,p)
@inline nb_yp(i,j,k,g,ρ,p)=-ρ[i,j,k]*∂yᶜᶠᶜ(i,j,k,g,p)
nb_reduce(f)=Field(Average(f;dims=(3,)))
nb_host(f)=(compute!(f);dropdims(Array(interior(f));dims=3))
function nb_curl(uf,vf,dx,dy)
 u=(uf.+circshift(uf,(-1,0)))./2;v=(vf[:,1:end-1].+vf[:,2:end])./2
 py=hcat(u[:,1],u,u[:,end])
 (circshift(v,(-1,0)).-circshift(v,(1,0)))./(2dx) .- (py[:,3:end].-py[:,1:end-2])./(2dy)
end
mutable struct NativeVorticityBudget
 terms::Any; total::Any; momentum::Any; pressure::Any; pressure_field::Any
 rho_mean::Float64; dx::Float64;dy::Float64; mask::Any;output::String
 iteration::Int; t0::Float64; z0::Any; stages::Any; stage_times::Any
 origin::Any; origin_time::Float64; cumulative::Any; third::Any;checks::Int;max_relative::Float64;max_partition::Float64
end
function NativeVorticityBudget(m,out,cfg,rho)
 isnothing(m.timestepper.implicit_solver)||error("Budget requires explicit physics")
 string(nameof(typeof(m.timestepper)))=="SSPRungeKutta3"||error("Budget requires SSPRK3")
 cfg["rotation_model"]=="f-plane"||error("This audit is f-plane only")
 cfg["jet_speed"]==0||error("No imposed jet forcing")
 string(nameof(typeof(m.closure)))=="ScalarDiffusivity"||error("Unaudited closure")
 args(force)=(AM.continuity_density(m.dynamics),m.advection.momentum,m.velocities,m.closure,m.closure_fields,AM.advecting_momentum(m),m.coriolis,m.clock,Oceananigans.fields(m),force,AM.slow_dynamics(m.timestepper,m.dynamics))
 xf=(nb_xadv,nb_xcor,nb_xdif,nb_force);yf=(nb_yadv,nb_ycor,nb_ydif,nb_force)
 terms=[(nb_reduce(KernelFunctionOperation{Face,Center,Center}(xf[i],m.grid,args(m.forcing.ρu)...)),nb_reduce(KernelFunctionOperation{Center,Face,Center}(yf[i],m.grid,args(m.forcing.ρv)...))) for i in 1:4]
 pressure=(nb_reduce(KernelFunctionOperation{Face,Center,Center}(nb_xp,m.grid,args(m.forcing.ρu)[1],m.dynamics.pressure_anomaly)),nb_reduce(KernelFunctionOperation{Center,Face,Center}(nb_yp,m.grid,args(m.forcing.ρv)[1],m.dynamics.pressure_anomaly)))
 y=((1:cfg["ny"]).-.5).*cfg["ly"]./cfg["ny"].-cfg["ly"]/2;mask=repeat(reshape(abs.(y).<=.4cfg["ly"],1,:),cfg["nx"],1)
 path=joinpath(out,"native_vorticity");mkpath(path)
 open(joinpath(path,"checks.csv"),"w") do io;println(io,"iteration,start_time,end_time,dt,absolute_residual,relative_residual,rhs_partition_residual,pressure_curl_max");end
 NativeVorticityBudget(terms,(nb_reduce(m.timestepper.Gⁿ.ρu),nb_reduce(m.timestepper.Gⁿ.ρv)),(nb_reduce(m.momentum.ρu),nb_reduce(m.momentum.ρv)),pressure,m.dynamics.pressure_anomaly,mean(rho),cfg["lx"]/cfg["nx"],cfg["ly"]/cfg["ny"],mask,path,-1,0.,nothing,Any[nothing,nothing,nothing],fill(NaN,3),nothing,NaN,nothing,zeros(2,4),0,0.,0.)
end
function nb_pair(b,pair;normal_zero=true)
 u=nb_host(pair[1])./b.rho_mean;v=nb_host(pair[2])./b.rho_mean
 normal_zero && (v[:,[1,end]].=0)
 nb_curl(u,v,b.dx,b.dy)
end
function nb_terms(b)
 parts=cat([nb_pair(b,p) for p in b.terms]...;dims=3);total=nb_pair(b,b.total)
 residual=maximum(abs,(dropdims(sum(parts;dims=3);dims=3).-total)[b.mask]);scale=max(maximum(abs,total[b.mask]),maximum(abs,parts),1e-25)
 residual<1e-10scale+1e-20||error("Native vorticity RHS partition failed: $residual / $scale")
 fill_halo_regions!(b.pressure_field);pg=maximum(abs,nb_pair(b,b.pressure)[b.mask])
 pg<1e-15||error("Pressure curl failed to cancel: $pg")
 (;parts,residual,pg)
end
function nb_finish!(b,m)
 all(s->s!==nothing,b.stages)||error("Missing RK stage")
 dt=Float64(m.clock.time)-b.t0;dt>0||error("Nonpositive sampled step")
 isapprox(b.stage_times[2],b.t0+dt;atol=1e-7)||error("RK stage2 time")
 isapprox(b.stage_times[3],b.t0+dt/2;atol=1e-7)||error("RK stage3 time")
 z=nb_pair(b,b.momentum);parts=dt.*(b.stages[1].parts./6 .+ b.stages[2].parts./6 .+ 2b.stages[3].parts./3)
 observed=z.-b.z0;pred=dropdims(sum(parts;dims=3);dims=3);err=maximum(abs,(observed.-pred)[b.mask]);scale=max(maximum(abs,observed[b.mask]),maximum(abs,parts),1e-25)
 allowance=8192eps(Float64)*max(maximum(abs,z),maximum(abs,b.z0),1e-8)
 err<=allowance+1e-8scale||error("Native vorticity step closure failed: $err vs $scale allowance $allowance")
 rel=err/max(scale,allowance)
 # Exact finite-increment decomposition of Z^3, not sampled 3Z² dZ/dt.
 weight=z.^2 .+ z.*b.z0 .+ b.z0.^2
 for i in 1:4
  q=weight.*parts[:,:,i];b.third[1,i]+=mean(q);b.third[2,i]+=mean(q[b.mask])
 end
 b.cumulative .+= parts;b.checks+=1;b.max_relative=max(b.max_relative,rel);b.max_partition=max(b.max_partition,maximum(s.residual for s in b.stages))
 open(joinpath(b.output,"checks.csv"),"a") do io;println(io,join((b.iteration,b.t0,m.clock.time,dt,err,rel,maximum(s.residual for s in b.stages),maximum(s.pg for s in b.stages)),','));end
 b.iteration=-1;nothing
end
function (b::NativeVorticityBudget)(m)
 it=m.clock.iteration;stage=m.clock.stage
 if stage==1
  if b.iteration>=0 && it==b.iteration+1;nb_finish!(b,m)
  elseif b.iteration>=0 && it!=b.iteration;error("Skipped native vorticity step");end
  b.iteration=it;b.t0=Float64(m.clock.time);b.z0=nb_pair(b,b.momentum)
  if b.origin===nothing;b.origin=copy(b.z0);b.origin_time=b.t0;b.cumulative=zeros(size(b.z0)...,4);end
  b.stages=Any[nothing,nothing,nothing];b.stage_times.=NaN
 end
 if it==b.iteration
  stage in 1:3||error("Unknown RK stage");b.stages[stage]=nb_terms(b);b.stage_times[stage]=Float64(m.clock.time)
 end
 nothing
end
Oceananigans.prognostic_state(::NativeVorticityBudget)=nothing
Oceananigans.restore_prognostic_state!(b::NativeVorticityBudget,::Nothing)=b
function save_native_vorticity!(b,m)
 if b.iteration>=0 && m.clock.iteration==b.iteration+1;nb_finish!(b,m)
 elseif b.iteration>=0
  m.clock.iteration==b.iteration && Float64(m.clock.time)==b.t0 && b.stages[2]===nothing && b.stages[3]===nothing||error("Incomplete RK step at output")
 end
 z=nb_pair(b,b.momentum);origin=b.origin===nothing ? z : b.origin;parts=b.cumulative===nothing ? zeros(size(z)...,4) : b.cumulative
 residual=z.-origin.-dropdims(sum(parts;dims=3);dims=3)
 thirds=[mean(z.^3 .- origin.^3),mean((z.^3 .- origin.^3)[b.mask])]
 third_residual=thirds.-vec(sum(b.third;dims=2))
 jldsave(joinpath(b.output,@sprintf("budget_%012d.jld2",round(Int,m.clock.time)));t=Float64(m.clock.time),origin_time=b.origin_time,Z=z,origin_Z=origin,cumulative_parts=parts,cumulative_residual=residual,third_moment_storage=thirds,third_moment_parts=copy(b.third),third_moment_residual=third_residual,mask=b.mask,checks=b.checks,max_relative_step_residual=b.max_relative,max_rhs_partition_residual=b.max_partition,term_names=["native_advection","native_coriolis","native_diffusion","native_external_forcing"])
 nothing
end
