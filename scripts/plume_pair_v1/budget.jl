using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans: TendencyCallsite
using Oceananigans.Operators: V⁻¹ᶠᶜᶜ, V⁻¹ᶜᶠᶜ, δzᵃᵃᶜ, ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.Advection: _advective_momentum_flux_Wu, _advective_momentum_flux_Wv
using Oceananigans.BoundaryConditions: fill_halo_regions!

@inline pxadv(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.x_momentum_flux_divergence(i,j,k,g,a,M,v,dyn)
@inline pyadv(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.y_momentum_flux_divergence(i,j,k,g,a,M,v,dyn)
@inline pxvert(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-V⁻¹ᶠᶜᶜ(i,j,k,g)*δzᵃᵃᶜ(i,j,k,g,_advective_momentum_flux_Wu,a,M[3],v.u)
@inline pyvert(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-V⁻¹ᶜᶠᶜ(i,j,k,g)*δzᵃᵃᶜ(i,j,k,g,_advective_momentum_flux_Wv,a,M[3],v.v)
@inline pxcor(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.x_f_cross_U(i,j,k,g,f,M)
@inline pycor(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.y_f_cross_U(i,j,k,g,f,M)
@inline pxdif(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,ρ,c,cf,t,F,nothing)
@inline pydif(i,j,k,g,ρ,a,v,c,cf,M,f,t,F,dyn)=-AM.∂ⱼ_𝒯₂ⱼ(i,j,k,g,ρ,c,cf,t,F,nothing)
@inline pxp(i,j,k,g,ρ,p)=-ρ[i,j,k]*∂xᶠᶜᶜ(i,j,k,g,p)
@inline pyp(i,j,k,g,ρ,p)=-ρ[i,j,k]*∂yᶜᶠᶜ(i,j,k,g,p)
@inline masked_op(i,j,k,g,lo,hi,fun,args...) = lo <= k <= hi ? fun(i,j,k,g,args...) : zero(eltype(g))
@inline masked_field(i,j,k,g,lo,hi,a) = lo <= k <= hi ? a[i,j,k] : zero(eltype(g))
host2(f)=(compute!(f);dropdims(Array(interior(f));dims=3))
dx(a,h)=(circshift(a,(-1,0)).-circshift(a,(1,0)))./(2h)
dy(a,h)=(circshift(a,(0,-1)).-circshift(a,(0,1)))./(2h)
pcurl(u,v,h)=dx((v.+circshift(v,(0,-1)))./2,h).-dy((u.+circshift(u,(-1,0)))./2,h)
const PLUME_TERMS=["horizontal_advection","vertical_advection_interface_transport","coriolis","viscosity"]

mutable struct PlumeBudget
    terms::Any; totals::Any; momenta::Any; pressures::Any; pressure_field::Any
    norms::Vector{Float64}; h::Float64; origin::Any; cumulative::Any
    iteration::Int; t0::Float64; z0::Any; stages::Any; stage_times::Vector{Float64}
    checks::Int; max_residual::Float64; max_partition::Float64; max_pressure::Float64
    path::String; layer_indices::Any
end
function PlumeBudget(m,out,rho,z,h)
    isnothing(m.timestepper.implicit_solver)||error("Explicit SSPRK3 required")
    string(nameof(typeof(m.timestepper)))=="SSPRungeKutta3"||error("SSPRK3 required")
    bounds=((70000.,80000.),(80000.,84000.),(84000.,110000.))
    ranges=[findall(a .<= z .< b) for (a,b) in bounds]
    all(!isempty(r) for r in ranges)||error("Empty budget layer")
    args=(AM.continuity_density(m.dynamics),m.advection.momentum,m.velocities,m.closure,m.closure_fields,AM.advecting_momentum(m),m.coriolis,m.clock,Oceananigans.fields(m),AM.slow_dynamics(m.timestepper,m.dynamics))
    reduce_op(fun,loc,r,args...)=Field(Average(KernelFunctionOperation{loc...}(masked_op,m.grid,first(r),last(r),fun,args...);dims=(3,)))
    reduce_field(a,loc,r)=Field(Average(KernelFunctionOperation{loc...}(masked_field,m.grid,first(r),last(r),a);dims=(3,)))
    xf=(pxadv,pxvert,pxcor,pxdif);yf=(pyadv,pyvert,pycor,pydif)
    terms=[[(reduce_op(xf[i],(Face,Center,Center),r,args...),reduce_op(yf[i],(Center,Face,Center),r,args...)) for i in 1:4] for r in ranges]
    pairs(x,y)=[(reduce_field(x,(Face,Center,Center),r),reduce_field(y,(Center,Face,Center),r)) for r in ranges]
    pressures=[(reduce_op(pxp,(Face,Center,Center),r,args[1],m.dynamics.pressure_anomaly),reduce_op(pyp,(Center,Face,Center),r,args[1],m.dynamics.pressure_anomaly)) for r in ranges]
    norms=[sum(rho[r])/length(rho) for r in ranges]
    path=joinpath(out,"budget");mkpath(path)
    open(joinpath(path,"checks.csv"),"w") do io;println(io,"iteration,time,absolute_step_residual,relative_step_residual,rhs_partition,pressure_curl");end
    PlumeBudget(terms,pairs(m.timestepper.Gⁿ.ρu,m.timestepper.Gⁿ.ρv),pairs(m.momentum.ρu,m.momentum.ρv),pressures,m.dynamics.pressure_anomaly,norms,h,nothing,nothing,-1,0.,nothing,Any[nothing,nothing,nothing],fill(NaN,3),0,0.,0.,0.,path,ranges)
end
paircurl(b,p,n)=pcurl(host2(p[1]),host2(p[2]),b.h)./n
budget_z(b)=cat([paircurl(b,p,n) for (p,n) in zip(b.momenta,b.norms)]...;dims=3)
function budget_stage(b)
    parts=cat([cat([paircurl(b,p,n) for p in ps]...;dims=3) for (ps,n) in zip(b.terms,b.norms)]...;dims=4)
    # First operation is total advection; split it without discarding vertical flux.
    parts[:,:,1,:].-=parts[:,:,2,:]
    total=cat([paircurl(b,p,n) for (p,n) in zip(b.totals,b.norms)]...;dims=3)
    pred=dropdims(sum(parts;dims=3);dims=3)
    err=maximum(abs,total.-pred);scale=max(maximum(abs,total),maximum(abs,parts),1e-25)
    err<=1e-10scale+1e-20||error("Layer RHS partition failed: $err / $scale")
    fill_halo_regions!(b.pressure_field)
    pg=maximum(maximum(abs,paircurl(b,p,n)) for (p,n) in zip(b.pressures,b.norms))
    pg<=1e-15||error("Layer pressure curl failed: $pg")
    b.max_partition=max(b.max_partition,err);b.max_pressure=max(b.max_pressure,pg)
    parts
end
function finish_step!(b,m)
    all(s->s!==nothing,b.stages)||error("Missing SSPRK stage")
    dt=Float64(m.clock.time)-b.t0
    isapprox(b.stage_times[2],b.t0+dt;atol=1e-7)&&isapprox(b.stage_times[3],b.t0+dt/2;atol=1e-7)||error("RK stage times")
    dz=dt.*(b.stages[1]./6 .+ b.stages[2]./6 .+ 2b.stages[3]./3)
    z=budget_z(b);pred=dropdims(sum(dz;dims=3);dims=3)
    err=maximum(abs,z.-b.z0.-pred);scale=max(maximum(abs,z.-b.z0),maximum(abs,dz),1e-25)
    allowance=8192eps(Float64)*max(maximum(abs,z),maximum(abs,b.z0),1e-8)
    err<=allowance+1e-8scale||error("Layer step closure failed: $err / $scale")
    b.cumulative.+=dz;b.checks+=1;b.max_residual=max(b.max_residual,err)
    if b.checks%30==0
        open(joinpath(b.path,"checks.csv"),"a") do io;println(io,join((m.clock.iteration,m.clock.time,err,err/max(scale,allowance),b.max_partition,b.max_pressure),','));end
    end
    b.iteration=-1
end
function (b::PlumeBudget)(m)
    it=m.clock.iteration;s=m.clock.stage
    if s==1
        if b.iteration>=0;it==b.iteration+1||error("Skipped budget step");finish_step!(b,m);end
        b.iteration=it;b.t0=Float64(m.clock.time);b.z0=budget_z(b)
        if b.origin===nothing;b.origin=copy(b.z0);b.cumulative=zeros(size(b.z0,1),size(b.z0,2),4,3);end
        b.stages=Any[nothing,nothing,nothing];b.stage_times.=NaN
    end
    b.stages[s]=budget_stage(b);b.stage_times[s]=Float64(m.clock.time)
    nothing
end
Oceananigans.prognostic_state(::PlumeBudget)=nothing
Oceananigans.restore_prognostic_state!(b::PlumeBudget,::Nothing)=b
function save_budget!(b,m)
    if b.iteration>=0 && m.clock.iteration==b.iteration+1;finish_step!(b,m);end
    z=budget_z(b);origin=b.origin===nothing ? z : b.origin
    parts=b.cumulative===nothing ? zeros(size(z,1),size(z,2),4,3) : b.cumulative
    jldsave(joinpath(b.path,@sprintf("budget_%06d.jld2",round(Int,m.clock.time)));t=Float64(m.clock.time),z,origin,cumulative=parts,residual=z.-origin.-dropdims(sum(parts;dims=3);dims=3),terms=PLUME_TERMS,norms=b.norms,checks=b.checks,max_step_absolute_residual=b.max_residual,max_partition=b.max_partition,max_pressure=b.max_pressure,layer_indices=b.layer_indices)
end
