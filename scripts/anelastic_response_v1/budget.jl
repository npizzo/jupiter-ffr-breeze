import Breeze.AtmosphereModels as AM
# Native accounting and a separately identified centered stress decomposition.
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
@inline nb_uu(i,j,k,g,ρ,v)=ρ[i,j,k]*((v.u[i,j,k]+v.u[i+1,j,k])/2)^2
@inline nb_uv(i,j,k,g,ρ,v)=ρ[i,j,k]*((v.u[i,j,k]+v.u[i+1,j,k])/2)*((v.v[i,j,k]+v.v[i,j+1,k])/2)
@inline nb_vv(i,j,k,g,ρ,v)=ρ[i,j,k]*((v.v[i,j,k]+v.v[i,j+1,k])/2)^2
nb_reduce(f)=Field(Average(f;dims=(3,)))
nb_host(f)=(compute!(f);dropdims(Array(interior(f));dims=3))
# All accounting is for P_{0<|k|<=2} H of the centered density-weighted
# depth-mean velocity. H removes centered gradients; it is not a replacement
# for the model's native mass-flux projection. Centered-stress/WENO discrepancy
# remains a separate, measured term.
function lowproject(u,v,dx,dy)
    nx,ny=size(u);ii=[i<=nx÷2 ? i : i-nx for i=0:nx-1];jj=[j<=ny÷2 ? j : j-ny for j=0:ny-1]
    kx=reshape(sin.(2pi .* ii ./ nx) ./ dx,:,1);ky=reshape(sin.(2pi .* jj ./ ny) ./ dy,1,:)
    k2=kx .^ 2  .+  ky .^ 2;mask=(reshape(ii,:,1) .^ 2  .+  reshape(jj,1,:) .^ 2 .<=4).&(k2.>1e-24)
    a=FFTW.fft(u);b=FFTW.fft(v);dot=kx .* a .+ ky .* b;den=max.(k2,1e-24)
    real.(FFTW.ifft(mask .* (a .- kx .* dot ./ den))),real.(FFTW.ifft(mask .* (b .- ky .* dot ./ den)))
end
rdx(a,h)=(circshift(a,(-1,0)) .- circshift(a,(1,0))) ./ (2h)
rdy(a,h)=(circshift(a,(0,-1)) .- circshift(a,(0,1))) ./ (2h)
rcenter(u,v)=((u .+ circshift(u,(-1,0))) ./ 2,(v .+ circshift(v,(0,-1))) ./ 2)
mutable struct ResponseBudget
    terms::Any;total::Any;momentum::Any;products::Any;rho::Float64;dx::Float64;dy::Float64
    stride::Int;iteration::Int;t0::Float64;u0::Any;stages::Any;output::String;checks::Int
end
function ResponseBudget(m,out,c,stride)
    string(nameof(typeof(m.timestepper)))=="SSPRungeKutta3"||error("Requires SSPRK3")
    isnothing(m.timestepper.implicit_solver)||error("Explicit physics required")
    args(force)=(AM.continuity_density(m.dynamics),m.advection.momentum,m.velocities,m.closure,m.closure_fields,AM.advecting_momentum(m),m.coriolis,m.clock,Oceananigans.fields(m),force,AM.slow_dynamics(m.timestepper,m.dynamics))
    xf=(nb_xadv,nb_xcor,nb_xdif,nb_force);yf=(nb_yadv,nb_ycor,nb_ydif,nb_force)
    terms=[(nb_reduce(KernelFunctionOperation{Face,Center,Center}(xf[i],m.grid,args(m.forcing.ρu)...)),nb_reduce(KernelFunctionOperation{Center,Face,Center}(yf[i],m.grid,args(m.forcing.ρv)...))) for i=1:4]
    prod=[nb_reduce(KernelFunctionOperation{Center,Center,Center}(fun,m.grid,AM.continuity_density(m.dynamics),m.velocities)) for fun in (nb_uu,nb_uv,nb_vv)]
    path=joinpath(out,"transfer.csv")
    write(path,"time,dt,energy_start,energy_end,energy_change,mean_flux_work,convective_stress_work,native_minus_centered_flux_work,coriolis_work,viscous_work,horizontal_forcing_work,closure_residual,rhs_partition_relative,projected_velocity_residual\n")
    ResponseBudget(terms,(nb_reduce(m.timestepper.Gⁿ.ρu),nb_reduce(m.timestepper.Gⁿ.ρv)),(nb_reduce(m.momentum.ρu),nb_reduce(m.momentum.ρv)),prod,mean(Array(interior(AM.continuity_density(m.dynamics)))),c["length"]/c["nx"],c["length"]/c["ny"],stride,-1,0.,nothing,Any[],path,0)
end
function rpair(b,p)
    rcenter(nb_host(p[1]) ./ b.rho,nb_host(p[2]) ./ b.rho)
end
rproj(b,p)=lowproject(p[1],p[2],b.dx,b.dy)
function rstage(b)
    native=[rpair(b,p) for p in b.terms];total=rpair(b,b.total)
    scale=max(maximum(abs,total[1]),maximum(abs,total[2]),1e-20)
    err=max(maximum(abs,sum(p[1] for p in native) .- total[1]),maximum(abs,sum(p[2] for p in native) .- total[2]))/scale
    err<1e-8||error("Native RHS partition failed $err")
    u,v=rpair(b,b.momentum)
    xx,xy,yy=[nb_host(p) ./ b.rho for p in b.products]
    force(xx,xy,yy)=(-rdx(xx,b.dx) .- rdy(xy,b.dy),-rdx(xy,b.dx) .- rdy(yy,b.dy))
    meanforce=force(u .^ 2,u .* v,v .^ 2);conv=force(xx .- u .^ 2,xy .- u .* v,yy .- v .^ 2)
    discrepancy=(native[1][1] .- meanforce[1] .- conv[1],native[1][2] .- meanforce[2] .- conv[2])
    (;parts=[rproj(b,p) for p in [meanforce,conv,discrepancy,native[2:4]...]],err)
end
function finish_response!(b,m)
    b.iteration<0 && return
    m.clock.iteration==b.iteration+1 || return
    length(b.stages)==3||error("Missing stage in response budget")
    dt=Float64(m.clock.time)-b.t0;u1=rproj(b,rpair(b,b.momentum));u0=b.u0
    increments=[(dt .* (b.stages[1].parts[j][1] ./ 6  .+  b.stages[2].parts[j][1] ./ 6  .+  2b.stages[3].parts[j][1] ./ 3),dt .* (b.stages[1].parts[j][2] ./ 6  .+  b.stages[2].parts[j][2] ./ 6  .+  2b.stages[3].parts[j][2] ./ 3)) for j=1:6]
    expected=(sum(q[1] for q in increments),sum(q[2] for q in increments))
    ve=max(maximum(abs,u1[1] .- u0[1] .- expected[1]),maximum(abs,u1[2] .- u0[2] .- expected[2]))
    scale=max(maximum(abs,u1[1]),maximum(abs,u1[2]),1e-10)
    ve<1e-8scale+1e-13||error("Projected native budget failed $ve / $scale")
    e(a)=.5mean(a[1] .^ 2  .+  a[2] .^ 2)
    work=[mean(0.5  .*  (u1[1]  .+  u0[1])  .*  q[1]  .+  0.5  .*  (u1[2]  .+  u0[2])  .*  q[2]) for q in increments]
    res=e(u1)-e(u0)-sum(work)
    abs(res)<1e-8max(e(u1),e(u0),1e-20)+1e-20||error("Energy closure failed")
    open(b.output,"a") do io;println(io,join((Float64(m.clock.time),dt,e(u0),e(u1),e(u1)-e(u0),work...,res,maximum(s.err for s in b.stages),ve),','));end
    b.iteration=-1;b.checks+=1
end
function (b::ResponseBudget)(m)
    stage=m.clock.stage;it=m.clock.iteration
    if stage==1
        finish_response!(b,m)
        if it%b.stride==0
            b.iteration=it;b.t0=Float64(m.clock.time);b.u0=rproj(b,rpair(b,b.momentum));empty!(b.stages)
        end
    end
    if b.iteration==it
        stage==length(b.stages)+1||error("Unexpected stage order")
        push!(b.stages,rstage(b))
    end
end
Oceananigans.prognostic_state(::ResponseBudget)=nothing
Oceananigans.restore_prognostic_state!(b::ResponseBudget,::Nothing)=b
