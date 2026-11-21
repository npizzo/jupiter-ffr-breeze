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
nb_dx(a,h)=(circshift(a,(-1,0)).-circshift(a,(1,0)))./(2h)
function nb_dy(a,h,parity=1)
    p=hcat(parity.*a[:,1],a,parity.*a[:,end])
    (p[:,3:end].-p[:,1:end-2])./(2h)
end
nb_ccurl(u,v,dx,dy)=nb_dx(v,dx).-nb_dy(u,dy)
nb_center(u,v)=((u.+circshift(u,(-1,0)))./2,(v[:,1:end-1].+v[:,2:end])./2)
function nb_curl(uf,vf,dx,dy)
    u,v=nb_center(uf,vf);nb_ccurl(u,v,dx,dy)
end
function nb_stress(xx,xy,yy,dx,dy)
    fx=-nb_dx(xx,dx).-nb_dy(xy,dy,-1)
    fy=-nb_dx(xy,dx).-nb_dy(yy,dy)
    nb_ccurl(fx,fy,dx,dy)
end
const SIGN_TERMS=["mean_transport","mean_flux_minus_transport","convective_stress","native_advection_minus_centered_flux","native_coriolis","native_diffusion","native_horizontal_forcing"]
mutable struct NativeVorticityBudget
    terms::Any;total::Any;momentum::Any;pressure::Any;pressure_field::Any;products::Any
    rho_mean::Float64;dx::Float64;dy::Float64;f::Float64;masks::Any;output::String
    iteration::Int;t0::Float64;z0::Any;stages::Any;stage_times::Any
    origin::Any;origin_time::Float64;cumulative::Any;third::Any;central_third::Any
    bin_edges::Any;conditional_sum::Any;conditional_count::Any
    checks::Int;max_relative::Float64;max_partition::Float64;max_pressure::Float64;max_mass_div::Float64
end
function NativeVorticityBudget(m,out,cfg,rho)
    isnothing(m.timestepper.implicit_solver)||error("Budget requires explicit physics")
    string(nameof(typeof(m.timestepper)))=="SSPRungeKutta3"||error("Requires SSPRK3")
    cfg["rotation_model"]=="f-plane" && cfg["beta"]==0 && cfg["jet_speed"]==0 || error("Requires f-plane, no jet")
    string(nameof(typeof(m.closure)))=="ScalarDiffusivity"||error("Unaudited closure")
    args(force)=(AM.continuity_density(m.dynamics),m.advection.momentum,m.velocities,m.closure,m.closure_fields,AM.advecting_momentum(m),m.coriolis,m.clock,Oceananigans.fields(m),force,AM.slow_dynamics(m.timestepper,m.dynamics))
    xf=(nb_xadv,nb_xcor,nb_xdif,nb_force);yf=(nb_yadv,nb_ycor,nb_ydif,nb_force)
    terms=[(nb_reduce(KernelFunctionOperation{Face,Center,Center}(xf[i],m.grid,args(m.forcing.ρu)...)),nb_reduce(KernelFunctionOperation{Center,Face,Center}(yf[i],m.grid,args(m.forcing.ρv)...))) for i in 1:4]
    pressure=(nb_reduce(KernelFunctionOperation{Face,Center,Center}(nb_xp,m.grid,args(m.forcing.ρu)[1],m.dynamics.pressure_anomaly)),nb_reduce(KernelFunctionOperation{Center,Face,Center}(nb_yp,m.grid,args(m.forcing.ρv)[1],m.dynamics.pressure_anomaly)))
    products=[nb_reduce(KernelFunctionOperation{Center,Center,Center}(fun,m.grid,AM.continuity_density(m.dynamics),m.velocities)) for fun in (nb_uu,nb_uv,nb_vv)]
    y=((1:cfg["ny"]).-.5)./cfg["ny"].-.5
    interior_mask=repeat(reshape(abs.(y).<=.4,1,:),cfg["nx"],1)
    masks=[trues(cfg["nx"],cfg["ny"]),interior_mask,.!interior_mask]
    path=joinpath(out,"native_vorticity");mkpath(path)
    open(joinpath(path,"checks.csv"),"w") do io
        println(io,"iteration,start_time,end_time,dt,absolute_residual_full,relative_residual_full,rhs_partition_residual_full,pressure_curl_full,mass_divergence_full")
    end
    edges=vcat(-Inf,-reverse([.002,.005,.01,.02,.05,.1,.2,.5,1.,2.,5.]),0.,[.002,.005,.01,.02,.05,.1,.2,.5,1.,2.,5.],Inf)
    NativeVorticityBudget(terms,(nb_reduce(m.timestepper.Gⁿ.ρu),nb_reduce(m.timestepper.Gⁿ.ρv)),(nb_reduce(m.momentum.ρu),nb_reduce(m.momentum.ρv)),pressure,m.dynamics.pressure_anomaly,products,mean(rho),cfg["lx"]/cfg["nx"],cfg["ly"]/cfg["ny"],get(cfg,"f_sign",1)*2cfg["rotation_rate"]*sind(cfg["latitude"]),masks,path,-1,0.,nothing,Any[nothing,nothing,nothing],fill(NaN,3),nothing,NaN,nothing,zeros(3,7),zeros(3,7),edges,zeros(length(edges)-1,7,3),zeros(length(edges)-1,3),0,0.,0.,0.,0.)
end
function nb_faces(b,pair;normal_zero=true)
    u=nb_host(pair[1])./b.rho_mean;v=nb_host(pair[2])./b.rho_mean
    normal_zero && (v[:,[1,end]].=0)
    u,v
end
function nb_pair(b,pair)
    u,v=nb_faces(b,pair);nb_curl(u,v,b.dx,b.dy)
end
function nb_terms(b)
    native=cat([nb_pair(b,p) for p in b.terms]...;dims=3)
    total=nb_pair(b,b.total)
    residual=maximum(abs,dropdims(sum(native;dims=3);dims=3).-total)
    scale=max(maximum(abs,total),maximum(abs,native),1e-25)
    residual<1e-10scale+1e-20||error("Native RHS partition failed full domain: $residual / $scale")
    fill_halo_regions!(b.pressure_field)
    # Do not zero pressure-gradient normal faces: test the actual gradient curl.
    pu,pv=nb_faces(b,b.pressure;normal_zero=false)
    pg=maximum(abs,nb_curl(pu,pv,b.dx,b.dy))
    pg<1e-15||error("Full-domain pressure curl failed: $pg")
    uf,vf=nb_faces(b,b.momentum)
    massdiv=maximum(abs,(circshift(uf,(-1,0)).-uf)./b.dx.+diff(vf;dims=2)./b.dy)
    massdiv<1e-10||error("Depth mean divergence: $massdiv")
    u,v=nb_center(uf,vf);z=nb_ccurl(u,v,b.dx,b.dy)
    rxx=nb_host(b.products[1])./b.rho_mean.-u.^2
    rxy=nb_host(b.products[2])./b.rho_mean.-u.*v
    ryy=nb_host(b.products[3])./b.rho_mean.-v.^2
    conv=nb_stress(rxx,rxy,ryy,b.dx,b.dy)
    meanflux=nb_stress(u.^2,u.*v,v.^2,b.dx,b.dy)
    transport=-u.*nb_dx(z,b.dx).-v.*nb_dy(z,b.dy,-1)
    discrepancy=native[:,:,1].-meanflux.-conv
    parts=cat(transport,meanflux.-transport,conv,discrepancy,native[:,:,2:4];dims=3)
    (;parts,z,residual,pg,massdiv)
end
function nb_finish!(b,m)
    all(s->s!==nothing,b.stages)||error("Missing RK stage")
    dt=Float64(m.clock.time)-b.t0;dt>0||error("Nonpositive step")
    isapprox(b.stage_times[2],b.t0+dt;atol=1e-7)||error("RK stage 2 time")
    isapprox(b.stage_times[3],b.t0+dt/2;atol=1e-7)||error("RK stage 3 time")
    z=nb_pair(b,b.momentum)
    parts=dt.*(b.stages[1].parts./6 .+ b.stages[2].parts./6 .+ 2b.stages[3].parts./3)
    observed=z.-b.z0;pred=dropdims(sum(parts;dims=3);dims=3)
    err=maximum(abs,observed.-pred);scale=max(maximum(abs,observed),maximum(abs,parts),1e-25)
    allowance=8192eps(Float64)*max(maximum(abs,z),maximum(abs,b.z0),1e-8)
    err<=allowance+1e-8scale||error("Native full-domain closure failed: $err vs $scale allowance $allowance")
    weight=z.^2 .+ z.*b.z0 .+ b.z0.^2
    for (r,mask) in enumerate(b.masks)
        zc=z[mask].-mean(z[mask]);zc0=b.z0[mask].-mean(b.z0[mask]);cw=zc.^2 .+ zc.*zc0 .+ zc0.^2
        for k in 1:7
            q=parts[:,:,k][mask]
            b.third[r,k]+=mean(weight[mask].*q)
            b.central_third[r,k]+=mean(cw.*(q.-mean(q)))
        end
    end
    # Stage quadrature for conditional drift; no sparse-output reconstruction.
    for (s,wt) in zip(b.stages,(1/6,1/6,2/3))
        zz=vec(s.z)./abs(b.f);flat=reshape(s.parts,:,7)
        for n in eachindex(zz)
            bin=clamp(searchsortedlast(b.bin_edges,zz[n]),1,length(b.bin_edges)-1)
            for (r,mask) in enumerate(b.masks)
                mask[n]||continue
                b.conditional_count[bin,r]+=dt*wt
                for k in 1:7;b.conditional_sum[bin,k,r]+=dt*wt*flat[n,k];end
            end
        end
    end
    b.cumulative .+= parts;b.checks+=1
    b.max_relative=max(b.max_relative,err/max(scale,allowance))
    b.max_partition=max(b.max_partition,maximum(s.residual for s in b.stages))
    b.max_pressure=max(b.max_pressure,maximum(s.pg for s in b.stages))
    b.max_mass_div=max(b.max_mass_div,maximum(s.massdiv for s in b.stages))
    open(joinpath(b.output,"checks.csv"),"a") do io
        println(io,join((m.clock.iteration,b.t0,m.clock.time,dt,err,err/max(scale,allowance),maximum(s.residual for s in b.stages),maximum(s.pg for s in b.stages),maximum(s.massdiv for s in b.stages)),','))
    end
    b.iteration=-1;nothing
end
function (b::NativeVorticityBudget)(m)
    it=m.clock.iteration;stage=m.clock.stage
    if stage==1
        if b.iteration>=0 && it==b.iteration+1;nb_finish!(b,m)
        elseif b.iteration>=0 && it!=b.iteration;error("Skipped step");end
        b.iteration=it;b.t0=Float64(m.clock.time);b.z0=nb_pair(b,b.momentum)
        if b.origin===nothing;b.origin=copy(b.z0);b.origin_time=b.t0;b.cumulative=zeros(size(b.z0)...,7);end
        b.stages=Any[nothing,nothing,nothing];b.stage_times.=NaN
    end
    if it==b.iteration
        stage in 1:3||error("Unknown RK stage")
        b.stages[stage]=nb_terms(b);b.stage_times[stage]=Float64(m.clock.time)
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
    z=nb_pair(b,b.momentum);origin=b.origin===nothing ? z : b.origin
    parts=b.cumulative===nothing ? zeros(size(z)...,7) : b.cumulative
    residual=z.-origin.-dropdims(sum(parts;dims=3);dims=3)
    raw=[mean(z[mask].^3 .- origin[mask].^3) for mask in b.masks]
    central=[mean((z[mask].-mean(z[mask])).^3)-mean((origin[mask].-mean(origin[mask])).^3) for mask in b.masks]
    jldsave(joinpath(b.output,@sprintf("budget_%012d.jld2",round(Int,m.clock.time)));
        t=Float64(m.clock.time),origin_time=b.origin_time,Z=z,origin_Z=origin,cumulative_parts=parts,cumulative_residual=residual,
        third_moment_storage=raw,third_moment_parts=copy(b.third),third_moment_residual=raw.-vec(sum(b.third;dims=2)),
        central_third_storage=central,central_third_parts=copy(b.central_third),central_third_residual=central.-vec(sum(b.central_third;dims=2)),
        conditional_edges=b.bin_edges,conditional_sum=copy(b.conditional_sum),conditional_count=copy(b.conditional_count),
        regions=["full","interior","wall_strips"],term_names=SIGN_TERMS,checks=b.checks,f=b.f,
        max_relative_step_residual=b.max_relative,max_rhs_partition_residual=b.max_partition,max_pressure_curl=b.max_pressure,max_mass_divergence=b.max_mass_div,
        diagnostic_note="Centered stress and mean-flux decomposition is distinct from native WENO; discrepancy retained. Interior transport is not discarded as zero. Full-domain and region moments include mean corrections.")
    nothing
end
