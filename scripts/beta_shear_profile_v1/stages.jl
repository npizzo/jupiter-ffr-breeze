using Oceananigans.Utils: KernelAbstractions
using .KernelAbstractions: @kernel, @index, get_backend
@kernel function point_tendencies!(out,coords,g,F,p,V,adv,rho,theta,q,G,U0)
    n=@index(Global,Linear);i=coords[1,n];j=coords[2,n];k=coords[3,n]
    div=divergence(i,j,k,F,p)
    hv=2p.nu[i,j,k]*rho[i,j,k]*(ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,F.u,F.v,F.w)-div^2/3)
    hs=p.entropy_diffusion_factor*sgs_heating(i,j,k,F,p)
    he=external_heating(i,j,k,F,p)
    fac=theta[i,j,k]/(2.5F.T[i,j,k])
    out[1,n]=q[i,j,k];out[2,n]=G[i,j,k];out[3,n]=U0[i,j,k]
    out[4,n]=-AM.div_ρUc(i,j,k,g,adv,rho,V,theta)
    out[5,n]=fac*hv;out[6,n]=fac*hs;out[7,n]=fac*he
    out[8,n]=rho[i,j,k];out[9,n]=F.T[i,j,k]
end
@kernel function point_states!(out,coords,q,G,U0,rho,T)
    n=@index(Global,Linear);i=coords[1,n];j=coords[2,n];k=coords[3,n]
    out[1,n]=q[i,j,k];out[2,n]=G[i,j,k];out[3,n]=U0[i,j,k]
    out[4,n]=rho[i,j,k];out[5,n]=T[i,j,k]
end
function stage_recorder(sim,p,out)
    m=sim.model;nx,ny,nz=size(m.grid)
    # Fixed patch around the previously observed first failure; small-model preflight uses its corner.
    xyz=hcat([[i,j,k] for k=max(1,nz-1):nz,j=max(1,min(200,ny)):min(204,ny),i=max(1,min(83,nx)):min(87,nx)]...)
    arch=m.architecture;on(a)=Oceananigans.Architectures.on_architecture(arch,a)
    coords=on(xyz);buf=on(zeros(9,size(xyz,2)));backend=get_backend(buf)
    # Forcing parameters are GPU-safe before the host-only initialization report is added.
    pp=NamedTuple{Tuple(k for k in keys(p) if k!=:seed_info)}(Tuple(p[k] for k in keys(p) if k!=:seed_info))
    io=open(joinpath(out,"stage_terms.csv"),"w")
    println(io,"event,time,iteration,next_stage,last_stage_dt,last_dt,i,j,k,q,G,U0,advection,viscous_source,entropy_source,external_source,rho,T")
    function write_rows(kind,model,a)
        for n=1:size(xyz,2)
            vals=kind=="tendency" ? a[:,n] : [a[1,n],a[2,n],a[3,n],NaN,NaN,NaN,NaN,a[4,n],a[5,n]]
            println(io,join((kind,model.clock.time,model.clock.iteration,model.clock.stage,model.clock.last_stage_Δt,model.clock.last_Δt,xyz[:,n]...,vals...),','))
        end
        flush(io)
    end
    function tend(model)
        point_tendencies!(backend,64)(buf,coords,model.grid,Oceananigans.fields(model),pp,AM.thermodynamic_transport_velocities(model),model.advection.ρθ,model.dynamics.dry_density,model.formulation.potential_temperature,model.formulation.potential_temperature_density,model.timestepper.Gⁿ.ρθ,model.timestepper.U⁰.ρθ;ndrange=size(xyz,2))
        KernelAbstractions.synchronize(backend);write_rows("tendency",model,Array(buf))
    end
    function state(model)
        point_states!(backend,64)(buf,coords,model.formulation.potential_temperature_density,model.timestepper.Gⁿ.ρθ,model.timestepper.U⁰.ρθ,model.dynamics.dry_density,model.temperature;ndrange=size(xyz,2))
        KernelAbstractions.synchronize(backend);write_rows("state",model,Array(buf))
    end
    add_callback!(sim,state,IterationInterval(1);name=:stage_state,callsite=Oceananigans.UpdateStateCallsite())
    add_callback!(sim,tend,IterationInterval(1);name=:stage_tendency,callsite=Oceananigans.TendencyCallsite())
    io
end
