# Native face-momentum tendencies divided by face density, not material acceleration.
@kernel function hotspot_terms!(out,coords,g,F,p,V,M,dyn,adv,closure,cf,clock,P,Pr,rr,oldrot)
    n=@index(Global,Linear);i=coords[1,n];j=coords[2,n];k=coords[3,n]
    rx=.5*(F.ρᵈ[i-1,j,k]+F.ρᵈ[i,j,k]);ry=.5*(F.ρᵈ[i,j-1,k]+F.ρᵈ[i,j,k]);rz=.5*(F.ρᵈ[i,j,k-1]+F.ρᵈ[i,j,k])
    out[1,n]=-((M.ρu[i+1,j,k]-M.ρu[i,j,k])/p.dx+(M.ρv[i,j+1,k]-M.ρv[i,j,k])/p.dy+(M.ρw[i,j,k+1]-M.ρw[i,j,k])/p.dz)/F.ρᵈ[i,j,k]
    out[2,n]=-(P[i,j,k]-P[i-1,j,k])/p.dx/rx
    out[3,n]=rotation_x(i,j,k,F,p)/rx
    out[4,n]=-AM.x_momentum_flux_divergence(i,j,k,g,adv,M,V,dyn)/rx
    out[5,n]=(-AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,dyn.dry_density,closure,cf,clock,F,nothing)+∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p))/rx
    out[6,n]=-OC.x_f_cross_U(i,j,k,g,oldrot,M)/rx
    if j>1
        out[7,n]=-(P[i,j,k]-P[i,j-1,k])/p.dy/ry
        out[8,n]=rotation_y(i,j,k,F,p)/ry
        out[9,n]=-AM.y_momentum_flux_divergence(i,j,k,g,adv,M,V,dyn)/ry
        out[10,n]=(-AM.∂ⱼ_𝒯₂ⱼ(i,j,k,g,dyn.dry_density,closure,cf,clock,F,nothing)+∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p))/ry
        out[11,n]=-OC.y_f_cross_U(i,j,k,g,oldrot,M)/ry
    else
        for q=7:11;out[q,n]=0.;end
    end
    if k>1
        out[12,n]=(-((P[i,j,k]-Pr[i,j,k])-(P[i,j,k-1]-Pr[i,j,k-1]))/p.dz-GRAV*(rz-.5*(rr[i,j,k]+rr[i,j,k-1])))/rz
        out[13,n]=rotation_z(i,j,k,F,p)/rz
        out[14,n]=-AM.z_momentum_flux_divergence(i,j,k,g,adv,M,V,dyn)/rz
        out[15,n]=(-AM.∂ⱼ_𝒯₃ⱼ(i,j,k,g,dyn.dry_density,closure,cf,clock,F,nothing)+∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p))/rz
        out[16,n]=-OC.z_f_cross_U(i,j,k,g,oldrot,M)/rz
    else
        for q=12:16;out[q,n]=0.;end
    end
    out[17,n]=p.nu[i,j,k]
    out[18,n]=external_heating(i,j,k,F,p)/F.ρᵈ[i,j,k]
end
function make_hotspot_probe(m,p,out)
    on(a)=Oceananigans.Architectures.on_architecture(m.architecture,a)
    # Fixed top-wall probes, nearby interior probes, and a moving maximum.
    fixed=[(1,1,p.nz),(1,p.ny,p.nz),(1,2,p.nz-2),(1,p.ny-1,p.nz-2),(1,div(p.ny,2),p.nz),(1,p.ny,round(Int,.7/p.dz))]
    coords=on(zeros(Int,3,length(fixed)+1));buf=on(zeros(18,length(fixed)+1));backend=get_backend(buf)
    pp=NamedTuple{Tuple(k for k in keys(p) if k!=:seed_info)}(Tuple(p[k] for k in keys(p) if k!=:seed_info))
    oldrot=ChannelCoriolis(BetaPlane(;f₀=p.f-p.ly*p.beta/2,β=p.beta,scheme=Oceananigans.Coriolis.EnergyConserving()),ConstantCartesianCoriolis(;fx=0.,fy=p.fh,fz=0.))
    path=joinpath(out,"hotspots.csv")
    open(path,"w") do io
        println(io,"time,probe,i,j,k,x,y,z,u,v,w,rho,T,mach,along_axis,mass_convergence_over_rho,px,cx,ax,vx,cx_old,py,cy,ay,vy,cy_old,pz_gravity,cz,az,vz,cz_old,nu,external_heat_over_rho")
    end
    summary=joinpath(out,"hotspot_extent.csv")
    open(summary,"w") do io;println(io,"time,mach_max,outside_mach,bulk_mach,fraction_above04,north_ke,south_ke");end
    z=reshape([zcell(k,p) for k=1:p.nz],1,1,p.nz);y=reshape([(j-.5)*p.dy-p.ly/2 for j=1:p.ny],1,p.ny,1)
    corner=(z.>.9).&(abs.(y).>3.5);bulk=(z.<.78).&(abs.(y).<3.5)
    function record(model,a,t)
        mach=sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T)
        points=[fixed;[Tuple(argmax(mach))]];xyz=hcat([collect(q) for q in points]...);copyto!(coords,xyz)
        hotspot_terms!(backend,32)(buf,coords,model.grid,Oceananigans.fields(model),pp,model.velocities,model.momentum,model.dynamics,model.advection.momentum,model.closure,model.closure_fields,model.clock,model.dynamics.pressure,model.dynamics.reference_state.pressure,model.dynamics.reference_state.density,oldrot;ndrange=length(points))
        KernelAbstractions.synchronize(backend);b=Array(buf)
        open(path,"a") do io
            for (n,(i,j,k)) in enumerate(points)
                f=p.f+p.beta*((j-.5)*p.dy-p.ly/2);along=(p.fh*a.v[i,j,k]+f*a.w[i,j,k])/hypot(p.fh,f)
                vals=(t,n,i,j,k,(i-.5)*p.dx,(j-.5)*p.dy-p.ly/2,zcell(k,p),a.u[i,j,k],a.v[i,j,k],a.w[i,j,k],a.rho[i,j,k],a.T[i,j,k],mach[i,j,k],along,b[:,n]...)
                @assert all(isfinite,vals);println(io,join(vals,','))
            end
        end
        ke=.5 .*a.rho.*(a.u.^2 .+a.v.^2 .+a.w.^2);dv=p.dx*p.dy*p.dz
        open(summary,"a") do io;println(io,join((t,maximum(mach),maximum(mach.*(.!corner)),maximum(mach.*bulk),mean(mach.>.4),dv*sum(ke.*corner.*(y.>0)),dv*sum(ke.*corner.*(y.<0))),','));end
    end
end
