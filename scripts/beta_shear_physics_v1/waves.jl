function sponge_power(a,p,cfg)
    get(cfg,"sponge",false) || return 0.
    rate=get(cfg,"sponge_rate",8.);depth=get(cfg,"sponge_depth",.18)
    rz=0.5 .* (a.rho[:,:,1:end-1].+a.rho[:,:,2:end]);z=[k*p.dz for k=1:p.nz-1]
    q=clamp.((z.-(p.lz-depth))./depth,0.,1.);r=rate.*q.^2 .*(3 .-2 .*q)
    power=-p.dx*p.dy*p.dz*sum(rz.*a.wf[:,:,2:end-1].^2 .*reshape(r,1,1,:))
    if get(cfg,"horizontal_sponge",false)
        rx=0.5 .* (a.rho.+circshift(a.rho,(1,0,0)));ry=yface_density(a.rho)
        rates=reshape([upper_rate(k,p) for k=1:p.nz],1,1,p.nz)
        target=reshape([jet_target((j-.5)*p.dy-p.ly/2,zcell(k,p),p) for j=1:p.ny,k=1:p.nz],1,p.ny,p.nz)
        power-=p.dx*p.dy*p.dz*(sum(rates.*rx.*a.uf.*(a.uf.-target))+sum(rates.*ry.*a.vf.^2))
    end
    power
end
function wave_output(m,p,cfg,out,t,dt,a)
    get(cfg,"wave_diagnostics",true) || return
    path=joinpath(out,"wave.csv")
    if !isfile(path)
        open(path,"w") do io;println(io,"time,dt,mach,wmax,deep_mach,cap_w_rms,top_w_rms,vertical_grid_tail,hotspot_w,hotspot_T,sponge_power,substeps_stage1,substeps_stage2,substeps_stage3");end
    end
    z=[zcell(k,p) for k=1:p.nz];mach=sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T)
    i=clamp(round(Int,2.1875/p.dx+.5),1,p.nx);j=clamp(round(Int,(-2.1875+p.ly/2)/p.dy+.5),1,p.ny)
    k=argmin(abs.(z.-1.33125));deep=findall(z.<DC);cap=findall(z.>.82);top=findall(z.>1.42)
    q=dct(a.w,3);spec=vec(sum(abs2.(q);dims=(1,2)));tail=sum(spec[ceil(Int,2p.nz/3):end])/max(sum(spec),1e-99)
    # Record actual counts with the same installed function and substepper constants.
    st=m.timestepper.substepper
    counts=[Breeze.CompressibleEquations.stage_substep_count_and_size(st.substep_distribution,b,dt,st.minimum_horizontal_spacing,m.thermodynamic_constants,st.acoustic_cfl,st.substeps)[1] for b in (1/3,1/2,1.)]
    open(path,"a") do io;println(io,join((t,dt,maximum(mach),maximum(abs,a.w),maximum(mach[:,:,deep]),sqrt(mean(a.w[:,:,cap].^2)),sqrt(mean(a.w[:,:,top].^2)),tail,a.w[i,j,k],a.T[i,j,k],sponge_power(a,p,cfg),counts...),','));end
    # Small columns every sample. Whole-model fields remain on Unity.
    P=Array(interior(m.dynamics.pressure))
    jldsave(joinpath(out,@sprintf("wave_%07d.jld2",round(Int,10000t)));t,z,x=(i-.5)*p.dx,y=(j-.5)*p.dy-p.ly/2,w=a.w[i,j,:],T=a.T[i,j,:],rho=a.rho[i,j,:],P=P[i,j,:],
        mean_T=vec(mean(a.T;dims=(1,2))),mean_P=vec(mean(P;dims=(1,2))),mean_rho=vec(mean(a.rho;dims=(1,2))),
        wave_pressure_flux=vec(mean((P.-mean(P;dims=(1,2))).*a.w;dims=(1,2))),w_profile_rms=vec(sqrt.(mean(a.w.^2;dims=(1,2)))),vertical_spectrum=spec)
end
function initial_audit(cfg,out)
    mkpath(out);m,p=build_model(cfg;arch=CPU());a=arrays(m);F=Oceananigans.fields(m);g=m.grid;M=m.momentum;V=m.velocities
    P=Array(interior(m.dynamics.pressure));Pr=Array(interior(m.dynamics.reference_state.pressure));rr=Array(interior(m.dynamics.reference_state.density))
    mc=Array(interior(M.ρu));mn=Array(interior(M.ρv));mw=Array(interior(M.ρw))
    dr=-(circshift(mc,(-1,0,0)).-mc)./p.dx.-diff(mn;dims=2)./p.dy.-diff(mw;dims=3)./p.dz
    ax=zeros(p.nx,p.ny,p.nz);ay=zeros(p.nx,p.ny+1,p.nz);az=zeros(p.nx,p.ny,p.nz+1);dq=zero(ax)
    interior(p.umean).=mean(a.uf;dims=1)
    for k=1:p.nz,j=1:p.ny,i=1:p.nx
        im=i==1 ? p.nx : i-1
        ax[i,j,k]=(-(P[i,j,k]-P[im,j,k])/p.dx+fx(i,j,k,g,m.clock,F,p)-AM.x_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)-AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing))/(0.5*(a.rho[i,j,k]+a.rho[im,j,k]))
        if j>1
            ay[i,j,k]=(-(P[i,j,k]-P[i,j-1,k])/p.dy+fy(i,j,k,g,m.clock,F,p)-AM.y_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)-AM.∂ⱼ_𝒯₂ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing))/(0.5*(a.rho[i,j,k]+a.rho[i,j-1,k]))
        end
        if k>1
            rz=0.5*(a.rho[i,j,k]+a.rho[i,j,k-1]);rref=0.5*(rr[1,1,k]+rr[1,1,k-1])
            az[i,j,k]=(-((P[i,j,k]-Pr[1,1,k])-(P[i,j,k-1]-Pr[1,1,k-1]))/p.dz-GRAV*(rz-rref)+fz(i,j,k,g,m.clock,F,p)-AM.z_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)-AM.∂ⱼ_𝒯₃ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing))/rz
        end
        theta=m.formulation.potential_temperature[i,j,k]
        dq[i,j,k]=-AM.div_ρUc(i,j,k,g,m.advection.ρθ,m.dynamics.dry_density,V,m.formulation.potential_temperature)+theta/(2.5*a.T[i,j,k])*heat_source(i,j,k,g,m.clock,F,p)
    end
    rho_departure=a.rho./reshape(vec(rr),1,1,p.nz).-1;zeta=channel_vorticity(a.uf,a.vf,p)
    jldsave(joinpath(out,"initial_audit.jld2");ax,ay,az,dr,dq,rho_departure,zeta,rho=a.rho,T=a.T,P,Pr,rr,dx=p.dx,dy=p.dy,dz=p.dz)
    write_toml(joinpath(out,"initial_audit.toml"),Dict("native_energy"=>audit(m,p),"max_vertical_momentum_tendency_over_density"=>maximum(abs,az),"max_horizontal_momentum_tendency_over_density"=>max(maximum(abs,ax),maximum(abs,ay)),"max_relative_mass_tendency"=>maximum(abs,dr./a.rho),"max_density_departure"=>maximum(abs,rho_departure),"note"=>"ax/ay/az are momentum tendencies divided by face density, not full material acceleration; include native advection, pressure, rotation, viscosity and imposed forcing."))
end
