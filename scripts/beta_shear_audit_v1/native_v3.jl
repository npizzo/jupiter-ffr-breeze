# Read-only native-operator audit; all execution and source fields stay on Unity.
include("../beta_shear_diagnosis_v2/model.jl")
include("../beta_shear_diagnosis_v2/initialization.jl")
const AM=Breeze.AtmosphereModels
const OC=Oceananigans.Coriolis
function audit(m,p)
    F=Oceananigans.fields(m);g=m.grid;M=m.momentum;V=m.velocities
    a=arrays(m);rho=a.rho;P=Array(interior(m.dynamics.pressure));theta=Array(interior(m.formulation.potential_temperature))
    Pr=Array(interior(m.dynamics.reference_state.pressure));rr=Array(interior(m.dynamics.reference_state.density))
    interior(p.umean).=mean(a.uf;dims=1)
    nx,ny,nz=size(rho);dv=p.dx*p.dy*p.dz
    # ρu,ρv,ρw are native face momenta. Every bounded normal wall velocity is zero.
    mc=Array(interior(M.ρu));mn=Array(interior(M.ρv));mw=Array(interior(M.ρw))
    dr=-(circshift(mc,(-1,0,0)).-mc)./p.dx.-diff(mn;dims=2)./p.dy.-diff(mw;dims=3)./p.dz
    cov=0.25 .*(a.uf.^2 .+circshift(a.uf,(-1,0,0)).^2 .+a.vf[:,1:end-1,:].^2 .+a.vf[:,2:end,:].^2 .+a.wf[:,:,1:end-1].^2 .+a.wf[:,:,2:end].^2)
    r=Dict{String,Float64}("mass_tendency"=>sum(dr)*dv,"advective_mass_kinetic"=>-sum(cov.*dr)*dv,"potential_energy_advection"=>sum(reshape([GRAV*zcell(k,p) for k=1:nz],1,1,nz).*dr)*dv)
    keyslist=["momentum_advection_work","thermodynamic_advection","coriolis_work","horizontal_coriolis_work","vertical_coriolis_work","pressure_gravity_work","native_viscosity_work","trace_viscosity_work","viscous_heat","entropy_heat","jet_work","external_heat","thermal_source_identity_error","top_entropy_heat","top_viscous_heat"]
    for key in keyslist;r[key]=0.;end
    for k=1:nz,j=1:ny,i=1:nx
        adv=-AM.div_ρUc(i,j,k,g,m.advection.ρθ,m.dynamics.dry_density,V,m.formulation.potential_temperature)
        r["thermodynamic_advection"]+=2.5*a.T[i,j,k]/theta[i,j,k]*adv
        div=divergence(i,j,k,F,p)
        hv=2p.nu[i,j,k]*rho[i,j,k]*(ΣᵢⱼΣᵢⱼᶜᶜᶜ(i,j,k,g,V.u,V.v,V.w)-div^2/3)
        he=p.entropy_diffusion_factor*sgs_heating(i,j,k,F,p)
        hx=p.heat*(conduction(i,j,k,F.T,p)-rho[i,j,k]*(a.T[i,j,k]-p.tt)*cool_rate(k,p))
        r["viscous_heat"]+=hv;r["entropy_heat"]+=he;r["external_heat"]+=hx
        r["thermal_source_identity_error"]+=heat_source(i,j,k,g,m.clock,F,p)-(hv+he+hx)
        if k==nz;r["top_entropy_heat"]+=he;r["top_viscous_heat"]+=hv;end
        u=a.uf[i,j,k]
        r["momentum_advection_work"]-=u*AM.x_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)
        r["coriolis_work"]-=u*OC.x_f_cross_U(i,j,k,g,m.coriolis,M)
        r["horizontal_coriolis_work"]-=u*OC.x_f_cross_U(i,j,k,g,m.coriolis.horizontal,momentum_aliases(M))
        r["vertical_coriolis_work"]-=u*OC.x_f_cross_U(i,j,k,g,m.coriolis.vertical,momentum_aliases(M))
        r["native_viscosity_work"]-=u*AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing)
        r["trace_viscosity_work"]+=u*∂xᶠᶜᶜ(i,j,k,g,trace_stress,F,p)
        r["jet_work"]+=u*jet_momentum(i,j,k,F,p)
        im=i==1 ? nx : i-1
        r["pressure_gravity_work"]-=u*(P[i,j,k]-P[im,j,k])/p.dx
        if j>1
            v=a.vf[i,j,k]
            r["momentum_advection_work"]-=v*AM.y_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)
            r["coriolis_work"]-=v*OC.y_f_cross_U(i,j,k,g,m.coriolis,M)
            r["horizontal_coriolis_work"]-=v*OC.y_f_cross_U(i,j,k,g,m.coriolis.horizontal,momentum_aliases(M))
            r["vertical_coriolis_work"]-=v*OC.y_f_cross_U(i,j,k,g,m.coriolis.vertical,momentum_aliases(M))
            r["native_viscosity_work"]-=v*AM.∂ⱼ_𝒯₂ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing)
            r["trace_viscosity_work"]+=v*∂yᶜᶠᶜ(i,j,k,g,trace_stress,F,p)
            r["pressure_gravity_work"]-=v*(P[i,j,k]-P[i,j-1,k])/p.dy
        end
        if k>1
            w=a.wf[i,j,k]
            r["momentum_advection_work"]-=w*AM.z_momentum_flux_divergence(i,j,k,g,m.advection.momentum,M,V,m.dynamics)
            r["coriolis_work"]-=w*OC.z_f_cross_U(i,j,k,g,m.coriolis,M)
            r["horizontal_coriolis_work"]-=w*OC.z_f_cross_U(i,j,k,g,m.coriolis.horizontal,momentum_aliases(M))
            r["vertical_coriolis_work"]-=w*OC.z_f_cross_U(i,j,k,g,m.coriolis.vertical,momentum_aliases(M))
            r["native_viscosity_work"]-=w*AM.∂ⱼ_𝒯₃ⱼ(i,j,k,g,m.dynamics.dry_density,m.closure,m.closure_fields,m.clock,F,nothing)
            r["trace_viscosity_work"]+=w*∂zᶜᶜᶠ(i,j,k,g,trace_stress,F,p)
            dp=((P[i,j,k]-Pr[1,1,k])-(P[i,j,k-1]-Pr[1,1,k-1]))/p.dz
            buoy=-GRAV*.5*(rho[i,j,k]+rho[i,j,k-1]-rr[1,1,k]-rr[1,1,k-1])
            r["pressure_gravity_work"]+=w*(-dp+buoy)
        end
    end
    for key in keyslist;r[key]*=dv;end
    r["viscous_energy_defect"]=r["native_viscosity_work"]+r["trace_viscosity_work"]+r["viscous_heat"]
    r["inviscid_non_coriolis_defect"]=r["momentum_advection_work"]+r["advective_mass_kinetic"]+r["thermodynamic_advection"]+r["potential_energy_advection"]+r["pressure_gravity_work"]
    r["semidiscrete_energy_defect"]=r["viscous_energy_defect"]+r["inviscid_non_coriolis_defect"]+r["coriolis_work"]+r["entropy_heat"]
    r["physical_source_power"]=r["jet_work"]+r["external_heat"]
    r["predicted_energy_derivative"]=r["semidiscrete_energy_defect"]+r["physical_source_power"]
    r["Tmin"]=minimum(a.T)
    @assert all(isfinite,values(r))
    return r
end
function main()
    out=ARGS[1];mkpath(out)
    cases=[("baseline","data/beta-shear-diag-v2-2136871/baseline",[0,10,20,30]),("no_entropy_diffusion","data/beta-shear-diag-v2-2136871/no_entropy_diffusion",[20,30]),("no_horizontal_rotation","data/beta-shear-diag-v2-2136871/no_horizontal_rotation",[20,30]),("original_3d","data/beta-shear-v1-pilot-2134379/parent",[20,100,140])]
    for (name,root,times) in cases
        cfg=TOML.parsefile(joinpath(root,"metadata.toml"));m,p=build_model(cfg;arch=CPU(),perturb=false)
        for t in times
            a=load(joinpath(root,"fields3d",@sprintf("state_%07d.jld2",t)))
            rho=Float64.(a["density"]);T=Float64.(a["temperature"]);theta=(rho.*T).^.6 ./rho
            set!(m;ρ=rho,θ=theta,u=Float64.(a["u_face"]),v=Float64.(a["v_face"]),w=Float64.(a["w_face"]),qᵗ=0.)
            Oceananigans.TimeSteppers.update_state!(m)
            r=audit(m,p);r["time"]=t
            write_toml(joinpath(out,name*"_"*string(t)*".toml"),r)
            println("AUDIT ",name," ",t," ",r);flush(stdout)
        end
        m=nothing;GC.gc()
    end
end
main()
