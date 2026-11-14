# Observer only. All arrays below remain on Unity. Native finite-volume fluxes
# are decomposed on their own faces, then differenced onto momentum/scalar cells.
import Oceananigans.Advection as OA
import Oceananigans.Operators as OP
import Breeze.Advection as BA

@inline my(i,j,k,g,a,M,u) = OA._advective_momentum_flux_Vu(i,j,k,g,a,M[2],u)
@inline mz(i,j,k,g,a,M,u) = OA._advective_momentum_flux_Wu(i,j,k,g,a,M[3],u)
@inline vy(i,j,k,g,M) = OP.ℑxᶠᵃᵃ(i,j,k,g,OP.Ay_qᶜᶠᶜ,M[2])
@inline wz(i,j,k,g,M) = OP.ℑxᶠᵃᵃ(i,j,k,g,OP.Az_qᶜᶜᶠ,M[3])
@inline uy(i,j,k,g,u) = OP.ℑyᵃᶠᵃ(i,j,k,g,u)
@inline uz(i,j,k,g,u) = OP.ℑzᵃᵃᶠ(i,j,k,g,u)
@inline puy(i,j,k,g,M,u) = vy(i,j,k,g,M)*uy(i,j,k,g,u)
@inline puz(i,j,k,g,M,u) = wz(i,j,k,g,M)*uz(i,j,k,g,u)
@inline hy(i,j,k,g,r,a,U,t) = BA.tracer_mass_flux_y(i,j,k,g,r,a,U.v,t)
@inline hz(i,j,k,g,r,a,U,t) = BA.tracer_mass_flux_z(i,j,k,g,r,a,U.w,t)
@inline hmy(i,j,k,g,r,U) = OP.Ayᶜᶠᶜ(i,j,k,g)*OP.ℑyᵃᶠᵃ(i,j,k,g,r)*U.v[i,j,k]
@inline hmz(i,j,k,g,r,U) = OP.Azᶜᶜᶠ(i,j,k,g)*OP.ℑzᵃᵃᶠ(i,j,k,g,r)*U.w[i,j,k]
@inline phy(i,j,k,g,r,U,t) = hmy(i,j,k,g,r,U)*uy(i,j,k,g,t)
@inline phz(i,j,k,g,r,U,t) = hmz(i,j,k,g,r,U)*uz(i,j,k,g,t)
@inline hadv(i,j,k,g,r,a,U,t) = -AM.div_ρUc(i,j,k,g,a,r,U,t)
@inline hdiff(i,j,k,g,r,cl,cf,t,clock,fields) = -AM.∇_dot_Jᶜ(i,j,k,g,r,cl,cf,Val(1),t,clock,fields,nothing)
@inline hforce(i,j,k,g,F,clock,fields) = F(i,j,k,g,clock,fields)
@inline hbalance(i,j,k,g,r,t,rate) = -r[i,j,k]*t[i,j,k]*rate[1]

mutable struct FaceBudget
    flux_fields::Any
    heat_fields::Any
    heat_total::Any
    heat_state::Any
    volume::Float64
    bottom::Matrix{Float64}
    origin::Any
    integrated_momentum::Any
    integrated_heat::Any
    integrated_heat_transport::Any
    origin_time::Float64
    end_time::Float64
    max_momentum_flux_error::Float64
    max_heat_partition_error::Float64
    max_heat_increment_error::Float64
    steps::Int
end
const FACE_MOMENTUM_NAMES=["mean_y","mean_z","eddy_y","eddy_z","reconstruction_y","reconstruction_z"]
const FACE_HEAT_NAMES=["advection","diffusion","stable_cap","thermal_feedback","bottom_flux"]
# 'reconstruction' means native WENO minus a specified linear face product,
# not an independently measured dissipative stress.
function FaceBudget(m,cfg,rate)
    r=AM.continuity_density(m.dynamics); M=AM.advecting_momentum(m)
    u=m.velocities.u;t=m.formulation.potential_temperature
    U=AM.thermodynamic_transport_velocities(m); a=m.advection.momentum; ah=m.advection.ρθ
    fields=Oceananigans.fields(m)
    reduce_at(loc,func,args...)=zonal_reduction(KernelFunctionOperation{loc[1],loc[2],loc[3]}(func,m.grid,args...))
    ym=(Face,Face,Center);zm=(Face,Center,Face);yh=(Center,Face,Center);zh=(Center,Center,Face);cc=(Center,Center,Center)
    ff=(reduce_at(ym,my,a,M,u),reduce_at(zm,mz,a,M,u),
        reduce_at(ym,vy,M),reduce_at(zm,wz,M),reduce_at(ym,uy,u),reduce_at(zm,uz,u),
        reduce_at(ym,puy,M,u),reduce_at(zm,puz,M,u),
        reduce_at(yh,hy,r,ah,U,t),reduce_at(zh,hz,r,ah,U,t),
        reduce_at(yh,hmy,r,U),reduce_at(zh,hmz,r,U),reduce_at(yh,uy,t),reduce_at(zh,uz,t),
        reduce_at(yh,phy,r,U,t),reduce_at(zh,phz,r,U,t))
    hf=(reduce_at(cc,hadv,r,ah,U,t),reduce_at(cc,hdiff,r,m.closure,m.closure_fields,t,m.clock,fields),
        reduce_at(cc,hforce,m.forcing.ρθ,m.clock,fields),reduce_at(cc,hbalance,r,t,rate))
    ny,nz=cfg["ny"],cfg["nz"];dz=cfg["height"]/nz
    R=8.314462618/cfg["molar_mass"];exner=(cfg["base_pressure"]/1e5)^(R/cfg["heat_capacity"])
    bottom=zeros(ny,nz);bottom[:,1].=cfg["heat_flux"]/(cfg["heat_capacity"]*exner*dz)
    FaceBudget(ff,hf,zonal_reduction(m.timestepper.Gⁿ.ρθ),zonal_reduction(m.formulation.potential_temperature_density),
        (cfg["lx"]/cfg["nx"])*(cfg["ly"]/ny)*dz,bottom,nothing,zeros(ny,nz,6),zeros(ny,nz,5),zeros(ny,nz,6),NaN,NaN,0,0,0,0)
end
function flux_partition(f,volume)
    # Each f contains native_y,z; mass_y,z; velocity/scalar_y,z; linear_product_y,z.
    mean_y=f[3].*f[5];mean_z=f[4].*f[6]
    eddy_y=f[7].-mean_y;eddy_z=f[8].-mean_z
    remainder_y=f[1].-f[7];remainder_z=f[2].-f[8]
    dy(q)=-diff(q;dims=1)./volume
    dz(q)=-diff(q;dims=2)./volume
    return cat(dy(mean_y),dz(mean_z),dy(eddy_y),dz(eddy_z),dy(remainder_y),dz(remainder_z);dims=3)
end
function face_terms(b)
    fs=zonal_host.(b.flux_fields);hs=zonal_host.(b.heat_fields)
    momentum=flux_partition(fs[1:8],b.volume)
    heat_transport=flux_partition(fs[9:16],b.volume)
    heat=cat(hs[1],hs[2],hs[3].-hs[4],hs[4],b.bottom;dims=3)
    # TendencyCallsite is before compute_flux_bc_tendencies!: add the known,
    # uniform bottom rho-theta flux explicitly, and validate against increments.
    htotal=zonal_host(b.heat_total).+b.bottom
    herr=maximum(abs,dropdims(sum(heat;dims=3);dims=3).-htotal)
    scale=max(maximum(abs,heat),maximum(abs,htotal),1e-15)
    herr<=2e-10*scale+1e-13 || error("Heat RHS partition failed: $herr, $scale")
    transport_error=maximum(abs,dropdims(sum(heat_transport;dims=3);dims=3).-hs[1])
    transport_error<=2e-10*max(maximum(abs,hs[1]),1e-15)+1e-12 || error("Native heat face flux mismatch: $transport_error")
    state=zonal_host(b.heat_state)
    return (;momentum,heat,heat_transport,state,herr)
end
function integrate_face!(b,stages,dt,m)
    mm=map(s->s.face,stages)
    quad(k)=dt.*(getproperty(mm[1],k)./6 .+ getproperty(mm[2],k)./6 .+ 2getproperty(mm[3],k)./3)
    mom=quad(:momentum);heat=quad(:heat);ht=quad(:heat_transport)
    for s in stages
        err=maximum(abs,dropdims(sum(s.face.momentum;dims=3);dims=3).-s.parts[:,:,1])
        scale=max(maximum(abs,s.parts[:,:,1]),1e-15)
        err<=2e-10*scale+1e-13 || error("Native momentum face flux mismatch: $err, $scale")
        b.max_momentum_flux_error=max(b.max_momentum_flux_error,err)
        b.max_heat_partition_error=max(b.max_heat_partition_error,s.face.herr)
    end
    current=zonal_host(b.heat_state);initial=mm[1].state
    predicted=dropdims(sum(heat;dims=3);dims=3)
    heat_error=maximum(abs,current.-initial.-predicted)
    allowance=8192eps(Float64)*max(maximum(abs,current),1.0)+1e-8*maximum(abs,heat)
    heat_error<=allowance || error("Native heat increment mismatch: $heat_error > $allowance")
    if b.origin===nothing
        b.origin=copy(initial);b.origin_time=Float64(m.clock.time)-dt
    end
    b.end_time=Float64(m.clock.time);b.steps+=1
    b.integrated_momentum .+= mom;b.integrated_heat .+= heat;b.integrated_heat_transport .+= ht
    b.max_heat_increment_error=max(b.max_heat_increment_error,heat_error)
    nothing
end
function face_snapshot(b,m)
    current=zonal_host(b.heat_state)
    origin=b.origin===nothing ? current : b.origin
    residual=current.-origin.-dropdims(sum(b.integrated_heat;dims=3);dims=3)
    return (;face_origin_time=b.origin_time,face_end_time=b.end_time,face_steps=b.steps,
        face_integrated_momentum=copy(b.integrated_momentum),face_momentum_names=FACE_MOMENTUM_NAMES,
        face_integrated_heat=copy(b.integrated_heat),face_heat_names=FACE_HEAT_NAMES,
        face_integrated_heat_transport=copy(b.integrated_heat_transport),face_heat_origin=origin,face_heat_current=current,
        face_heat_residual=residual,face_max_momentum_flux_error=b.max_momentum_flux_error,
        face_max_heat_partition_error=b.max_heat_partition_error,face_max_heat_increment_error=b.max_heat_increment_error)
end
