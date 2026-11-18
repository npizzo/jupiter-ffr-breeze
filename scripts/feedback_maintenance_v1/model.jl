# Matched maintenance ablation; never modifies the frozen polynomial campaign.
include("../viscous2d_feedback512_v1/model.jl")
using Oceananigans.Advection: div_𝐯u, div_𝐯v
const MAINTENANCE_CASES=("full","no_low_quadratic")
mutable struct MaintenanceTerms{B}
    base::B
    quadratic::Matrix{ComplexF64}
    low::BitMatrix
    remove_low::Bool
end
function setup_maintenance(u,v;case,t=400.,dt=.0025)
    case in MAINTENANCE_CASES || error("Unsupported maintenance case")
    m,b=setup_polynomial(u,v;case="biased_plus",t,dt)
    m,MaintenanceTerms(b,zeros(ComplexF64,size(u)),BitMatrix((b.op.k2.>0).&(b.op.k2.<=4.)),case=="no_low_quadratic")
end
@noinline function quadratic_buffer!(r,q,bf)
    @inbounds for k in eachindex(q);z=real(q[k]);r[k]=complex(bf*z*z,0.);end
    nothing
end
@noinline function remove_low_quadratic!(reaction,quadratic,low)
    @inbounds for k in eachindex(reaction);low[k] && (reaction[k]-=quadratic[k]);end
    nothing
end
function (b::MaintenanceTerms)(m)
    b.remove_low || return b.base(m)
    p=b.base;o=p.op;nx,ny=o.nx,o.ny;idx=nx/o.lx;idy=ny/o.ly
    u=interior(m.velocities.u);v=interior(m.velocities.v)
    guard_cfl(u,v,nx,ny,idx,idy,p.dt)
    curl_to_buffer!(p.q,u,v,nx,ny,idx,idy)
    polynomial_buffer!(p.reaction,p.q,p.b_f,p.c,p.dt)
    quadratic_buffer!(b.quadratic,p.q,p.b_f)
    mul!(p.q,p.forward,p.q);mul!(p.reaction,p.forward,p.reaction);mul!(b.quadratic,p.forward,b.quadratic)
    remove_low_quadratic!(p.reaction,b.quadratic,b.low)
    combine_sources!(p.q,p.reaction,p.linear,p.invlap)
    mul!(p.q,p.backward,p.q)
    add_streamfunction_tendency!(interior(m.timestepper.Gⁿ.u),interior(m.timestepper.Gⁿ.v),p.q,nx,ny,idx,idy)
    nothing
end
maintenance_step!(sim,b;ceiling=.0025)=set_stable_step!(sim,b.base;ceiling)
function maintenance_hashes()
    paths=["model.jl","run.jl","check.jl","../viscous2d_feedback512_v1/model.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(x=>bytes2hex(sha256(read(joinpath(@__DIR__,x)))) for x in paths)
end
function check_source_checkpoint(d)
    @assert size(d["u"])==size(d["v"])==(512,512) && d["t"]==400.
    c=d["config"]
    @assert c["case"]=="biased_plus" && c["grid"]==[512,512] && c["forcing_band"]==[12.,16.]
    for (k,v) in [("a",.0002),("b",.002),("c",.0002),("f",1.),("viscosity",.0001),("dt_ceiling",.0025),("beta",0.)]
        @assert c[k]==v
    end
    @assert !c["allow_zonal_forcing"]
    nothing
end
function maintenance_sources(o,q;remove_low=false)
    qb=.002 .*q.^2;qb.-=mean(qb)
    if remove_low
        bh=FFTW.fft(qb);bh[(o.k2.>0).&(o.k2.<=4.)].=0.;qb=real.(FFTW.ifft(bh))
    end
    qc=-.0002 .*q.^3;qc.-=mean(qc)
    ql=band_curl(o,q,.0002)
    qv=spectral_curl(o,q,-.0001 .*o.lambda)
    (;quadratic=qb,cubic=qc,linear=ql,viscous=qv)
end
# Native advection evaluated at actual staggered faces, not u·grad(zeta).
# Called only at sparse saved times; hot timestep path is unchanged.
@noinline function native_advection_arrays!(gu,gv,grid,advection,velocities,nx,ny)
    @inbounds for j=1:ny,i=1:nx
        gu[i,j]=-div_𝐯u(i,j,1,grid,advection,velocities,velocities.u)
        gv[i,j]=-div_𝐯v(i,j,1,grid,advection,velocities,velocities.v)
    end
    nothing
end
function native_advection_curl(m,o)
    gu=zeros(o.nx,o.ny);gv=similar(gu)
    native_advection_arrays!(gu,gv,m.grid,m.advection.momentum,m.velocities,o.nx,o.ny)
    vorticity(o,gu,gv)
end
function maintenance_diagnostics(m,b)
    o=b.base.op;u,v=state(m);q=vorticity(o,u,v);qh=FFTW.fft(q);n=length(q)^2
    il=map(x->x==0 ? 0. : 1/x,o.lambda)
    ek=.5 .*abs2.(qh).*il./n;low=b.low
    s=maintenance_sources(o,q;remove_low=b.remove_low)
    adv=native_advection_curl(m,o)
    out=Dict{String,Float64}("energy"=>.5*(mean(abs2,u)+mean(abs2,v)),"energy_k_le_2"=>sum(ek[low]),
       "enstrophy"=>.5mean(abs2,q),"mean_zeta"=>mean(q),"skewness"=>mean(q.^3)/mean(abs2,q)^1.5,
       "enstrophy_k_ge_128"=>.5sum(abs2.(qh)[o.k2.>=128^2])/n)
    for (name,source) in pairs(merge((advection=adv,),s))
        work=real.(conj.(qh).*FFTW.fft(source)).*il./n
        out[string(name)*"_energy_rate"]=sum(work)
        out[string(name)*"_lowk_energy_rate"]=sum(work[low])
    end
    out["lowk_total_rate"]=sum(out[string(name)*"_lowk_energy_rate"] for name in (:advection,:linear,:quadratic,:cubic,:viscous))
    q,out
end
