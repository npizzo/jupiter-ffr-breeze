include("../instability2d_v1/core.jl")
using .Instability2D, Oceananigans, Statistics, LinearAlgebra, Random, TOML, Printf
using Oceananigans: TendencyCallsite
using Oceananigans.Solvers: FFTW

const TERMS = ("advection_coriolis", "instability", "viscosity", "cubic_drag")
mutable struct Audit
    op::BandOperator
    amplitude::Float64
    nu::Float64
    alpha::Float64
    remove_mean_drag::Bool
    dt::Float64
    stages::Vector{Any}
    iteration::Int
    finished::Int
    third::Vector{Float64}
    absolute_third::Vector{Float64}
    max_increment_error::Float64
    max_moment_error::Float64
    max_divergence::Float64
    max_cfl::Float64
    max_drag_work::Float64
    min_injection_work::Float64
    max_viscous_work::Float64
    energy::Vector{Float64}
    absolute_energy::Vector{Float64}
    max_energy_error::Float64
end
Audit(o;a=0.009375,nu=0.001,alpha=1.0,dt=0.01,remove_mean_drag=true)=
    Audit(o,a,nu,alpha,remove_mean_drag,dt,Any[nothing,nothing,nothing],-1,-1,zeros(4),zeros(4),0.,0.,0.,0.,-Inf,Inf,-Inf,zeros(4),zeros(4),0.)

function (b::Audit)(m)
    o=b.op; s=m.clock.stage; it=m.clock.iteration
    if s==1
        b.iteration==it-1 && b.finished!=b.iteration && finish!(b,m)
        b.iteration=it; b.stages.=nothing
    end
    u,v=state(m);q=vorticity(o,u,v)
    gu=copy(dropdims(Array(interior(m.timestepper.Gⁿ.u));dims=3))
    gv=copy(dropdims(Array(interior(m.timestepper.Gⁿ.v));dims=3))
    native=vorticity(o,gu,gv)
    injection=band_curl(o,q,b.amplitude)
    viscous=spectral_curl(o,q,-b.nu.*o.lambda)
    fu,fv=velocity_from_curl(o,injection.+viscous)
    vc=(v.+circshift(v,(0,-1)))./2;vu=(vc.+circshift(vc,(1,0)))./2
    uc=(u.+circshift(u,(-1,0)))./2;uv=(uc.+circshift(uc,(0,1)))./2
    du=-b.alpha .* (u.^2 .+ vu.^2) .* u
    dv=-b.alpha .* (v.^2 .+ uv.^2) .* v
    # Keep the zero-mean subspace: spatially varying drag otherwise creates a
    # uniform velocity that undergoes f-dependent inertial oscillations.
    # For zero-mean u,v this projection does not change the drag's total work.
    if b.remove_mean_drag;du.-=mean(du);dv.-=mean(dv);end
    drag=vorticity(o,du,dv)
    interior(m.timestepper.Gⁿ.u)[:,:,1].+=fu.+du
    interior(m.timestepper.Gⁿ.v)[:,:,1].+=fv.+dv
    # Work computed using compatible velocity/curl relation in the periodic box.
    psi=spectral_curl(o,q,map(x->x==0 ? 0.0 : -1/x,o.lambda))
    b.min_injection_work=min(b.min_injection_work,-mean(psi.*injection))
    b.max_viscous_work=max(b.max_viscous_work,-mean(psi.*viscous))
    b.max_drag_work=max(b.max_drag_work,mean(u.*du.+v.*dv))
    dx=o.lx/o.nx;dy=o.ly/o.ny
    div=(circshift(u,(-1,0)).-u)./dx.+(circshift(v,(0,-1)).-v)./dy
    b.max_divergence=max(b.max_divergence,maximum(abs,div))
    b.max_cfl=max(b.max_cfl,b.dt*(maximum(abs,u)/dx+maximum(abs,v)/dy))
    b.max_cfl<0.5 || error("CFL guard exceeded; no automatic timestep change")
    all(isfinite,q)||error("Nonfinite vorticity")
    b.stages[s]=(q=q,psi=psi,parts=(native,injection,viscous,drag))
    nothing
end

function finish!(b,m)
    b.iteration<0 && return
    b.finished==b.iteration && return
    m.clock.iteration==b.iteration+1 || return
    all(x->x!==nothing,b.stages)||error("Missing RK stage")
    u,v=state(m);qend=vorticity(b.op,u,v)
    gamma=(8/15,5/12,3/4);zeta=(0.,-17/60,-5/12)
    increments=zeros(4); eincrements=zeros(4)
    psiend=spectral_curl(b.op,qend,map(x->x==0 ? 0.0 : -1/x,b.op.lambda))
    for s=1:3
        q=b.stages[s].q;qn=s<3 ? b.stages[s+1].q : qend
        deltas=[b.dt.*(gamma[s].*b.stages[s].parts[k].+
                (s==1 ? zeros(size(q)) : zeta[s].*b.stages[s-1].parts[k])) for k=1:4]
        predicted=reduce(+,deltas)
        err=maximum(abs,qn.-q.-predicted)/max(maximum(abs,q),maximum(abs,qn),1e-12)
        b.max_increment_error=max(b.max_increment_error,err)
        err<1e-9||error("Native vorticity increment did not close: $err")
        weight=qn.^2 .+qn.*q.+q.^2
        for k=1:4;increments[k]+=mean(weight.*deltas[k]);end
        psin=s<3 ? b.stages[s+1].psi : psiend
        eweight=-(psin.+b.stages[s].psi)./2
        for k=1:4;eincrements[k]+=mean(eweight.*deltas[k]);end
    end
    storage=mean(qend.^3)-mean(b.stages[1].q.^3)
    scale=max(mean(abs.(qend).^3),mean(abs.(b.stages[1].q).^3),1e-18)
    err=abs(storage-sum(increments))/scale
    b.max_moment_error=max(b.max_moment_error,err)
    err<1e-9||error("Third moment did not close: $err")
    b.third.+=increments;b.absolute_third.+=abs.(increments)
    e0=-mean(b.stages[1].psi.*b.stages[1].q)/2
    e1=-mean(psiend.*qend)/2
    eerr=abs(e1-e0-sum(eincrements))/max(e0,e1,1e-20)
    b.max_energy_error=max(b.max_energy_error,eerr)
    eerr<1e-9||error("Discrete energy increments did not close: $eerr")
    b.energy.+=eincrements;b.absolute_energy.+=abs.(eincrements)
    b.finished=b.iteration
    nothing
end

function attach_audit!(sim,b)
    add_callback!(sim,b,IterationInterval(1);name=:audited_forcing,callsite=TendencyCallsite())
    add_callback!(sim,s->finish!(b,s.model),IterationInterval(1);name=:finish_budget)
end

function random_seed(o,seed;urms=0.02)
    rng=MersenneTwister(seed)
    q=spectral_curl(o,randn(rng,o.nx,o.ny),(o.k2.>=4^2).&(o.k2.<=10^2))
    u,v=velocity_from_curl(o,q);a=urms/sqrt(mean(u.^2 .+v.^2))
    a.*u,a.*v
end

function report(b)
    physical=sum(b.absolute_third[2:4])
    ratio=b.absolute_third[1]/max(physical,1e-30)
    Dict("max_relative_energy_increment_error"=>b.max_energy_error,
         "energy_increments"=>b.energy,"absolute_step_energy_increments"=>b.absolute_energy,
         "energy_terms"=>collect(TERMS),
         "max_relative_vorticity_increment_error"=>b.max_increment_error,
         "max_relative_third_moment_error"=>b.max_moment_error,
         "max_divergence"=>b.max_divergence,"max_cfl"=>b.max_cfl,
         "max_drag_work"=>b.max_drag_work,"max_viscous_work"=>b.max_viscous_work,
         "min_injection_work"=>b.min_injection_work,
         "third_moment_terms"=>collect(TERMS),"third_moment_increments"=>b.third,
         "absolute_step_third_moment_increments"=>b.absolute_third,
         "numerical_advection_fraction"=>ratio,
         "physical_attribution_gate"=>(ratio<0.1 ? "PASS" : "FAIL_NO_MECHANISM_CLAIM"))
end
