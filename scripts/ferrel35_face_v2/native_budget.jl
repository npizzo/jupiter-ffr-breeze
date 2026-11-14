using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans: TendencyCallsite
using Oceananigans.Operators: ∂xᶠᶜᶜ
using Oceananigans.BoundaryConditions: fill_halo_regions!
import Breeze.AtmosphereModels as AM

# These are the exact explicit native operators used by the pinned Breeze model.
# Projection is separate; its x mean vanishes for periodic x and rho_ref(z).
@inline budget_advection(i,j,k,g,ρ,a,v,cl,cf,M,f,clock,fields,F,dyn) = -AM.x_momentum_flux_divergence(i,j,k,g,a,M,v,dyn)
@inline budget_coriolis(i,j,k,g,ρ,a,v,cl,cf,M,f,clock,fields,F,dyn) = -AM.x_f_cross_U(i,j,k,g,f,M)
@inline budget_diffusion(i,j,k,g,ρ,a,v,cl,cf,M,f,clock,fields,F,dyn) = -AM.∂ⱼ_𝒯₁ⱼ(i,j,k,g,ρ,cl,cf,clock,fields,nothing)
@inline budget_forcing(i,j,k,g,ρ,a,v,cl,cf,M,f,clock,fields,F,dyn) = F(i,j,k,g,clock,fields)
@inline budget_pressure(i,j,k,g,ρ,p) = -ρ[i,j,k]*∂xᶠᶜᶜ(i,j,k,g,p)

zonal_host(field) = (compute!(field); dropdims(Array(interior(field));dims=1))
zonal_reduction(field) = Field(Average(field;dims=(1,)))

mutable struct NativeBudget
    terms::Any
    total::Any
    momentum::Any
    pressure::Any
    pressure_source::Any
    output::String
    interval::Int
    first_steps::Int
    active_iteration::Int
    start_time::Float64
    start_state::Any
    stages::Vector{Any}
    stage_times::Vector{Float64}
    checks::Int
    max_absolute_residual::Float64
    max_relative_residual::Float64
    max_roundoff_units::Float64
    origin_state::Any
    origin_time::Float64
    integrated_parts::Any
    face::Any
end

function NativeBudget(m, out, cfg)
    isnothing(m.timestepper.implicit_solver) || error("Native budget requires fully explicit diffusion/advection")
    string(nameof(typeof(m.timestepper)))=="SSPRungeKutta3" || error("Budget quadrature requires SSPRK3")
    cfg["meridional_boundary"]=="channel" || error("Ferrel budget requires a channel")
    cfg["jet_speed"]==0 || error("No surrogate jet body force permitted")
    # This campaign has only ScalarDiffusivity and default zero-flux momentum BCs.
    # A different closure, drag, immersed boundary or momentum BC needs a new budget audit.
    string(nameof(typeof(m.closure)))=="ScalarDiffusivity" || error("Unaudited closure")
    args=(AM.continuity_density(m.dynamics), m.advection.momentum, m.velocities,
          m.closure,m.closure_fields,AM.advecting_momentum(m),m.coriolis,m.clock,
          Oceananigans.fields(m),m.forcing.ρu,AM.slow_dynamics(m.timestepper,m.dynamics))
    funcs=(budget_advection,budget_coriolis,budget_diffusion,budget_forcing)
    terms=map(f->zonal_reduction(KernelFunctionOperation{Face,Center,Center}(f,m.grid,args...)),funcs)
    pressure=zonal_reduction(KernelFunctionOperation{Face,Center,Center}(budget_pressure,m.grid,args[1],m.dynamics.pressure_anomaly))
    path=joinpath(out,"native_budget");mkpath(path)
    open(joinpath(path,"checks.csv"),"w") do io
        println(io,"iteration,start_time,end_time,dt,absolute_residual,relative_to_leading_increment,roundoff_units,max_rhs_partition_residual,max_zonal_pressure_gradient")
    end
    NativeBudget(terms,zonal_reduction(m.timestepper.Gⁿ.ρu),zonal_reduction(m.momentum.ρu),pressure,m.dynamics.pressure_anomaly,
        path,get(cfg,"native_budget_interval",250),get(cfg,"native_budget_first_steps",8),
        -1,0.0,nothing,Any[nothing,nothing,nothing],fill(NaN,3),0,0.0,0.0,0.0,nothing,NaN,nothing,nothing)
end

function native_terms(b::NativeBudget)
    parts=cat(zonal_host.(b.terms)...;dims=3)
    total=zonal_host(b.total)
    residual=maximum(abs,dropdims(sum(parts;dims=3);dims=3).-total)
    scale=max(maximum(abs,parts),maximum(abs,total),1e-15)
    residual <= 2e-11*scale+1e-14 || error("Native RHS partition mismatch: $residual / $scale")
    # p'/rho is not a prognostic field; direct initialization does not enter the
    # model's ordinary prognostic-halo refresh. Projection fills it after a solve.
    # Refresh it explicitly before any diagnostic stencil, including t=0/restart.
    fill_halo_regions!(b.pressure_source)
    pg=maximum(abs,zonal_host(b.pressure))
    pg < 1e-10 || error("Zonal pressure projection no longer cancels: $pg")
    face = b.face===nothing ? nothing : face_terms(b.face)
    return (;parts,total,residual,pg,face)
end

function initialize_pressure_diagnostic!(m)
    pressure=m.dynamics.pressure_anomaly
    density=AM.continuity_density(m.dynamics)
    gradient=zonal_reduction(KernelFunctionOperation{Face,Center,Center}(budget_pressure,m.grid,density,pressure))
    before=maximum(abs,zonal_host(gradient))
    fill_halo_regions!(pressure)
    after=maximum(abs,zonal_host(gradient))
    after<1e-10 || error("Initialized native zonal pressure gradient remains nonzero after halo fill: $after")
    return (;before,after)
end

function finish_budget_step!(b,m)
    all(x->x!==nothing,b.stages) || error("Missing native SSPRK3 stage")
    final_state=zonal_host(b.momentum)
    dt=Float64(m.clock.time)-b.start_time
    dt>0 || error("Nonpositive sampled step duration")
    isapprox(b.stage_times[2],b.start_time+dt;atol=1e-7,rtol=1e-12) || error("Unexpected SSPRK3 second-stage time")
    isapprox(b.stage_times[3],b.start_time+dt/2;atol=1e-7,rtol=1e-12) || error("Unexpected SSPRK3 third-stage time")
    integrated_parts=dt.*(b.stages[1].parts./6 .+ b.stages[2].parts./6 .+ 2b.stages[3].parts./3)
    predicted=dropdims(sum(integrated_parts;dims=3);dims=3)
    observed=final_state.-b.start_state
    residual=observed.-predicted
    absolute=maximum(abs,residual)
    leading=max(maximum(abs,observed),maximum(abs,predicted),maximum(abs,integrated_parts),1e-30)
    state_scale=max(maximum(abs,b.start_state),maximum(abs,final_state),1.0)
    roundoff=eps(Float64)*state_scale
    relative=absolute/max(leading,4096roundoff)
    units=absolute/roundoff
    # Independent state-increment check, including the effect of projection and RK mixing.
    absolute <= 4096roundoff+1e-9leading || error("Native zonal momentum closure failed: residual=$absolute, leading=$leading, roundoff=$roundoff")
    partition=maximum(s.residual for s in b.stages)
    pressure=maximum(s.pg for s in b.stages)
    open(joinpath(b.output,"checks.csv"),"a") do io
        println(io,join((b.active_iteration,b.start_time,m.clock.time,dt,absolute,relative,units,partition,pressure),','))
    end
    if b.active_iteration % b.interval == 0 || b.checks < 8
    jldsave(joinpath(b.output,@sprintf("step_%09d.jld2",b.active_iteration));
        iteration=b.active_iteration,start_time=b.start_time,end_time=Float64(m.clock.time),dt,
        momentum_start=b.start_state,momentum_end=final_state,observed,predicted,residual,
        integrated_parts,term_names=["advection","coriolis","diffusion","external_forcing"],
        absolute_residual=absolute,relative_to_leading_increment=relative,roundoff_units=units,
        stage_times=copy(b.stage_times),max_rhs_partition_residual=partition,max_zonal_pressure_gradient=pressure)
    end
    b.face===nothing || integrate_face!(b.face,b.stages,dt,m)
    b.integrated_parts .+= integrated_parts
    b.checks+=1;b.max_absolute_residual=max(b.max_absolute_residual,absolute)
    b.max_relative_residual=max(b.max_relative_residual,relative);b.max_roundoff_units=max(b.max_roundoff_units,units)
    b.active_iteration=-1
    return nothing
end

function (b::NativeBudget)(m)
    iteration=m.clock.iteration;stage=m.clock.stage
    if stage==1
        if b.active_iteration>=0 && iteration==b.active_iteration+1
            finish_budget_step!(b,m)
        elseif b.active_iteration>=0 && iteration!=b.active_iteration
            error("Skipped sampled momentum-budget step")
        end
        if true # Audit and integrate every complete step; save detailed stage samples sparsely.
            b.active_iteration=iteration;b.start_time=Float64(m.clock.time)
            b.start_state=zonal_host(b.momentum)
            if b.origin_state===nothing
                b.origin_state=copy(b.start_state);b.origin_time=b.start_time
                b.integrated_parts=zeros(size(b.start_state)...,4)
            end
            b.stages=Any[nothing,nothing,nothing];b.stage_times.=NaN
        end
    end
    if iteration==b.active_iteration
        stage in 1:3 || error("Unexpected SSPRK3 stage $stage")
        b.stages[stage]=native_terms(b);b.stage_times[stage]=Float64(m.clock.time)
    end
    return nothing
end

# Diagnostic observers deliberately restart with an empty sample. No physical state is altered.
Oceananigans.prognostic_state(::NativeBudget)=nothing
Oceananigans.restore_prognostic_state!(b::NativeBudget,::Nothing)=b


function flush_native_budget!(b,m)
    if b.active_iteration>=0 && m.clock.iteration==b.active_iteration+1
        finish_budget_step!(b,m)
    elseif b.active_iteration>=0
        # update_state! may already have evaluated the next stage-1 tendency at
        # this boundary. It has zero duration and must not enter the integral.
        m.clock.iteration==b.active_iteration && Float64(m.clock.time)==b.start_time &&
            b.stages[2]===nothing && b.stages[3]===nothing ||
            error("Cannot save budget during an unfinished RK step")
    end
    return nothing
end
function cumulative_native_budget(b,m)
    flush_native_budget!(b,m)
    current=zonal_host(b.momentum)
    initial=b.origin_state===nothing ? current : b.origin_state
    parts=b.integrated_parts===nothing ? zeros(size(current)...,4) : b.integrated_parts
    t0=b.origin_state===nothing ? Float64(m.clock.time) : b.origin_time
    residual=current.-initial.-dropdims(sum(parts;dims=3);dims=3)
    scale=max(maximum(abs,parts),maximum(abs,current.-initial),1e-30)
    allowance=8192eps(Float64)*max(maximum(abs,current),maximum(abs,initial),1.0)*max(b.checks,1)
    maximum(abs,residual)<=allowance+1e-8scale || error("Cumulative native momentum budget failed")
    return (;native_budget_origin_time=t0,native_budget_end_time=Float64(m.clock.time),
        native_budget_steps=b.checks,native_rhou_origin=initial,native_rhou=current,
        native_integrated_parts=copy(parts),native_integrated_residual=residual,
        native_integrated_relative_residual=maximum(abs,residual)/max(scale,allowance),
        native_term_names=["advection","coriolis","diffusion","external_forcing"])
end
