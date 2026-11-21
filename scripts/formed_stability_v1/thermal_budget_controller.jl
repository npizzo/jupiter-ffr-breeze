# Feedback removes net imposed enthalpy input with a vertically uniform fractional
# theta tendency. That tendency alone leaves g*∂z(log(theta)) unchanged.
mutable struct ThermalBudgetController{A,T,P}
    rate::A
    theta::T
    profile::P
    rho::Vector{Float64}
    energy_weights::Vector{Float64}
    mask::Vector{Float64}
    target::Vector{Float64}
    dz::Float64
    tau::Float64
    bottom_heat::Float64
    bottom_theta::Float64
    enabled::Bool
    alpha::Float64
    previous_time::Float64
    previous_theta_rate::Float64
    theta_integral::Float64
    net_enthalpy_integral::Float64
    preupdate_net::Float64
    max_preupdate_net::Float64
end
function (c::ThermalBudgetController)(sim)
    compute!(c.theta);compute!(c.profile)
    th=vec(Array(interior(c.profile)))
    relaxation=-c.mask.*(th.-c.target)./c.tau
    cap_power=c.dz*sum(c.energy_weights.*relaxation)
    thermal_content=c.dz*sum(c.energy_weights.*th)
    theta_content=c.dz*sum(c.rho.*th)
    cap_theta=c.dz*sum(c.rho.*relaxation)
    now=Float64(sim.model.clock.time)
    net=c.bottom_heat+cap_power-c.alpha*thermal_content
    before_theta=c.bottom_theta+cap_theta-c.alpha*theta_content
    if c.previous_time>=0 && now>c.previous_time
        dt=now-c.previous_time
        c.theta_integral += dt*(c.previous_theta_rate+before_theta)/2
        # The post-update external enthalpy source is zero when feedback is enabled.
        c.net_enthalpy_integral += dt*net/2
        c.max_preupdate_net=max(c.max_preupdate_net,abs(net))
    end
    c.preupdate_net=net
    c.alpha=c.enabled ? (c.bottom_heat+cap_power)/thermal_content : 0.0
    copyto!(c.rate,[c.alpha])
    c.previous_time=now
    c.previous_theta_rate=c.bottom_theta+cap_theta-c.alpha*theta_content
    return nothing
end
function Oceananigans.prognostic_state(c::ThermalBudgetController)
    return (; alpha=c.alpha,previous_time=c.previous_time,previous_theta_rate=c.previous_theta_rate,
              theta_integral=c.theta_integral,net_enthalpy_integral=c.net_enthalpy_integral,
              preupdate_net=c.preupdate_net,max_preupdate_net=c.max_preupdate_net)
end
function Oceananigans.restore_prognostic_state!(c::ThermalBudgetController,s::NamedTuple)
    for key in keys(s);setproperty!(c,key,getproperty(s,key));end
    copyto!(c.rate,[c.alpha]);return c
end
