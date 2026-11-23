function guard_spec(cfg)
    Dict("mach_max"=>get(cfg,"mach_stop",.5),"temperature_min"=>get(cfg,"temperature_stop",.5),"theta_log_departure_max"=>get(cfg,"theta_log_departure_stop",log(2.)),"mean_temperature_drift_max"=>get(cfg,"profile_drift_stop",.05),"mass_error_max"=>1e-7)
end
function guard_reason(s,g)
    all(isfinite,values(s)) || return "NONFINITE"
    minimum(s[k] for k in ["density_min","rhotheta_min","pressure_min","temperature_min"])>0 || return "POSITIVITY"
    s["temperature_min"]>=g["temperature_min"] || return "TEMPERATURE_FLOOR"
    s["mach_max"]<g["mach_max"] || return "MACH"
    s["theta_log_departure"]<=g["theta_log_departure_max"] || return "THETA_DEPARTURE"
    s["mean_temperature_drift"]<=g["mean_temperature_drift_max"] || return "THERMAL_DRIFT"
    s["mass_error"]<=g["mass_error_max"] || return "MASS"
    "PASS"
end
function guard_checks(cfg)
    g=guard_spec(cfg);s=Dict("density_min"=>1.,"rhotheta_min"=>1.,"pressure_min"=>1.,"temperature_min"=>1.,"mach_max"=>.1,"theta_log_departure"=>0.,"mean_temperature_drift"=>0.,"mass_error"=>0.)
    @assert guard_reason(s,g)=="PASS"
    cases=[("density_min",-1.,"POSITIVITY"),("rhotheta_min",-1.,"POSITIVITY"),("pressure_min",-1.,"POSITIVITY"),("temperature_min",-1.,"POSITIVITY"),("temperature_min",.49,"TEMPERATURE_FLOOR"),("mach_max",.5,"MACH"),("theta_log_departure",1.,"THETA_DEPARTURE"),("mean_temperature_drift",.051,"THERMAL_DRIFT"),("mass_error",1e-6,"MASS"),("temperature_min",NaN,"NONFINITE")]
    for (key,value,expected) in cases;@assert guard_reason(merge(s,Dict(key=>value)),g)==expected;end
    Dict("status"=>"PASS","triggered_checks"=>length(cases),"limits"=>g)
end
function health(m,theta0,Tprofile0,mass0)
    a=arrays(m);rt=Array(interior(m.formulation.potential_temperature_density));P=Array(interior(m.dynamics.pressure))
    valid=minimum(a.rho)>0 && minimum(rt)>0 && minimum(a.T)>0
    Dict("density_min"=>minimum(a.rho),"rhotheta_min"=>minimum(rt),"pressure_min"=>minimum(P),"temperature_min"=>minimum(a.T),"mach_max"=>(valid ? maximum(sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T)) : Inf),"theta_log_departure"=>(valid ? maximum(abs,log.((rt./a.rho)./theta0)) : Inf),"mean_temperature_drift"=>maximum(abs,vec(mean(a.T;dims=(1,2)))./Tprofile0 .-1),"mass_error"=>abs(mean(a.rho)/mass0-1))
end
