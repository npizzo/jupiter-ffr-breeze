function guard_spec(cfg)
    Dict("mach_max"=>get(cfg,"mach_stop",.5),"temperature_min"=>get(cfg,"temperature_stop",.5),"theta_log_departure_max"=>get(cfg,"theta_log_departure_stop",log(2.)),"mean_temperature_drift_max"=>get(cfg,"profile_drift_stop",.05),"mass_error_max"=>1e-7,"outside_mach_max"=>get(cfg,"outside_mach_stop",.5))
end
function guard_reason(s,g)
    all(isfinite,values(s)) || return "NONFINITE"
    minimum(s[k] for k in ["density_min","rhotheta_min","pressure_min","temperature_min"])>0 || return "POSITIVITY"
    s["temperature_min"]>=g["temperature_min"] || return "TEMPERATURE_FLOOR"
    s["mach_max"]<g["mach_max"] || return "MACH"
    get(s,"outside_mach",0.)<g["outside_mach_max"] || return "OUTSIDE_MACH"
    s["theta_log_departure"]<=g["theta_log_departure_max"] || return "THETA_DEPARTURE"
    s["mean_temperature_drift"]<=g["mean_temperature_drift_max"] || return "THERMAL_DRIFT"
    s["mass_error"]<=g["mass_error_max"] || return "MASS"
    "PASS"
end
function guard_checks(cfg)
    g=guard_spec(cfg);s=Dict("density_min"=>1.,"rhotheta_min"=>1.,"pressure_min"=>1.,"temperature_min"=>1.,"mach_max"=>.1,"theta_log_departure"=>0.,"mean_temperature_drift"=>0.,"mass_error"=>0.)
    @assert guard_reason(s,g)=="PASS"
    cases=[("density_min",-1.,"POSITIVITY"),("rhotheta_min",-1.,"POSITIVITY"),("pressure_min",-1.,"POSITIVITY"),("temperature_min",0.,"POSITIVITY"),("mach_max",g["mach_max"],"MACH"),("theta_log_departure",g["theta_log_departure_max"]+1.,"THETA_DEPARTURE"),("mean_temperature_drift",g["mean_temperature_drift_max"]+1.,"THERMAL_DRIFT"),("mass_error",1e-6,"MASS"),("temperature_min",NaN,"NONFINITE")]
    @assert guard_reason(merge(s,Dict("mach_max"=>g["mach_max"])),g)=="MACH"

    for (key,value,expected) in cases;@assert guard_reason(merge(s,Dict(key=>value)),g)==expected;end
    Dict("status"=>"PASS","triggered_checks"=>length(cases),"limits"=>g)
end
function health(m,theta0,Tprofile0,mass0,outside_mask=nothing)
    rho=interior(density_field(m));rt=interior(m.formulation.potential_temperature_density)
    T=interior(m.temperature);P=interior(pressure_field(m))
    uf=interior(m.velocities.u);vf=interior(m.velocities.v);wf=interior(m.velocities.w)
    u=.5 .* (uf.+circshift(uf,(-1,0,0)))
    v=.5 .* (vf[:,1:end-1,:].+vf[:,2:end,:]);w=.5 .* (wf[:,:,1:end-1].+wf[:,:,2:end])
    rmin=minimum(rho);rtmin=minimum(rt);Tmin=minimum(T);valid=rmin>0 && rtmin>0 && Tmin>0
    mach=valid ? sqrt.(u.^2 .+v.^2 .+w.^2)./sqrt.((5/3).*T) : zero(T)
    Dict("outside_mach"=>(outside_mask===nothing ? 0. : (valid ? maximum(mach.*outside_mask) : Inf)),"density_min"=>rmin,"rhotheta_min"=>rtmin,"pressure_min"=>minimum(P),"temperature_min"=>Tmin,
         "mach_max"=>(valid ? maximum(mach) : Inf),
         "theta_log_departure"=>(valid ? maximum(abs,log.((rt./rho)./theta0)) : Inf),
         "mean_temperature_drift"=>maximum(abs,vec(sum(T;dims=(1,2)))./(size(T,1)*size(T,2))./Tprofile0 .-1),"mass_error"=>abs((sum(rho)/length(rho))/mass0-1))
end
