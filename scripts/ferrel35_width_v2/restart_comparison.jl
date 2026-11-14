using LinearAlgebra

# Psi = -integral(rho*v dz) has one extra length dimension relative to rho*v.
# Keep the existing native-flux tolerance; propagate it through the integral.
function streamfunction_restart_report(a,b,H;rtol=1e-9,atol=1e-10)
    H>0 || error("Positive depth required")
    fa,fb=a["native_rhov"],b["native_rhov"]
    pa,pb=a["mass_streamfunction"],b["mass_streamfunction"]
    nz=size(fa,2)
    size(fa)==size(fb) || error("Restart flux shape changed")
    size(pa)==size(pb)==(size(fa,1),nz+1) || error("Restart streamfunction shape changed")
    dflux=norm(fa-fb);dpsi=norm(pa-pb)
    # The cumulative sum operator has 2-norm <= nz; dz*nz = H.
    # Include a conservative floating-point summation allowance, not a fitted tolerance.
    roundoff=8eps(Float64)*nz*H*max(norm(fa),norm(fb))
    propagated_bound=H*dflux+roundoff
    native_pass=isapprox(fa,fb;rtol,atol)
    normalized_pass=isapprox(pa/H,pb/H;rtol,atol)
    propagated_pass=dpsi<=propagated_bound
    return Dict("passed"=>native_pass && normalized_pass && propagated_pass,
      "native_flux_pass"=>native_pass,"depth_normalized_pass"=>normalized_pass,
      "propagated_error_pass"=>propagated_pass,"depth_m"=>H,
      "native_flux_l2_difference"=>dflux,"streamfunction_l2_difference"=>dpsi,
      "streamfunction_relative_l2_difference"=>dpsi/max(norm(pa),norm(pb),eps()),
      "propagated_bound"=>propagated_bound,"summation_roundoff_allowance"=>roundoff,
      "native_flux_atol"=>atol,"rtol"=>rtol)
end
