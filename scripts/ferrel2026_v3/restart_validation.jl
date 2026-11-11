using Test, Printf

"""Validate a derived mass streamfunction after strict physical-field restart checks.

The caller independently retains rtol=1e-9, atol=1e-10 on all twelve physical/native
fields, including rho*v. Psi is an integral in different units, so its allowed
difference is determined from the measured mass-flux difference and Float64
summation roundoff. This function does not replace the physical-field checks.
"""
function derived_streamfunction_restart_report(a,b,metadata)
    va,vb=a["native_rhov"],b["native_rhov"]
    pa,pb=a["mass_streamfunction"],b["mass_streamfunction"]
    H=Float64(metadata["height"]);nz=Int(metadata["nz"])
    isfinite(H) && H>0 && nz>0 || error("Invalid streamfunction integration geometry")
    dz=H/nz
    size(va)==size(vb) || error("Native mass-flux shapes differ")
    size(va,2)==nz || error("Native mass-flux shape disagrees with metadata")
    size(pa)==size(pb)==(size(va,1),nz+1) || error("Incorrect streamfunction staggering")
    all(x->eltype(x)==Float64,(va,vb,pa,pb)) || error("Roundoff envelope requires Float64 fields")
    all(x->all(isfinite,x),(va,vb,pa,pb)) || error("Nonfinite native mass flux or streamfunction")

    delta_v=vb.-va
    observed_delta=pb.-pa
    predicted_delta=zeros(size(pa))
    predicted_delta[:,2:end].=-dz.*cumsum(delta_v;dims=2)
    identity_residual=maximum(abs,observed_delta.-predicted_delta)

    # Units: each term in this envelope is kg/(m*s), the units of Psi.
    # The absolute mass-flux throughput guards against cancellation inside a
    # column. 8*(Nz+2)*epsilon covers two cumulative sums, products and differences.
    mass_throughput=dz*maximum(sum(abs.(va).+abs.(vb);dims=2))
    psi_scale=max(maximum(abs,pa),maximum(abs,pb),mass_throughput,floatmin(Float64))
    roundoff_bound=8*(nz+2)*eps(Float64)*psi_scale
    max_delta_v=maximum(abs,delta_v)
    integration_bound=H*max_delta_v
    observed_max=maximum(abs,observed_delta)

    pointwise_bound=zeros(size(pa))
    pointwise_bound[:,2:end].=dz.*cumsum(abs.(delta_v);dims=2)
    roundtrip_errors=Float64[]
    for (v,p) in ((va,pa),(vb,pb))
        reconstructed=zeros(size(p))
        reconstructed[:,2:end].=-dz.*cumsum(v;dims=2)
        push!(roundtrip_errors,maximum(abs,p.-reconstructed))
    end

    checks=Dict{String,Bool}(
        "difference_identity_pass"=>(identity_residual<=roundoff_bound),
        "whole_depth_bound_pass"=>(observed_max<=integration_bound+roundoff_bound),
        "pointwise_integral_bound_pass"=>all(abs.(observed_delta).<=pointwise_bound.+roundoff_bound),
        "individual_roundtrip_pass"=>(maximum(roundtrip_errors)<=roundoff_bound),
        "bottom_boundary_pass"=>(maximum(abs,pa[:,1])==0 && maximum(abs,pb[:,1])==0))
    return Dict{String,Any}(
        "checks"=>checks,"all_conditions_pass"=>all(values(checks)),
        "criterion"=>"Strict physical-field checks plus measured native-mass-flux integral consistency; no physical tolerance relaxed",
        "height_m"=>H,"dz_m"=>dz,"nz"=>nz,
        "physical_fields_rtol"=>1e-9,"physical_fields_atol"=>1e-10,
        "max_mass_flux_difference_kg_m2_s"=>max_delta_v,
        "max_streamfunction_difference_kg_m_s"=>observed_max,
        "integration_bound_H_max_delta_mass_flux"=>integration_bound,
        "roundoff_bound_kg_m_s"=>roundoff_bound,
        "derived_streamfunction_bound_kg_m_s"=>integration_bound+roundoff_bound,
        "difference_identity_max_residual_kg_m_s"=>identity_residual,
        "individual_roundtrip_max_errors"=>roundtrip_errors)
end


function check_derived_streamfunction_restart(a,b,metadata)
    report=derived_streamfunction_restart_report(a,b,metadata)
    for name in sort(collect(keys(report["checks"])))
        passed=report["checks"][name]
        @testset "$name" begin
            @test passed
        end
    end
    @printf("Derived restart Psi: max delta=%.5e; H*max delta(rhov)=%.5e; roundoff=%.5e; integral residual=%.5e\n",
        report["max_streamfunction_difference_kg_m_s"],report["integration_bound_H_max_delta_mass_flux"],
        report["roundoff_bound_kg_m_s"],report["difference_identity_max_residual_kg_m_s"])
    return report
end
