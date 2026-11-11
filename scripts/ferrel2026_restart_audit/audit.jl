# Read-only analysis of existing Unity output. Run inside a Unity CPU allocation.
using JLD2, TOML, LinearAlgebra, Printf, SHA

length(ARGS)==2 || error("usage: audit.jl CPU_SMOKE_ROOT NEW_AUDIT_OUTPUT")
root,out=abspath.(ARGS)
isdir(out) && !isempty(readdir(out)) && error("Audit output must be fresh")
mkpath(out)
first_file=joinpath(root,"continuous","frame_00002.jld2")
restart_file=joinpath(root,"resumed","frame_00002.jld2")
a,b=load(first_file),load(restart_file)
metadata=TOML.parsefile(joinpath(root,"continuous","metadata.toml"))
H=Float64(metadata["height"]);nz=Int(metadata["nz"]);dz=H/nz
keys_to_check=("u","v","w","theta_profile","vorticity_slices","tracer_slices","ubar","vbar",
    "native_rhou","native_rhov","native_rhow","mass_streamfunction","native_rhou_rhs")
rtol,atol=1e-9,1e-10
physical_pass=Ref(true)
comparison=Dict{String,Any}()
open(joinpath(out,"field_errors.csv"),"w") do io
    println(io,"field,max_abs_reference,max_abs_restart,max_abs_error,norm_reference,norm_restart,norm_error,relative_norm_error,original_norm_limit,original_isapprox")
    for key in keys_to_check
        x,y=a[key],b[key]
        size(x)==size(y) || error("Shape mismatch for $key")
        error_norm=norm(y.-x);scale=max(norm(x),norm(y))
        max_error=maximum(abs,y.-x)
        pass=isapprox(x,y;rtol,atol)
        key=="mass_streamfunction" || (physical_pass[] &= pass)
        info=Dict("max_abs_reference"=>maximum(abs,x),"max_abs_restart"=>maximum(abs,y),
            "max_abs_error"=>max_error,"norm_reference"=>norm(x),"norm_restart"=>norm(y),
            "norm_error"=>error_norm,"relative_norm_error"=>error_norm/max(scale,floatmin(Float64)),
            "original_norm_limit"=>max(atol,rtol*scale),"original_isapprox"=>pass)
        comparison[key]=info
        println(io,join((key,info["max_abs_reference"],info["max_abs_restart"],max_error,norm(x),norm(y),
            error_norm,info["relative_norm_error"],info["original_norm_limit"],pass),','))
        @printf("%-23s maxerr=%12.5e normerr=%12.5e relnorm=%12.5e original=%s\n",
            key,max_error,error_norm,info["relative_norm_error"],pass)
    end
end

va,vb=a["native_rhov"],b["native_rhov"]
pa,pb=a["mass_streamfunction"],b["mass_streamfunction"]
size(va,2)==nz || error("Metadata/grid mismatch")
size(pa)==(size(va,1),nz+1) || error("Incorrect streamfunction staggering")
delta_v=vb.-va
predicted_delta=zeros(size(pa))
predicted_delta[:,2:end].=-dz.*cumsum(delta_v;dims=2)
observed_delta=pb.-pa
identity_residual=maximum(abs,observed_delta.-predicted_delta)

# All terms below have streamfunction units kg/(m*s). The envelope covers
# independent cumulative sums, multiplications and subtractions in Float64.
# Absolute mass-flux throughput also covers cancellation within each column.
mass_throughput=dz*maximum(sum(abs.(va).+abs.(vb);dims=2))
psi_scale=max(maximum(abs,pa),maximum(abs,pb),mass_throughput,floatmin(Float64))
roundoff_bound=8*(nz+2)*eps(Float64)*psi_scale
max_delta_v=maximum(abs,delta_v)
integration_bound=H*max_delta_v
derived_bound=integration_bound+roundoff_bound
observed_max=maximum(abs,observed_delta)
identity_pass=identity_residual<=roundoff_bound
derived_pass=observed_max<=derived_bound
per_column_bound=zeros(size(pa))
per_column_bound[:,2:end].=dz.*cumsum(abs.(delta_v);dims=2)
pointwise_bound_pass=all(abs.(observed_delta).<=per_column_bound.+roundoff_bound)

roundtrip=Float64[]
for (v,p) in ((va,pa),(vb,pb))
    reconstructed=zeros(size(p));reconstructed[:,2:end].=-dz.*cumsum(v;dims=2)
    push!(roundtrip,maximum(abs,p.-reconstructed))
end
roundtrip_pass=maximum(roundtrip)<=roundoff_bound
passed=physical_pass[] && identity_pass && derived_pass && pointwise_bound_pass && roundtrip_pass
summary=Dict{String,Any}(
    "status"=>(passed ? "PASS" : "FAIL"),"source_directory"=>root,
    "reference_file"=>first_file,"restart_file"=>restart_file,
    "height_m"=>H,"dz_m"=>dz,"nz"=>nz,
    "physical_fields_original_tolerance_pass"=>physical_pass[],
    "original_rtol"=>rtol,"original_atol"=>atol,
    "streamfunction_original_norm_test_pass"=>comparison["mass_streamfunction"]["original_isapprox"],
    "max_mass_flux_difference_kg_m2_s"=>max_delta_v,
    "max_streamfunction_difference_kg_m_s"=>observed_max,
    "integration_bound_H_max_delta_mass_flux"=>integration_bound,
    "roundoff_bound_kg_m_s"=>roundoff_bound,
    "derived_streamfunction_bound_kg_m_s"=>derived_bound,
    "difference_identity_max_residual_kg_m_s"=>identity_residual,
    "difference_identity_pass"=>identity_pass,"derived_bound_pass"=>derived_pass,
    "pointwise_integral_bound_pass"=>pointwise_bound_pass,
    "individual_roundtrip_max_errors"=>roundtrip,"individual_roundtrip_pass"=>roundtrip_pass,
    "audit_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "field_comparison"=>comparison)
open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,summary);end
@printf("Streamfunction observed=%12.5e H*max(delta_rhov)=%12.5e roundoff=%12.5e identity_residual=%12.5e\n",
    observed_max,integration_bound,roundoff_bound,identity_residual)
println("Restart derived-field audit: ",summary["status"])
passed || error("Restart discrepancy is not explained by the measured native mass-flux differences")
