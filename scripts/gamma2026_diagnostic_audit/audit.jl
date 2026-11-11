# Read-only diagnostic audit. Run only in an allocated Unity CPU job.
using TOML, SHA, JLD2, FFTW, Statistics, LinearAlgebra, Printf
haskey(ENV,"SLURM_JOB_ID") || error("An allocated Unity job is required")
length(ARGS)==4 || error("usage: audit.jl ORIGINAL_CONFIG COMPLETE_CASE PROBE_CSV NEW_OUTPUT")
config_path,case_root,probe_csv,out=abspath.(ARGS)
ispath(out) && error("Use a fresh audit output directory")
repo=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(repo,"scripts","gamma2026_v3","core.jl"))
using .Gamma2026V3, GeophysicalFlows
const QG=Gamma2026V3.QG
sha(path)=bytes2hex(sha256(read(path)))
const EXPECTED_CORE="4f42089927e37b6ecb75de4d4b0bf337930550b2a7f0b83cfa9f89fb5bbbc6cc"
Gamma2026V3.source_hash()==EXPECTED_CORE || error("Unexpected frozen model core")

function read_numeric_csv(path)
    lines=readlines(path)
    header=String.(split(lines[1],','))
    rows=Dict{String,Float64}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values=parse.(Float64,split(line,','))
        length(values)==length(header) || error("Malformed diagnostic CSV")
        all(isfinite,values) || error("Nonfinite CSV diagnostic")
        push!(rows,Dict(zip(header,values)))
    end
    isempty(rows) && error("Empty diagnostic CSV")
    return rows
end

function norms(a)
    Dict("l2"=>norm(a),"max_abs"=>maximum(abs,a),"sum_abs2"=>sum(abs2,a))
end
function difference(a,b)
    size(a)==size(b) || error("Array shape mismatch")
    delta=a.-b
    Dict("max_abs"=>maximum(abs,delta),"l2"=>norm(delta),
         "relative_l2"=>norm(delta)/max(norm(b),1e-300),
         "relative_max"=>maximum(abs,delta)/max(maximum(abs,b),1e-300))
end
relative(a,b)=abs(a-b)/max(abs(b),1e-300)

# Independent real-FFT Parseval formula: first/last kx rows once, interior twice.
# The sum over ky is complete. This routine does not call FourierFlows.parsevalsum2.
function half_mean_square_spectral(h,nx,ny)
    iseven(nx) || error("Audit currently expects an even x grid")
    size(h)==(nx÷2+1,ny) || error("Unexpected real-FFT shape")
    weighted=sum(abs2,@view(h[1,:]))+sum(abs2,@view(h[end,:]))+
             2sum(abs2,@view(h[2:end-1,:]))
    return 0.5weighted/(nx*ny)^2
end

function diagnostic_at_state(p,c;label)
    s=Gamma2026V3.scales(c); factor=s.U/s.L
    before=copy(Array(p.sol));time_before=p.clock.t;step_before=p.clock.step
    QG.updatevars!(p)
    after=copy(Array(p.sol))
    q=Float64.(Array(p.vars.q));eta=Float64.(Array(p.params.eta))
    etah=copy(Array(p.params.etah))
    q_inverse=irfft(after,p.grid.nx)
    eta_inverse=irfft(etah,p.grid.nx)
    eta_direct_transform=rfft(eta)
    native_pv=q.+eta
    before_enstrophy=copy(Array(p.sol));before_eta=copy(Array(p.params.etah))
    model_pv=Float64(QG.enstrophy(p))*factor^2
    after_enstrophy=copy(Array(p.sol))
    model_energy=Float64(QG.energy(p))*s.U^2
    pv_real=0.5mean(abs2,native_pv)*factor^2
    pv_fourier=half_mean_square_spectral(after.+etah,p.grid.nx,p.grid.ny)*factor^2
    relative_real=0.5mean(abs2,q)*factor^2
    relative_fourier=half_mean_square_spectral(after,p.grid.nx,p.grid.ny)*factor^2
    eta_real=0.5mean(abs2,eta)*factor^2
    eta_fourier=half_mean_square_spectral(etah,p.grid.nx,p.grid.ny)*factor^2
    before_metrics=copy(Array(p.sol))
    repeated=Gamma2026V3.metrics(p)
    after_metrics=copy(Array(p.sol))
    p.clock.t==time_before && p.clock.step==step_before || error("Audit unexpectedly advanced the model")
    result=Dict{String,Any}(
        "label"=>label,"time_days"=>p.clock.t*s.T/Gamma2026V3.DAY,"step"=>p.clock.step,
        "state_before_updatevars"=>norms(before),"state_after_updatevars"=>norms(after),
        "updatevars_state_difference"=>difference(after,before),
        "enstrophy_call_state_difference"=>difference(after_enstrophy,before_enstrophy),
        "enstrophy_call_background_difference"=>difference(Array(p.params.etah),before_eta),
        "metrics_call_state_difference"=>difference(after_metrics,before_metrics),
        "vars_q_vs_independent_irfft"=>difference(q,q_inverse),
        "eta_real_vs_independent_irfft"=>difference(eta,eta_inverse),
        "eta_saved_transform_vs_independent_rfft"=>difference(etah,eta_direct_transform),
        "model_energy_m2_s2"=>model_energy,
        "model_pv_enstrophy_s2"=>model_pv,"independent_real_pv_enstrophy_s2"=>pv_real,
        "independent_spectral_pv_enstrophy_s2"=>pv_fourier,
        "independent_real_relative_enstrophy_s2"=>relative_real,
        "independent_spectral_relative_enstrophy_s2"=>relative_fourier,
        "background_real_enstrophy_s2"=>eta_real,"background_spectral_enstrophy_s2"=>eta_fourier,
        "model_vs_independent_spectral_pv_relative_difference"=>relative(model_pv,pv_fourier),
        "model_vs_independent_real_pv_relative_difference"=>relative(model_pv,pv_real),
        "independent_real_vs_spectral_pv_relative_difference"=>relative(pv_real,pv_fourier),
        "background_real_vs_spectral_relative_difference"=>relative(eta_real,eta_fourier),
        "metrics_pv_enstrophy_s2"=>repeated.Z*factor^2,
        "metrics_relative_enstrophy_s2"=>repeated.relative_Z*factor^2)
    return result,q,eta
end

function saved_snapshot_analysis(path,c,eta,solver_rows;native_q=nothing)
    d=JLD2.load(path)
    Int(d["sampling_stride"])==1 || error("Native snapshots required")
    zeta=Float64.(d["vorticity_s_inv"]);pv=Float64.(d["pv_s_inv"])
    size(zeta)==size(eta)==size(pv) || error("Snapshot shape mismatch")
    all(isfinite,zeta) && all(isfinite,pv) || error("Nonfinite snapshot")
    day=Float64(d["time_seconds"])/Gamma2026V3.DAY
    match=findall(row->abs(row["time_days"]-day)<1e-7,solver_rows)
    length(match)==1 || error("No unique solver diagnostic at snapshot time")
    solver=solver_rows[only(match)]
    s=Gamma2026V3.scales(c);factor=s.U/s.L
    pv_saved=0.5mean(abs2,pv);zeta_saved=0.5mean(abs2,zeta)
    # This reconstructs the exact Float64 background but retains exported-zeta error.
    mixed_pv=zeta.+eta.*factor
    mixed_real=0.5mean(abs2,mixed_pv)
    mixed_spectral=half_mean_square_spectral(rfft(mixed_pv),c["n"],c["n"])
    report=Dict{String,Any}(
        "file"=>basename(path),"sha256"=>sha(path),"time_days"=>day,
        "solver_pv_enstrophy_s2"=>solver["pv_enstrophy_s2"],
        "solver_relative_enstrophy_s2"=>solver["relative_enstrophy_s2"],
        "saved_float32_pv_reduced_in_float64_s2"=>pv_saved,
        "saved_float32_relative_enstrophy_reduced_in_float64_s2"=>zeta_saved,
        "saved_pv_vs_solver_relative_difference"=>relative(pv_saved,solver["pv_enstrophy_s2"]),
        "saved_relative_enstrophy_vs_solver_relative_difference"=>relative(zeta_saved,solver["relative_enstrophy_s2"]),
        "float64_eta_plus_saved_zeta_real_enstrophy_s2"=>mixed_real,
        "float64_eta_plus_saved_zeta_spectral_enstrophy_s2"=>mixed_spectral,
        "mixed_real_vs_spectral_relative_difference"=>relative(mixed_real,mixed_spectral),
        "saved_pv_vs_mixed_fields"=>difference(pv,mixed_pv),
        "native_float64_state_available"=>native_q!==nothing)
    if native_q!==nothing
        # Recreate the exact expression in the frozen Float32 snapshot exporter.
        recreated_zeta=Float64.(Float32.(native_q).*Float32(factor))
        recreated_pv=Float64.((Float32.(native_q).+Float32.(eta)).*Float32(factor))
        native_pv=(native_q.+eta).*factor
        native_pv_z=0.5mean(abs2,native_pv)
        report["native_real_pv_enstrophy_s2"]=native_pv_z
        report["native_real_relative_enstrophy_s2"]=0.5mean(abs2,native_q)*factor^2
        report["recreated_export_pv_enstrophy_s2"]=0.5mean(abs2,recreated_pv)
        report["recreated_export_relative_enstrophy_s2"]=0.5mean(abs2,recreated_zeta)
        report["recreated_export_vs_native_pv_relative_difference"]=relative(0.5mean(abs2,recreated_pv),native_pv_z)
        report["saved_vs_recreated_pv_fields"]=difference(pv,recreated_pv)
        report["saved_vs_recreated_zeta_fields"]=difference(zeta,recreated_zeta)
        report["native_vs_solver_pv_relative_difference"]=relative(native_pv_z,solver["pv_enstrophy_s2"])
    end
    return report
end

c=TOML.parsefile(config_path)
metadata_path=joinpath(case_root,"metadata.toml")
m=TOML.parsefile(metadata_path)
m["status"]=="COMPLETE" || error("Audit requires an immutable completed parent case")
Gamma2026V3.config_hash(c)==m["config_hash"] || error("Config differs from parent metadata")
Gamma2026V3.source_hash()==m["source_hash"] || error("Core differs from parent metadata")
checkpoint=joinpath(case_root,"checkpoint.jld2")
before_hashes=Dict("checkpoint"=>sha(checkpoint),"metadata"=>sha(metadata_path),
                   "config"=>sha(config_path),"probe_csv"=>sha(probe_csv))
ck=JLD2.load(checkpoint)
ck["source_hash"]==Gamma2026V3.source_hash() || error("Checkpoint core mismatch")
ck["config_hash"]==Gamma2026V3.config_hash(c) || error("Checkpoint config mismatch")
ck["geophysicalflows_commit"]==Gamma2026V3.GF_COMMIT || error("Pinned solver differs")
solver_rows=read_numeric_csv(joinpath(case_root,"diagnostics.csv"))
probe_rows=read_numeric_csv(probe_csv)
worst_index=argmax([row["snapshot_pv_enstrophy_relative_error"] for row in probe_rows])
worst=probe_rows[worst_index]
println("Worst publication snapshot PV discrepancy at day ",worst["time_days"],
        ": ",worst["snapshot_pv_enstrophy_relative_error"]);flush(stdout)

p=Gamma2026V3.problem(c,CPU();initial=false)
raw=copy(ck["sol"])
parent_E0=Gamma2026V3.restore!(p,c,checkpoint)
after_restore=copy(Array(p.sol))
checkpoint_report,q_final,eta=diagnostic_at_state(p,c;label="Completed parent checkpoint reconstructed on CPU")
checkpoint_report["raw_checkpoint_state"]=norms(raw)
checkpoint_report["after_restore_state"]=norms(after_restore)
checkpoint_report["restore_projection_difference"]=difference(after_restore,raw)
checkpoint_report["original_initial_energy_nondimensional"]=parent_E0
checkpoint_report["interpretation"]="Checkpoint spectral state is exact input; physical-space reconstruction uses CPU FFTW, while the original model ran on a GPU. No time integration was performed."
final_frame=joinpath(case_root,"frames",@sprintf("frame_%06d.jld2",round(Int,m["completed_days"]/c["output_days"])))
final_snapshot=saved_snapshot_analysis(final_frame,c,eta,solver_rows;native_q=q_final)
abs(final_snapshot["time_days"]-checkpoint_report["time_days"])<1e-7 || error("Checkpoint and final frame times differ")

# Always check the initial state as a second independent state, and reconstruct
# it explicitly if it is the worst mismatch time. This is a CPU reconstruction,
# not a claim to possess the original GPU Float64 initial physical field.
p0=Gamma2026V3.problem(c,CPU())
initial_report,q_initial,eta_initial=diagnostic_at_state(p0,c;label="Seeded initial state independently reconstructed on CPU")
initial_report["interpretation"]="Frozen config and initializer; CPU FFTW reconstruction may differ from original GPU arithmetic at roundoff."
initial_frame=joinpath(case_root,"frames","frame_000000.jld2")
initial_snapshot=saved_snapshot_analysis(initial_frame,c,eta_initial,solver_rows;native_q=q_initial)
worst_frame=joinpath(case_root,"frames",@sprintf("frame_%06d.jld2",round(Int,worst["time_days"]/c["output_days"])))
if abs(worst["time_days"])<1e-7
    worst_snapshot=initial_snapshot
elseif abs(worst["time_days"]-m["completed_days"])<1e-7
    worst_snapshot=final_snapshot
else
    worst_snapshot=saved_snapshot_analysis(worst_frame,c,eta,solver_rows)
    worst_snapshot["limitation"]="No intermediate Float64 checkpoint exists at this time; exact original-state attribution is unavailable. The mixed reconstruction uses exported Float32 relative vorticity with the exact background."
end
after_hashes=Dict("checkpoint"=>sha(checkpoint),"metadata"=>sha(metadata_path),
                  "config"=>sha(config_path),"probe_csv"=>sha(probe_csv))
before_hashes==after_hashes || error("A parent input changed during the audit")
report=Dict{String,Any}(
    "status"=>"COMPLETE","audit_kind"=>"Read-only diagnostic comparison; no model/tolerance changes",
    "compute_location"=>"Unity allocated CPU job","slurm_job_id"=>ENV["SLURM_JOB_ID"],
    "case"=>m["name"],"n"=>c["n"],"source_hash"=>Gamma2026V3.source_hash(),
    "geophysicalflows_commit"=>Gamma2026V3.GF_COMMIT,"audit_sha256"=>sha(@__FILE__),
    "input_hashes"=>before_hashes,"parent_inputs_unchanged"=>true,"model_steps_executed"=>0,
    "probe_rows_scanned"=>length(probe_rows),"worst_probe_row"=>worst,
    "checkpoint_state_diagnostics"=>checkpoint_report,"final_snapshot_comparison"=>final_snapshot,
    "initial_state_diagnostics"=>initial_report,"initial_snapshot_comparison"=>initial_snapshot,
    "worst_snapshot_comparison"=>worst_snapshot,
    "interpretation_limit"=>"Measured discrepancies are reported without presuming Float32 rounding is the sole cause. This audit does not modify existing publication tolerances or assert physical conservation.")
mkpath(out)
open(joinpath(out,"audit.toml"),"w") do io;TOML.print(io,report);end
println("READ_ONLY_DIAGNOSTIC_AUDIT_COMPLETE ",out);flush(stdout)
