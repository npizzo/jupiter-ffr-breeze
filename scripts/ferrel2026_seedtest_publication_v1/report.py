"""Unity-only four-case horizontal thermal-seed report; no simulation changes."""
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Run this observer/rendering job on an allocated Unity node")

from guards import CASES,sha,read,require,validate

PINS={"resolution":("scripts/ferrel2026_resolution/observe.py","10c489bead71a4f305e3d1f8c00a53e8398e05375896c405548c484da2daf7cd"),
      "publication":("scripts/ferrel2026_publication_v2/render_report.py","ad5fc63bf40dcccdcc7ce7ac88608f3e4545f52bdb043d35a25add181b407697")}


def load_frozen(name):
    path,digest=PINS[name]
    require(sha(path)==digest,f"Frozen {name} helper changed")
    spec=importlib.util.spec_from_file_location(f"ferrel_seed_{name}",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


R=load_frozen("resolution");P=load_frozen("publication")
np,h5py,plt=R.np,R.h5py,R.plt
LABELS={"zero512_b0":"Zero θ seed · 512² · unheated",
        "zero512_b162":"Zero θ seed · 512² · heated",
        "modes512_b162":"Horizontal θ modes · 512² · heated",
        "modes1024_b162":"Horizontal θ modes · 1024² · heated"}
# Presentation labels only; imported numerical and rendering functions unchanged.
P.LABELS.update(LABELS)
COLORS=dict(zip(CASES,("#222222","#bc4b35","#2878b5","#7b4ba5")))
WINDOW=(0.,21600.)
PAIRS=(("zero512_b0","zero512_b162","zero_seed_heating","Heating response without an imposed thermal perturbation; includes any adjustment/heating interaction."),
       ("zero512_b162","modes512_b162","heated_seed_sensitivity","Thermal-seed sensitivity at fixed512 grid and heating; not a heating effect."),
       ("modes512_b162","modes1024_b162","heated_resolution_sensitivity","Horizontal-grid sensitivity for a shared thermal seed; velocity discretization also changes. Not fieldwise or asymptotic convergence."))
CAVEATS=[
    "All results cover0–6 Earth hours of startup only. No mature convection, steady Ferrel-cell or inverse-cascade claim.",
    "No unheated modal case exists: modal-minus-zero-seed changes cannot be attributed to heating.",
    "A horizontal band-limited seed can delay shorter unstable modes. A cleaner response does not establish resolution of intrinsic convective scales.",
    "Both grids retain96 vertical levels and short vertical seed correlations. Five horizontal planes do not establish vertical or bottom-boundary-layer resolution.",
    "Differences are single-realization diagnostics, not statistical significance. Tail fractions use each case's own cell width; their physical cutoffs differ across grids.",
    "Float64 anomaly output removes the old Float32 output loss; it does not remove finite-precision model dynamics or supply a bound on accumulated solver error.",
    "Heat cospectra and centered profile heat flux are resolved covariance diagnostics, not native WENO face fluxes or a closed thermal transport budget.",
    "Centered eddy-to-mean KE transfer is approximate, with a reported native-advection remainder. Its sign alone does not prove an inverse cascade.",
    "Native momentum closure validates sampled SSPRK3 steps, not a continuous long-run integrated budget.",
    "Hourly snapshots undersample cap buoyancy waves. No temporal wave/convection classification is attempted.",
    "The cap relaxation, global thermal correction and velocity sponge remain idealized. This is not a calibrated Jovian cooling law.",
]


def half_ulp64(a):
    a=np.asarray(a,dtype=np.float64)
    return .5*np.maximum(np.nextafter(a,np.inf)-a,a-np.nextafter(a,-np.inf))


def precision_spectrum(w,anomaly,theta_mean,m,mode):
    # The stored arrays themselves have no additional output quantization.
    # This is a separately declared one-rounding representation sensitivity:
    # Float64 absolute-theta ulp plus the anomaly subtraction's ulp. A constant
    # mean error is removed by either orthogonal projection. It is NOT a bound
    # on all accumulated errors in the simulated state.
    qw=half_ulp64(w)
    qt=half_ulp64(anomaly+theta_mean)+half_ulp64(anomaly)
    row,shells=R.spectrum(w,anomaly,qw,qt,m,mode)
    row={k.replace("quantization","representation"):v for k,v in row.items()}
    row["input_error_semantics"]="One-rounding Float64 representation sensitivity, not output loss or accumulated dynamical-error bound"
    return row,shells


def analytic_checks():
    checks=R.self_checks()
    nx=ny=64;lx=2*np.pi;ly=np.pi
    x=(np.arange(nx)+.5)*lx/nx;y=(np.arange(ny)+.5)*ly/ny
    X,Y=np.meshgrid(x,y);w=np.cos(2*X)*np.cos(Y);a=1e-7*w
    m={"lx":lx,"ly":ly};zero=np.zeros_like(w)
    for mode in ("plane","zonal"):
        row,_=precision_spectrum(w,a,170.,m,mode)
        require(row["heat_cospectrum_status"]=="BOUNDED_NOT_MODEWISE_CONVERGED","Float64 anomaly lost a resolved weak signal")
        bound=row["heat_cospectrum_l1_representation_bound_K_m_s"]+row["heat_cospectrum_l1_arithmetic_bound_K_m_s"]
        require(abs(row["heat_cospectrum_sum_K_m_s"]-2.5e-8)<=bound,"Weak anomaly covariance normalization/sign failed")
        empty,_=precision_spectrum(w,zero,170.,m,mode)
        require(empty["heat_cospectrum_status"]=="INDETERMINATE_PRECISION_FLOOR","Zero thermal signal must not acquire an inferred heat spectrum")
    valid=R.profile_check(w,a,zero,zero,0.,1.,1.,1.,.125,float(np.mean(w*a)))
    try:R.profile_check(w,a,zero,zero,0.,1.,1.,1.,.125,1e-2)
    except ValueError:pass
    else:raise ValueError("Anomaly profile check accepted corrupted heat flux")
    checks["checks"]+= ["Float64 anomaly retains a weak signal erased by absolute Float32 output",
                        "zero anomaly remains precision-indeterminate; direct anomaly/profile check rejects corruption"]
    return checks


def observe(name,folder,out):
    m=read(folder/"metadata.toml");nx,ny,nz=(m[k] for k in ("nx","ny","nz"))
    frames=sorted(folder.glob("frame_*.jld2"));require(len(frames)==7,"Require exactly seven hourly frames")
    times=[];records=[];checks=[];activities=[];spectra={};profiles=[];frame_records=[];regime=None
    for fi,path in enumerate(frames):
        with h5py.File(path) as d:
            t=float(R.finite(d,"t"));require(t==fi*3600.,"Missing or shifted hourly frame");times.append(t)
            z=R.finite(d,"z",(nz,));rho=R.finite(d,"density",(nz,));pressure=R.finite(d,"pressure",(nz,))
            theta_mean=R.finite(d,"theta_profile",(nz,));ke=R.finite(d,"vertical_ke_profile",(nz,));heat=R.finite(d,"heat_flux_profile",(nz,))
            diffusion=R.finite(d,"diffusive_heat_flux_profile",(nz,));n2=R.finite(d,"mean_N2_profile",(nz,))
            heights=R.finite(d,"precision_heights",(5,));stored_indices=R.finite(d,"precision_indices_julia_1based",(5,))
            require(np.issubdtype(stored_indices.dtype,np.integer),"Precision indices must be stored as integers")
            indices=stored_indices.astype(int)-1
            require(np.array_equal(indices,[np.argmin(abs(z-h)) for h in m["precision_diagnostic_heights"]]),"Precision planes do not match requested locations")
            require(len(set(indices))==5 and np.array_equal(indices+1,m["precision_diagnostic_indices_julia_1based"]),"Precision-plane index mismatch")
            require(np.array_equal(heights,z[indices]) and np.array_equal(heights,m["precision_diagnostic_heights_actual_m"]),"Precision-plane height mismatch")
            video=R.finite(d,"video_heights",(3,));vi=np.asarray([np.argmin(abs(z-h)) for h in m["video_heights_requested_m"]])
            require(np.array_equal(video,z[vi]) and np.array_equal(video,m["video_heights_actual_m"]),"Video-plane index mismatch")
            W=R.finite(d,"precision_wcenter_slices",(5,ny,nx));T=R.finite(d,"precision_theta_anomaly_slices",(5,ny,nx))
            require(all(a.dtype==np.dtype("float64") for a in (W,T,z,rho,pressure,theta_mean,ke,heat)),"Precision diagnostics require native Float64 storage")
            W32=R.finite(d,"wcenter_slices",(3,ny,nx));T32=R.finite(d,"theta_slices",(3,ny,nx))
            for j,k in enumerate(vi):
                q=int(np.flatnonzero(indices==k)[0])
                require(np.array_equal(W[q].astype(np.float32),W32[j]),"Precision/movie w planes differ")
                require(np.array_equal((T[q]+theta_mean[k]).astype(np.float32),T32[j]),"Precision/movie theta planes differ")
                if j==0:
                    absolute=R.finite(d,"theta",(ny,nx))
                    require(absolute.dtype==np.dtype("float64") and np.array_equal(T[q],absolute-theta_mean[k]),"Lowest Float64 absolute/anomaly planes differ")
            exner=(pressure/1e5)**((8.314462618/m["molar_mass"])/m["heat_capacity"])
            if regime is None:regime=R.nominal(m,rho,theta_mean,pressure,z)
            w2=2*ke/rho;active=int(np.argmax(w2));thermal=int(np.argmax(abs(heat)))
            profiles.append({"z":z,"w2":w2,"heat":heat,"diffusion":diffusion,"n2":n2,"theta":theta_mean})
            activities.append({"time_seconds":t,"height_of_max_w_variance_m":float(z[active]),"max_w_rms_m_s":float(np.sqrt(max(w2[active],0))),
                "height_of_max_abs_heat_flux_m":float(z[thermal]),"max_abs_resolved_heat_flux_W_m2":float(abs(heat[thermal])),
                "min_mean_N2_s2":float(n2.min()),"max_positive_mean_N2_s2":float(max(n2.max(),0)),
                "N_profile_max_times_configured_dt":float(np.sqrt(max(n2.max(),0))*m["max_dt"]),
                "scope":"Horizontal-mean N only; not a local wave/CFL stability guarantee"})
            fingerprint=hashlib.sha256()
            for a in (W,T,z,rho,pressure,theta_mean,ke,heat):fingerprint.update(np.asarray(a,dtype="<f8",order="C").tobytes())
            frame_records.append({"file":path.name,"time_seconds":t,"precision_and_profile_fields_sha256":fingerprint.hexdigest(),
                "encoding":"C-order little-endian Float64; W,Tanomaly,z,rho,pressure,theta_mean,vertical_ke,heat_flux; shapes recorded separately"})
            for level,k in enumerate(indices):
                w,anomaly=W[level],T[level];zero=np.zeros_like(w)
                # The saved plane is ALREADY demeaned using the model mean.
                # Passing theta_mean=0 avoids subtracting the absolute mean twice.
                check=R.profile_check(w,anomaly,zero,zero,0.,rho[k],m["heat_capacity"],exner[k],ke[k],heat[k])
                check.update(time_seconds=t,level=level,height_m=float(z[k]),theta_storage="Float64 theta minus the model's same-plane mean",comparison_scope="Stored-field identities with Float64 arithmetic allowance; no output quantization")
                checks.append(check)
                for mode in ("plane","zonal"):
                    row,shells=precision_spectrum(w,anomaly,theta_mean[k],m,mode)
                    row.update(time_seconds=t,earth_day=t/86400,level=level,height_m=float(z[k]),theta_storage="Float64 anomaly",
                               centered_profile_heat_flux_W_m2=float(heat[k]))
                    records.append(row);spectra.setdefault((level,mode),[]).append(shells)
        print("SEED_SCREENED",name,f"{t/3600:g}h",flush=True)
    times=np.asarray(times);summary={};mean_spectra={};spectral_rows=[];profile_rows=[]
    for key,seq in spectra.items():
        level,mode=key;k=seq[0]["k"];dk=seq[0]["dk"]
        avg={field:R.average([s[field] for s in seq],times,*WINDOW) for field in ("w_power","heat_signed","heat_absolute","power_error","heat_error")}
        mean_spectra[key]=(k,dk,avg)
        summary[f"{level}_{mode}"]={"height_m":float(heights[level]),"w_variance_m2_s2":float(avg["w_power"].sum()),
            "w_power_l1_error_allowance_m2_s2":float(avg["power_error"]),"signed_heat_covariance_K_m_s":float(avg["heat_signed"].sum()),
            "absolute_mode_heat_throughput_K_m_s":float(avg["heat_absolute"].sum()),"heat_l1_error_allowance_K_m_s":float(avg["heat_error"]),
            "heat_status":"INDETERMINATE_PRECISION_FLOOR" if avg["heat_absolute"].sum()<=avg["heat_error"] else "BOUNDED_NOT_MODEWISE_CONVERGED"}
        for j in range(1,len(k)):
            spectral_rows.append({"window":"startup_0_6h","level":level,"height_m":float(heights[level]),"demeaning":mode,
                "k_lower_rad_m":j*dk,"k_upper_rad_m":(j+1)*dk,"k_center_rad_m":float(k[j]),"wavelength_m":float(2*np.pi/k[j]),"inverse_wavenumber_m":float(1/k[j]),
                "w_shell_variance_m2_s2":float(avg["w_power"][j]),"signed_heat_shell_K_m_s":float(avg["heat_signed"][j]),
                "absolute_mode_heat_in_shell_K_m_s":float(avg["heat_absolute"][j]),"any_subset_power_error_allowance":float(avg["power_error"]),"any_subset_heat_error_allowance":float(avg["heat_error"])})
    for field in ("w2","heat","diffusion","n2","theta"):
        avg=R.average([p[field] for p in profiles],times,*WINDOW)
        profile_rows.extend({"window":"startup_0_6h","height_m":float(z[i]),"quantity":field,"value":float(value)} for i,value in enumerate(avg))
    for suffix,rows in (("scale_timeseries",[R.flatten(r) for r in records]),("profile_agreement",checks),("vertical_activity",activities),("window_spectra",spectral_rows),("window_profiles",profile_rows)):
        R.write_csv(out/f"{name}_{suffix}.csv",rows)
    endpoint=[r for r in records if r["time_seconds"]==21600 and r["demeaning"]=="zonal"]
    return {"name":name,"m":m,"heights":heights,"times":times,"records":records,"spectra":mean_spectra,"profiles":profiles,
            "summary":{"case":name,"label":LABELS[name],"source":str(folder),"grid":[nx,ny,nz],"frames":7,"window_seconds":list(WINDOW),
                "precision_heights_actual_m":heights.tolist(),"precision_field_shape_h5py":[5,ny,nx],"movie_heights_actual_m":m["video_heights_actual_m"],"nominal_regime":regime,
                "mean_spectra":summary,"final_zonal_removed_scale_records":endpoint,"frame_manifest":frame_records,
                "tail_scope":"Every hourly snapshot has4/8/16-cell radial and component tails in CSV. Fractions are not averaged across time or grids.",
                "representation_scope":"One-rounding Float64 input sensitivity plus arithmetic-model allowance; not accumulated solver error or physical convergence."}}


def figures(cases,native,out):
    heights=cases[0]["heights"]
    require(all(np.array_equal(c["heights"],heights) for c in cases),"Cannot overlay different diagnostic heights")
    fig,axes=plt.subplots(5,2,figsize=(13,17),layout="constrained")
    for level,h in enumerate(heights):
        for case in cases:
            name=case["name"];k,dk,a=case["spectra"][(level,"zonal")];take=k>0
            wavelength=2*np.pi/k[take]/1000;power=a["w_power"][take]*k[take]/dk
            axes[level,0].loglog(wavelength,np.where(power>0,power,np.nan),color=COLORS[name],label=LABELS[name])
            axes[level,1].semilogx(wavelength,a["heat_signed"][take]*k[take]/dk,color=COLORS[name],label=LABELS[name])
        for ax in axes[level]:
            for n,style in ((512,":"),(1024,"--")):
                for cells in (4,8,16):ax.axvline(cells*cases[0]["m"]["lx"]/n/1000,color="gray",linestyle=style,alpha=.3,linewidth=.8)
            ax.set_xlabel("horizontal wavelength [km]");ax.grid(alpha=.2)
        axes[level,0].set(title=f"z={h/1000:.3f} km · all modes retained",ylabel="w variance per log k [m²/s²]")
        axes[level,1].set(title="Signed heat covariance · consult error status",ylabel="w–θ covariance per log k [K m/s]")
        axes[level,1].axhline(0,color="black",linewidth=.5)
    axes[0,0].legend(fontsize=7)
    fig.suptitle("Horizontal seed tests · 0–6 h mean · zonal mean removed\n4/8/16-cell markers: dotted512, dashed1024; single-seed sensitivity, not convergence",fontsize=11)
    fig.savefig(out/"five_height_spectra.png",dpi=110);plt.close(fig)

    fig,axes=plt.subplots(5,3,figsize=(15,16),layout="constrained")
    for level,h in enumerate(heights):
        for case in cases:
            name=case["name"];rows=[r for r in case["records"] if r["level"]==level and r["demeaning"]=="zonal"]
            hours=[r["time_seconds"]/3600 for r in rows]
            peak=[np.nan if r["peak_shell_wavelength_m"] is None else r["peak_shell_wavelength_m"]/1000 for r in rows]
            axes[level,0].plot(hours,peak,color=COLORS[name],label=LABELS[name])
            for col,metric in ((1,"w_variance"),(2,"absolute_heat_cospectrum")):
                tails=[r["tail_fractions"]["radial_lt_8cells"][metric] for r in rows]
                value=[np.nan if q["estimate"] is None or (col==2 and r["heat_cospectrum_status"].startswith("INDETERMINATE")) else q["estimate"] for q,r in zip(tails,rows)]
                axes[level,col].plot(hours,value,color=COLORS[name],label=LABELS[name])
                valid=np.isfinite(value)
                axes[level,col].fill_between(hours,[q["lower"] for q in tails],[q["upper"] for q in tails],where=valid,color=COLORS[name],alpha=.10)
        for ax in axes[level]:ax.set_xlabel("Earth hours");ax.grid(alpha=.2);ax.set_xticks([0,2,4,6])
        axes[level,0].set(title=f"z={h/1000:.3f} km",ylabel="peak wavelength [km]")
        axes[level,1].set(title="w variance at λ < 8 cells",ylabel="fraction",ylim=(-.02,1.02))
        axes[level,2].set(title="Absolute heat throughput at λ < 8 cells",ylabel="fraction",ylim=(-.02,1.02))
    axes[0,0].legend(fontsize=7)
    fig.suptitle("Hourly seed diagnostics · zonal mean removed\nEach grid has its own physical cutoff; gaps are indeterminate. Full4/8/16-cell radial/component tables are in CSV.",fontsize=11)
    fig.savefig(out/"five_height_scale_tails.png",dpi=110);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(14,5),layout="constrained")
    for case in cases:
        for ax,field,label in zip(axes,("w2","heat","n2"),("vertical velocity variance [m²/s²]","resolved heat flux [W/m²]","horizontal-mean N² [s⁻²]")):
            avg=R.average([p[field] for p in case["profiles"]],case["times"],*WINDOW)
            ax.plot(avg,case["profiles"][0]["z"]/1000,color=COLORS[case["name"]],label=LABELS[case["name"]])
            ax.set(xlabel=label,ylabel="height above bottom [km]");ax.grid(alpha=.2)
    for ax in axes:
        for h in heights:ax.axhline(h/1000,color="gray",alpha=.25,linewidth=.6)
    axes[0].legend(fontsize=7);fig.suptitle("0–6 h profiles · five spectral planes marked · no vertical-convergence claim",fontsize=11)
    fig.savefig(out/"full_depth_profiles.png",dpi=115);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout="constrained")
    fields=("mean_ke_J_m2","eddy_to_mean_ke_transfer_W_m2","max_abs_mass_streamfunction_kg_m_s")
    labels=("zonal-mean KE [J/m²]","centered eddy→mean KE transfer [W/m²]","maximum |Ψ| [kg/(m s)]")
    for case in native:
        for ax,key,label in zip(axes,fields,labels):
            ax.plot(case["times"]/3600,[r[key] for r in case["rows"]],color=COLORS[case["name"]],label=LABELS[case["name"]])
            ax.set(xlabel="Earth hours",ylabel=label);ax.grid(alpha=.2);ax.set_xticks([0,2,4,6])
    axes[0].legend(fontsize=7);axes[1].axhline(0,color="black",linewidth=.5)
    fig.suptitle("Startup mean-flow diagnostics · centered transfer is approximate; no steady-cell or cascade claim",fontsize=11)
    fig.savefig(out/"native_mean_flow_comparison.png",dpi=115);plt.close(fig)


def comparisons(native):
    by_name={c["name"]:c for c in native};result={}
    for left,right,label,scope in PAIRS:
        a,b=by_name[left],by_name[right]
        require(np.array_equal(a["times"],b["times"]),"Comparison requires identical saved times")
        metrics={}
        for key in ("mean_ke_J_m2","eddy_to_mean_ke_transfer_W_m2","max_abs_mass_streamfunction_kg_m_s","density_weighted_zonal_mean_v_rms_m_s","density_weighted_zonal_mean_w_rms_m_s"):
            delta=np.asarray([r[key] for r in b["rows"]])-np.asarray([r[key] for r in a["rows"]])
            metrics[key]={"mean_difference":float(P.time_mean(delta,a["times"],0)[0]),"final_difference":float(delta[-1])}
        result[label]={"difference":f"{right} minus {left}","window_seconds":list(WINDOW),"scope":scope,"metrics":metrics}
    return result


def publication_sources():
    here=Path(__file__).parent
    paths=[here/"report.py",here/"guards.py",Path("slurm/ferrel2026-seedtest-publication-v1.sbatch")]
    hashes={str(p):sha(p) for p in paths}
    for path,digest in PINS.values():
        actual=sha(path)
        require(actual==digest,f"Frozen helper changed during reporting: {path}")
        hashes[path]=actual
    return hashes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-root",type=Path,default=Path("data/ferrel2026-seedtest-gpu-v1"))
    parser.add_argument("--gate-root",type=Path,default=Path("data/ferrel2026-seedtest-preflight-v1"))
    parser.add_argument("--prior-root",type=Path,default=Path("data/ferrel2026-preflight-gpu-v1"))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--self-test-only",action="store_true")
    args=parser.parse_args()
    require(not args.output.exists(),"Use a fresh output directory; partial reports are never overwritten")
    args.output.mkdir(parents=True)
    sources_before=publication_sources();validation=analytic_checks()
    (args.output/"analytical_checks.json").write_text(json.dumps(validation,indent=2)+"\n")
    if args.self_test_only:
        (args.output/"self_test.json").write_text(json.dumps({"status":"PASS_ANALYTICAL_READER_TESTS_ONLY","slurm_job_id":os.environ["SLURM_JOB_ID"],"publication_sources":sources_before},indent=2)+"\n")
        print("ANALYTICAL_CHECKS_PASSED_NO_SIMULATION_RESULTS",flush=True);return
    folders,provenance=validate(args.production_root,args.gate_root,args.prior_root)
    native=[P.collect(name,folders[name],args.output) for name in CASES]
    observed=[observe(name,folders[name],args.output) for name in CASES]
    limits=np.maximum(np.max([c["limits"] for c in native],axis=0),1e-12)
    crange=[min(c["tracer_range"][0] for c in native),max(c["tracer_range"][1] for c in native)]
    media=[P.movie([case],args.output/f'{case["name"]}_vorticity_tracer.mp4',limits,crange,1) for case in native]
    native_by_name={c["name"]:c for c in native}
    for left,right,label,scope in PAIRS:
        movie=P.movie([native_by_name[left],native_by_name[right]],args.output/f"{label}.mp4",limits,crange,1)
        movie["comparison_scope"]=scope
        (args.output/f"{label}.json").write_text(json.dumps(movie,indent=2,allow_nan=False)+"\n")
        media.append(movie)
    figures(observed,native,args.output)
    # Recheck provenance and scheduler completion before declaring an immutable
    # report complete; this is read-only and never submits or cancels a job.
    _,after=validate(args.production_root,args.gate_root,args.prior_root)
    require(after==provenance and publication_sources()==sources_before,"Inputs or observer sources changed during reporting")
    report={"status":"COMPLETE_SEED_SCREEN_NOT_CONVERGENCE","compute_location":"Unity allocated CPU job","slurm_job_id":os.environ["SLURM_JOB_ID"],
        "generated_at_utc":datetime.now(timezone.utc).isoformat(),"raw_fields_downloaded":False,"publication_sources":sources_before,"provenance":provenance,
        "analytical_checks":validation,"window_seconds":list(WINDOW),"cases":{c["name"]:c["summary"] for c in observed},
        "native_diagnostics":{c["name"]:c["summary"] for c in native},"comparisons":comparisons(native),"media":media,"caveats":CAVEATS,
        "fourier_convention":"Frozen full complex FFT2 on x-periodic/even-wall-reflected y; no mode filtering; wavelength=2pi/k separately from1/k.",
        "precision_convention":"Float64 saved centered-w and theta anomaly, with direct stored-profile identity checks. Spectral sensitivity includes half-ULP absolute Float64 theta and anomaly/w representation plus the frozen arithmetic-model bounds; not accumulated solver-error bounds.",
        "arithmetic_assumption":"The reused FFT allowance assumes a stable Float64 FFT and uses gamma_(32*N_extended); it is not certified interval arithmetic for arbitrary FFT implementations or twiddle generation.",
        "normalization":{"vorticity":"Fixed asinh ζ/f0, width0.1, full per-depth extrema shared across all4 cases and7frames","tracer":"Fixed asinh C, width0.01, common full range across all4 cases, depths and times"},
        "reference_policy":"Old cellwise-noise pilot results are separately published references; no old fields are relabeled as five-plane data or treated as missing controls.",
        "software":{"python":sys.version,"numpy":np.__version__,"h5py":h5py.__version__,"matplotlib":R.matplotlib.__version__}}
    (args.output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    lines=["# Six-hour horizontal thermal-seed screen","","Complete0–6h startup diagnostics; not convergence or developed convection.",
           "The four cases share fixed movie scales. Each movie has three actual depth columns, vorticity above passive tracer.",
           "Five Float64 diagnostic planes are5.625,9.375,19.375,59.375,99.375km above the bottom; the three video planes remain19.375,59.375,99.375km.",
           "No unheated modal control exists. Modal-vs-zero changes are seed sensitivity; fine-vs-coarse changes are horizontal-grid sensitivity.","",
           "Watch on the [GitHub gallery landing page](https://npizzo.github.io/jupiter-ffr-breeze/).","", "## Movies",""]
    lines += [f'- [{item["file"]}]({item["file"]})' for item in media]
    lines += ["","## Interpretation limits",""]+[f"- {line}" for line in CAVEATS]
    (args.output/"README.md").write_text("\n".join(lines)+"\n")
    assets=[{"file":p.name,"bytes":p.stat().st_size,"sha256":sha(p)} for p in sorted(args.output.iterdir()) if p.is_file()]
    require(all(a["bytes"]<24_000_000 for a in assets),"Publication asset exceeds browser fallback limit")
    (args.output/"manifest.json").write_text(json.dumps({"status":"READY_FOR_REVIEW","assets":assets},indent=2)+"\n")
    print("READY_FOR_REVIEW",args.output,flush=True)


if __name__=="__main__":main()
