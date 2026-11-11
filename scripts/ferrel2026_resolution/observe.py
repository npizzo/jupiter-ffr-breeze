"""Allocated-only horizontal-resolution and heat-transport screen for Ferrel runs.

No model state is changed. All modes are retained. This screen cannot establish
convergence, resolve missing vertical scales, or separate waves from convection.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tomllib

# Fail before importing scientific libraries when invoked outside an allocation.
if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Run this observer inside an allocated Unity job")

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = np.finfo(float).eps
TINY = np.finfo(float).tiny
GRID_CACHE = {}
LABELS = {"b0":"Unheated", "b162":"Heated 162 W/m²", "boundary162":"Boundary sensitivity",
          "axis0":"Axisymmetric unheated", "axis162":"Axisymmetric heated", "fine162":"1024² heated"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gamma(operations):
    """Standard worst-case summation/product roundoff factor gamma_n."""
    n_eps = float(operations)*EPS
    if n_eps >= .01:
        raise ValueError("Arithmetic envelope inappropriate for requested size")
    return n_eps/(1-n_eps)


def rms(a):
    return float(np.sqrt(np.mean(np.asarray(a, dtype=float)**2)))


def finite(d, key, shape=None):
    a=np.asarray(d[key][()])
    if shape is not None and a.shape!=shape:
        raise ValueError(f"{key}: shape {a.shape} != {shape}")
    if not np.isfinite(a).all():
        raise ValueError(f"Nonfinite {key}")
    return a


def half_ulp32(a):
    if a.dtype!=np.dtype("float32"):
        raise ValueError("Expected saved Float32 plane")
    center=a.astype(float)
    up=np.nextafter(a,np.float32(np.inf)).astype(float)
    down=np.nextafter(a,np.float32(-np.inf)).astype(float)
    q=.5*np.maximum(up-center,center-down)
    if not np.isfinite(q).all():
        raise ValueError("Cannot bound quantization at Float32 range limit")
    return q


def projected(a, mode):
    axis=None if mode=="plane" else 1
    return a-np.mean(a,axis=axis,keepdims=True)


def geometry(nx,ny,lx,ly):
    key=(nx,ny,lx,ly)
    if key not in GRID_CACHE:
        kx=2*np.pi*np.fft.fftfreq(nx,lx/nx)[None,:]
        ky=2*np.pi*np.fft.fftfreq(2*ny,ly/ny)[:,None]
        radius=np.sqrt(kx*kx+ky*ky)
        dk=min(2*np.pi/lx,np.pi/ly)
        bins=np.floor(radius/dk+1e-10).astype(int)
        k=(np.arange(int(bins.max())+1)+.5)*dk
        masks={}
        for cells in (4,8,16):
            # Radial wavelength and direction-specific Nyquist fractions are
            # separate: diagonal waves need not have the same component scale.
            masks[f"radial_lt_{cells}cells"]=radius>2*np.pi/(cells*max(lx/nx,ly/ny))
            masks[f"component_lt_{cells}cells"]=(np.abs(kx)>2*np.pi/(cells*lx/nx)) | (np.abs(ky)>2*np.pi/(cells*ly/ny))
        GRID_CACHE[key]=(kx,ky,radius,dk,bins,k,masks)
    return GRID_CACHE[key]


def fraction_interval(numerator,total,error):
    if total<=error or total<=0:
        return {"estimate":None,"lower":0.,"upper":1.,"status":"INDETERMINATE_DENOMINATOR"}
    return {"estimate":float(numerator/total),
            "lower":float(max(0,numerator-error)/(total+error)),
            "upper":float(min(1,(numerator+error)/(total-error))),"status":"BOUNDED"}


def spectrum(w,t,qw,qt,m,mode):
    ny,nx=w.shape
    kx,ky,rad,dk,bins,k,masks=geometry(nx,ny,m["lx"],m["ly"])
    w0=projected(w,mode);t0=projected(t,mode)
    # Unweighted demeaning is an orthogonal projection, hence its L2 operator
    # norm is one: raw rms(q) bounds the projected quantization error.
    ew,et=rms(qw),rms(qt)
    # Absolute theta mean removal is cancellation-prone. Add a separate
    # conservative Float64 reduction/subtraction error, never hide it in q.
    reduction_n=w.size if mode=="plane" else nx
    arithmetic_w=gamma(reduction_n+4)*float(np.max(np.abs(w)))
    arithmetic_t=gamma(reduction_n+4)*float(np.max(np.abs(t)))
    extend=lambda a:np.concatenate((a,a[::-1,:]),axis=0)
    W=np.fft.fft2(extend(w0),norm="forward")
    T=np.fft.fft2(extend(t0),norm="forward")
    power=np.abs(W)**2
    cospectrum=(W.conj()*T).real
    # The zero mode is numerically retained in Parseval checks and explicitly
    # excluded from scale measures after those checks.
    sw,st=rms(w0),rms(t0)
    parseval_error=abs(float(power.sum())-sw*sw)
    cospectral_error=abs(float(cospectrum.sum())-float(np.mean(w0*t0)))
    # Conservative arithmetic-model allowance assuming a numerically stable
    # Float64 FFT: gamma_(32 N) exceeds the DFT/product/reduction operation depth.
    # This is not certified interval arithmetic for arbitrary FFT libraries or
    # twiddle generation. Keep it distinct from rigorous input half-ULP bounds.
    fft_gamma=gamma(32*W.size)
    parseval_bound=fft_gamma*max(sw*sw,TINY)
    cospectral_bound=fft_gamma*max(sw*st,TINY)
    if parseval_error>parseval_bound or cospectral_error>cospectral_bound:
        raise ValueError(f"FFT identity failure: {parseval_error}/{parseval_bound}, {cospectral_error}/{cospectral_bound}")
    power_quant=2*sw*ew+ew*ew
    heat_quant=ew*st+et*sw+ew*et
    ew_total,et_total=ew+arithmetic_w,et+arithmetic_t
    power_arithmetic=(2*sw*ew_total+ew_total**2)-power_quant+parseval_bound
    heat_arithmetic=(ew_total*st+et_total*sw+ew_total*et_total)-heat_quant+cospectral_bound
    power_bound=power_quant+power_arithmetic
    heat_bound=heat_quant+heat_arithmetic
    nonzero=rad>0
    power=np.where(nonzero,power,0)
    cospectrum=np.where(nonzero,cospectrum,0)
    shell=lambda a:np.bincount(bins.ravel(),weights=a.ravel(),minlength=len(k))
    P=shell(power);C=shell(cospectrum);A=shell(np.abs(cospectrum))
    P[0]=C[0]=A[0]=0
    total=float(power.sum());abs_heat=float(np.abs(cospectrum).sum())
    peak=int(np.argmax(P[1:])+1) if total>power_bound else None
    peak_k=None if peak is None else float(shell(rad*power)[peak]/P[peak])
    candidate=np.flatnonzero(P>=max(float(P.max())-2*power_bound,0))
    candidate=candidate[candidate>0]
    heat_valid=st>et+arithmetic_t and abs_heat>heat_bound and sw>ew+arithmetic_w
    record={"demeaning":mode,"w_variance_m2_s2":total,"theta_variance_K2":st*st,
            "theta_quantization_error_rms_K":et,"theta_arithmetic_error_rms_K":arithmetic_t,
            "w_quantization_error_rms_m_s":ew,"w_arithmetic_error_rms_m_s":arithmetic_w,
            "w_power_l1_quantization_bound_m2_s2":power_quant,
            "w_power_l1_arithmetic_bound_m2_s2":power_arithmetic,
            "heat_cospectrum_sum_K_m_s":float(cospectrum.sum()),
            "heat_cospectrum_positive_K_m_s":float(cospectrum[cospectrum>0].sum()),
            "heat_cospectrum_negative_K_m_s":float(cospectrum[cospectrum<0].sum()),
            "heat_absolute_mode_throughput_K_m_s":abs_heat,
            "heat_cospectrum_l1_quantization_bound_K_m_s":heat_quant,
            "heat_cospectrum_l1_arithmetic_bound_K_m_s":heat_arithmetic,
            "heat_cospectrum_status":"BOUNDED_NOT_MODEWISE_CONVERGED" if heat_valid else "INDETERMINATE_PRECISION_FLOOR",
            "parseval_absolute_residual_m2_s2":parseval_error,"parseval_arithmetic_bound_m2_s2":parseval_bound,
            "cospectral_identity_residual_K_m_s":cospectral_error,"cospectral_identity_arithmetic_bound_K_m_s":cospectral_bound,
            "peak_shell_k_rad_m":peak_k,
            "peak_shell_wavelength_m":None if peak is None else float(2*np.pi/peak_k),
            "peak_shell_inverse_wavenumber_m":None if peak is None else float(1/peak_k),
            "peak_shell_cells_per_wavelength":None if peak is None else float(2*np.pi/peak_k/max(m["lx"]/nx,m["ly"]/ny)),
            "peak_shell_k_lower_rad_m":None if peak is None else float(peak*dk),
            "peak_shell_k_upper_rad_m":None if peak is None else float((peak+1)*dk),
            "peak_uncertainty_candidate_shells":int(len(candidate)),
            "peak_precision_status":"NO_RESOLVED_POWER" if peak is None else ("MULTIPLE_SHELLS_WITHIN_GLOBAL_ERROR_BUDGET" if len(candidate)>1 else "ONE_SHELL_WITHIN_ARITHMETIC_MODEL_BUDGET"),
            "tail_fractions":{}}
    for key,mask in masks.items():
        record["tail_fractions"][key]={
            "w_variance":fraction_interval(float(power[mask].sum()),total,power_bound),
            "absolute_heat_cospectrum":fraction_interval(float(np.abs(cospectrum[mask]).sum()),abs_heat,heat_bound),
            "signed_heat_cospectrum_K_m_s":float(cospectrum[mask].sum())}
    # Absolute-cospectral errors apply to any shell/subset; this is conservative
    # and does not claim independent errors for different bins.
    return record,{"k":k,"dk":dk,"w_power":P,"heat_signed":C,"heat_absolute":A,
                   "power_error":power_bound,"heat_error":heat_bound}


def profile_check(w,t,qw,qt,theta_mean,rho,cp,exner,vertical_ke,heat_flux):
    raw_w2=float(np.mean(w*w));source_w2=float(2*vertical_ke/rho)
    var_quant=float(np.mean(2*np.abs(w)*qw+qw*qw))
    # Model and observer each reduce N values; enlarge true throughput by the
    # quantization contribution before covering both reductions/products.
    var_round=gamma(2*w.size+32)*max(raw_w2+var_quant,abs(source_w2),TINY)
    anomaly=t-theta_mean
    observed_cov=float(np.mean(w*anomaly));source_cov=float(heat_flux/(rho*cp*exner))
    cov_quant=float(np.mean(qw*np.abs(anomaly)+qt*np.abs(w)+qw*qt))
    # theta_mean is the saved Float64 value used by the model itself. This
    # avoids assuming NumPy and Julia computed identical cancellation-prone means.
    cov_round=gamma(2*w.size+64)*max(float(np.mean(np.abs(w*anomaly)))+cov_quant,abs(source_cov),TINY)
    # Subtracting theta_mean also incurs a Float64 subtraction error, even with
    # the exact Float64 theta plane. Bound its absolute input scales explicitly.
    cov_round+=2*EPS*float(np.mean((np.abs(w)+qw)*(np.abs(t)+qt+abs(theta_mean))))
    var_error=abs(raw_w2-source_w2);cov_error=abs(observed_cov-source_cov)
    if var_error>var_quant+var_round:
        raise ValueError("Saved w plane disagrees with Float64 vertical KE beyond quantization/arithmetic bounds")
    if cov_error>cov_quant+cov_round:
        raise ValueError("Saved w/theta plane disagrees with Float64 heat flux beyond quantization/arithmetic bounds")
    return {"raw_plane_w2_m2_s2":raw_w2,"profile_w2_m2_s2":source_w2,
            "w2_profile_residual_m2_s2":var_error,"w2_profile_quantization_bound_m2_s2":var_quant,
            "w2_profile_arithmetic_bound_m2_s2":var_round,
            "plane_covariance_using_saved_theta_mean_K_m_s":observed_cov,
            "profile_covariance_K_m_s":source_cov,"heat_profile_residual_K_m_s":cov_error,
            "heat_profile_quantization_bound_K_m_s":cov_quant,"heat_profile_arithmetic_bound_K_m_s":cov_round,
            "status":"WITHIN_DERIVED_BOUNDS_NOT_CONVERGENCE"}


def self_checks():
    nx=ny=64;lx=2*np.pi;ly=np.pi
    x=(np.arange(nx)+.5)*lx/nx;y=(np.arange(ny)+.5)*ly/ny
    X,Y=np.meshgrid(x,y)
    w=np.cos(2*X)*np.cos(Y)
    m={"lx":lx,"ly":ly};zero=np.zeros_like(w)
    known_k=np.sqrt(5.);known_variance=.25
    for mode in ("plane","zonal"):
        positive,_=spectrum(w,170+.25*w,zero,zero,m,mode)
        negative,_=spectrum(w,170-.25*w,zero,zero,m,mode)
        tolerance=gamma(32*2*nx*ny)
        if abs(positive["w_variance_m2_s2"]-known_variance)>tolerance*known_variance:
            raise ValueError("Analytical cosine variance/normalization test failed")
        if abs(positive["peak_shell_wavelength_m"]-2*np.pi/known_k)>tolerance*2*np.pi/known_k:
            raise ValueError("Analytical physical wavelength test failed")
        if abs(positive["peak_shell_inverse_wavenumber_m"]-1/known_k)>tolerance/known_k:
            raise ValueError("Inverse wavenumber versus wavelength test failed")
        for row,sign in ((positive,1),(negative,-1)):
            bound=row["heat_cospectrum_l1_arithmetic_bound_K_m_s"]
            if abs(row["heat_cospectrum_sum_K_m_s"]-sign*.25*known_variance)>bound:
                raise ValueError("Analytical covariance sign/normalization test failed")
    w32=w.astype(np.float32);qw=half_ulp32(w32)
    theta_true=170+1e-7*w
    theta32=theta_true.astype(np.float32);qt=half_ulp32(theta32)
    weak,_=spectrum(w32.astype(float),theta32.astype(float),qw,qt,m,"plane")
    if weak["heat_cospectrum_status"]!="INDETERMINATE_PRECISION_FLOOR":
        raise ValueError("Quantized weak-temperature case must be indeterminate")
    heat_true=float(np.mean(w*(theta_true-theta_true.mean())))
    profile_check(w32.astype(float),theta32.astype(float),qw,qt,float(theta_true.mean()),1.,1.,1.,known_variance/2,heat_true)
    try:
        profile_check(w32.astype(float),theta32.astype(float),qw,qt,float(theta_true.mean()),1.,1.,1.,known_variance/2,1e-2)
    except ValueError:
        pass
    else:
        raise ValueError("Quantization bound incorrectly accepted an inconsistent heat profile")
    return {"status":"PASS","checks":["even-wall single Fourier mode variance and physical wavelength",
        "plane/zonal projection and positive/negative covariance signs",
        "weak absolute Float32 theta correctly labeled indeterminate",
        "derived precision envelope accepts valid weak profile and rejects inconsistent heat flux"]}


def average(values,times,lo,hi):
    take=(times>=lo-1e-6)&(times<=hi+1e-6)
    selected=times[take]
    if len(selected)<2 or not np.isclose(selected[0],lo,rtol=0,atol=1e-6) or not np.isclose(selected[-1],hi,rtol=0,atol=1e-6):
        raise ValueError("Averaging window is not covered by saved frames")
    trap=np.trapezoid if hasattr(np,"trapezoid") else np.trapz
    return trap(np.asarray(values)[take],selected,axis=0)/(hi-lo)


def windows(end):
    choices=[("early_0_6h",0.,21600.)]
    if end>=172800:choices.append(("startup_excluded_12_48h",43200.,172800.))
    if end!=21600 and end!=172800:choices.append(("entire_complete_run",0.,float(end)))
    return choices


def nominal(m,rho,theta,pressure,z):
    h=m["cap_height"]
    bottom_exner=(m["base_pressure"]/1e5)**((8.314462618/m["molar_mass"])/m["heat_capacity"])
    # rho at the lowest cell center is an explicitly declared approximation to
    # bottom density; this is a nominal dimensionless analogue, not fitted physics.
    bottom_temperature=m["theta0"]*bottom_exner
    B=m["gravity"]*m["heat_flux"]/(rho[0]*m["heat_capacity"]*bottom_temperature)
    kappa=m["thermal_diffusivity"];nu=m["viscosity"];f=abs(m["f"])
    idx=int(np.argmin(abs(z-h)))
    return {"active_depth_assumption_m":h,"active_depth_basis":"configured neutral-layer cap height; not a measured mixing depth",
            "lower_Pr":nu/kappa,"cap_nominal_Pr":nu/m["cap_diffusivity"],
            "local_Ek_nu_over_f_H2":nu/(f*h*h),
            "bottom_buoyancy_flux_analogue_m2_s3":B,
            "flux_Rayleigh_analogue":B*h**4/(nu*kappa*kappa),
            "flux_velocity_analogue_m_s":float(np.cbrt(max(B,0)*h)),
            "flux_Rossby_analogue":float(np.cbrt(max(B,0)*h)/(f*h)),
            "reference_density_log_contrast_to_cap":float(np.log(rho[0]/rho[idx])),
            "density_reference_height_m":float(z[0]),"cap_density_height_m":float(z[idx]),
            "scope":"Explicit-coefficient, local-f nominal analogues only; tapered diffusivity, anelasticity, boundaries and numerical diffusion prevent direct RBC criticality claims"}


def gate(path):
    value=tomllib.loads(path.read_text())
    if value["status"]!="PASS":raise ValueError(f"Gate not passed: {path}")
    for name,expected in value["source_hashes"].items():
        if sha(Path("scripts/ferrel2026_v3")/name)!=expected:raise ValueError("Scientific source fingerprint changed")
    for name,expected in value["configuration_hashes"].items():
        if sha(Path("configs/ferrel2026")/f"{name}.toml")!=expected:raise ValueError("Configuration fingerprint changed")
    return {"file":str(path),"sha256":sha(path),"status":"PASS"}


def flatten(row):
    out={}
    for key,value in row.items():
        if key=="tail_fractions":
            for band,items in value.items():
                for item,content in items.items():
                    if isinstance(content,dict):
                        for sub,v in content.items():out[f"{band}_{item}_{sub}"]=v
                    else:out[f"{band}_{item}"]=content
        else:out[key]=value
    return out


def write_csv(path,rows):
    if not rows:raise ValueError("Refusing empty diagnostic table")
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def observe(name,folder,out):
    m=tomllib.loads((folder/"metadata.toml").read_text())
    if m["status"]!="COMPLETE":raise ValueError(f"Incomplete case: {folder}")
    if m["meridional_boundary"]!="channel" or m["coriolis_model"]!="beta-plane":raise ValueError("Wrong geometry")
    if not np.isfinite(m["f"]) or m["f"]==0:raise ValueError("Invalid f")
    nx,ny,nz=(m[key] for key in ("nx","ny","nz"))
    frames=sorted(folder.glob("frame_*.jld2"))
    times=[];records=[];checks=[];activities=[];all_spectra={};profile_series=[];regime=None
    for path in frames:
        with h5py.File(path) as d:
            t=float(finite(d,"t"));times.append(t)
            z=finite(d,"z",(nz,)).astype(float);rho=finite(d,"density",(nz,)).astype(float)
            pressure=finite(d,"pressure",(nz,)).astype(float)
            theta_profile=finite(d,"theta_profile",(nz,)).astype(float)
            ke=finite(d,"vertical_ke_profile",(nz,)).astype(float)
            heat=finite(d,"heat_flux_profile",(nz,)).astype(float)
            diffusion=finite(d,"diffusive_heat_flux_profile",(nz,)).astype(float)
            n2=finite(d,"mean_N2_profile",(nz,)).astype(float)
            heights=finite(d,"video_heights",(3,)).astype(float)
            indices=np.asarray([np.argmin(abs(z-h)) for h in m["video_heights_requested_m"]])
            if len(set(indices))!=3 or not np.array_equal(heights,z[indices]) or not np.array_equal(heights,m["video_heights_actual_m"]):
                raise ValueError("Diagnostic/video indices or heights do not match")
            W=finite(d,"wcenter_slices",(3,ny,nx));T=finite(d,"theta_slices",(3,ny,nx))
            first_theta=finite(d,"theta",(ny,nx))
            if first_theta.dtype!=np.dtype("float64") or not np.array_equal(first_theta.astype(np.float32),T[0]):
                raise ValueError("Cannot verify that Float64 theta matches the first video plane")
            exner=(pressure/1e5)**((8.314462618/m["molar_mass"])/m["heat_capacity"])
            if regime is None:regime=nominal(m,rho,theta_profile,pressure,z)
            profile_series.append({"z":z,"w2":2*ke/rho,"heat":heat,"diffusion":diffusion,"n2":n2,"theta":theta_profile})
            w2=2*ke/rho
            active=int(np.argmax(w2));thermal=int(np.argmax(abs(heat)))
            activities.append({"time_seconds":t,"earth_day":t/86400,
                "height_of_max_w_variance_m":float(z[active]),"max_w_rms_m_s":float(np.sqrt(max(w2[active],0))),
                "height_of_max_abs_heat_flux_m":float(z[thermal]),"max_abs_resolved_heat_flux_W_m2":float(abs(heat[thermal])),
                "max_positive_mean_N2_s2":float(max(n2.max(),0)),"min_mean_N2_s2":float(n2.min()),
                "N_profile_max_times_configured_dt":float(np.sqrt(max(n2.max(),0))*m["max_dt"]),
                "configured_dt_seconds":m["max_dt"],"wave_measure_scope":"horizontal-mean N only; not local Nmax or a stability guarantee"})
            for level,index in enumerate(indices):
                w=W[level].astype(float);qw=half_ulp32(W[level])
                theta=first_theta if level==0 else T[level].astype(float)
                qt=np.zeros_like(theta) if level==0 else half_ulp32(T[level])
                check=profile_check(w,theta,qw,qt,theta_profile[index],rho[index],m["heat_capacity"],exner[index],ke[index],heat[index])
                check.update(time_seconds=t,level=level,height_m=float(z[index]),
                             theta_storage="verified matching Float64 theta" if level==0 else "absolute Float32 theta_slices")
                checks.append(check)
                for mode in ("plane","zonal"):
                    row,shells=spectrum(w,theta,qw,qt,m,mode)
                    row.update(time_seconds=t,earth_day=t/86400,level=level,height_m=float(z[index]),theta_storage=check["theta_storage"])
                    records.append(row);all_spectra.setdefault((level,mode),[]).append(shells)
        print("SCREENED",name,f"t={t/3600:g}h",flush=True)
    times=np.asarray(times)
    expected=np.arange(0,m["stop_seconds"]+m["output_seconds"]/2,m["output_seconds"])
    if times.shape!=expected.shape or not np.allclose(times,expected,rtol=0,atol=1e-6):
        raise ValueError("Need complete chronological initial-to-final frames")
    rows=[];window_summaries={};plot_spectra={};profile_rows=[]
    for label,lo,hi in windows(m["stop_seconds"]):
        window_summaries[label]={"start_seconds":lo,"end_seconds":hi,"planes":{}}
        for (level,mode),seq in all_spectra.items():
            k=seq[0]["k"];dk=seq[0]["dk"]
            avg={key:average([s[key] for s in seq],times,lo,hi) for key in ("w_power","heat_signed","heat_absolute","power_error","heat_error")}
            plot_spectra[(label,level,mode)]=(k,dk,avg)
            P=avg["w_power"];A=avg["heat_absolute"];C=avg["heat_signed"]
            pb=float(avg["power_error"]);hb=float(avg["heat_error"])
            per={"w_variance_m2_s2":float(P.sum()),"w_power_l1_error_bound_m2_s2":pb,
                 "signed_heat_covariance_K_m_s":float(C.sum()),"absolute_heat_mode_throughput_K_m_s":float(A.sum()),
                 "heat_cospectrum_l1_error_bound_K_m_s":hb,
                 "heat_throughput_status":"INDETERMINATE_PRECISION_FLOOR" if A.sum()<=hb else "BOUNDED_NOT_MODEWISE_CONVERGED"}
            window_summaries[label]["planes"][f"{level}_{mode}"]=per
            for j in range(1,len(k)):
                rows.append({"window":label,"level":level,"height_m":float(heights[level]),"demeaning":mode,
                    "k_lower_rad_m":j*dk,"k_upper_rad_m":(j+1)*dk,"k_center_rad_m":float(k[j]),
                    "wavelength_m":float(2*np.pi/k[j]),"inverse_wavenumber_m":float(1/k[j]),
                    "w_shell_variance_m2_s2":float(P[j]),"signed_heat_shell_K_m_s":float(C[j]),
                    "absolute_mode_heat_in_shell_K_m_s":float(A[j]),
                    "any_subset_w_power_error_bound_m2_s2":pb,"any_subset_heat_error_bound_K_m_s":hb})
        for key in ("w2","heat","diffusion","n2","theta"):
            average_profile=average([item[key] for item in profile_series],times,lo,hi)
            for i,value in enumerate(average_profile):
                profile_rows.append({"window":label,"height_m":float(z[i]),"quantity":key,"value":float(value)})
    write_csv(out/f"{name}_scale_timeseries.csv",[flatten(r) for r in records])
    write_csv(out/f"{name}_profile_agreement.csv",checks)
    write_csv(out/f"{name}_vertical_activity.csv",activities)
    write_csv(out/f"{name}_window_spectra.csv",rows)
    write_csv(out/f"{name}_window_profiles.csv",profile_rows)
    plot(name,m,plot_spectra,profile_series,times,out)
    dt=np.genfromtxt(folder/"native_budget/checks.csv",delimiter=",",names=True)
    dt=np.atleast_1d(dt)
    post=dt[dt["start_time"]>=43200] if times[-1]>43200 else dt
    timestep={"scope":"sampled actual SSPRK3 clock increments, not continuous limiter telemetry",
              "sample_count":int(len(post)),"configured_cap_seconds":m["max_dt"],
              "minimum_sampled_dt_seconds":float(post["dt"].min()) if len(post) else None,
              "median_sampled_dt_seconds":float(np.median(post["dt"])) if len(post) else None,
              "maximum_sampled_dt_seconds":float(post["dt"].max()) if len(post) else None,
              "native_diffusive_limit_configured":False,"native_internal_wave_limit_configured":False,
              "warning":"Current wizard controls advection and max_dt; no inference of the active limiter from these samples alone."}
    return {"case":name,"label":LABELS[name],"source":str(folder),"grid":[nx,ny,nz],
            "dx_m":m["lx"]/nx,"dy_m":m["ly"]/ny,"dz_m":m["height"]/nz,
            "metadata_sha256":sha(folder/"metadata.toml"),"driver_sha256":m["script_sha256"],
            "configuration_sha256":m["configuration_sha256"],"frames":len(frames),"end_seconds":float(times[-1]),
            "heights_m":heights.tolist(),"nominal_regime":regime,"timestep":timestep,"windows":window_summaries,
            "heat_cospectrum_precision_indeterminate_records":sum(r["heat_cospectrum_status"].startswith("INDETERMINATE") for r in records),
            "total_plane_demeaning_records":len(records),"profile_checks":"WITHIN_DERIVED_QUANTIZATION_PLUS_ARITHMETIC_BOUNDS"}


def plot(name,m,spectra,profiles,times,out):
    fig,axes=plt.subplots(3,2,figsize=(12,11),layout="constrained")
    choices=windows(m["stop_seconds"]);colors={"plane":"tab:blue","zonal":"tab:orange"}
    for level in range(3):
        for wi,(label,lo,hi) in enumerate(choices):
            for mode in ("plane","zonal"):
                k,dk,a=spectra[(label,level,mode)];pick=k>0
                style="-" if wi==len(choices)-1 else "--"
                legend=f"{mode}-mean removed · {lo/3600:g}–{hi/3600:g}h"
                axes[level,0].loglog(2*np.pi/k[pick]/1000,np.maximum(a["w_power"][pick]*k[pick]/dk,TINY),style,color=colors[mode],label=legend)
                heat=a["heat_signed"][pick]*k[pick]/dk
                axes[level,1].semilogx(2*np.pi/k[pick]/1000,heat,style,color=colors[mode],label=legend)
        for ax in axes[level]:
            for cells in (4,8,16):ax.axvline(cells*m["lx"]/m["nx"]/1000,color="gray",alpha=.25,linewidth=.8)
            ax.grid(alpha=.2);ax.set_xlabel("horizontal wavelength [km]")
        axes[level,0].set_ylabel("w variance per log k [m²/s²]")
        axes[level,1].set_ylabel("signed w–theta covariance per log k [K m/s]")
        axes[level,1].axhline(0,color="black",linewidth=.6)
        axes[level,0].set_title(f'z={m["video_heights_actual_m"][level]/1000:.2f} km · all modes retained')
        axes[level,1].set_title("Consult precision status/bounds before interpreting weak heat flux")
    axes[0,0].legend(fontsize=7);fig.suptitle(f'{LABELS[name]} · horizontal-resolution screen, not convergence')
    fig.savefig(out/f"{name}_horizontal_scales.png",dpi=115);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,5),layout="constrained")
    z=profiles[0]["z"]/1000
    for label,lo,hi in choices:
        for ax,key,xlabel in zip(axes,("w2","heat","n2"),("vertical velocity variance [m²/s²]","resolved heat flux [W/m²]","horizontal-mean N² [s⁻²]")):
            value=average([p[key] for p in profiles],times,lo,hi)
            ax.plot(value,z,label=f"{lo/3600:g}–{hi/3600:g} h");ax.set(xlabel=xlabel,ylabel="height [km]")
            for height in m["video_heights_actual_m"]:ax.axhline(height/1000,color="gray",alpha=.3,linewidth=.7)
            ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle(f'{LABELS[name]} · full-depth aggregate profiles; horizontal spectra exist at three heights only')
    fig.savefig(out/f"{name}_vertical_activity.png",dpi=115);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",choices=("pilot","production"))
    p.add_argument("--production-root",type=Path,default=Path("data/ferrel2026-gpu-v1"))
    p.add_argument("--preflight-root",type=Path,default=Path("data/ferrel2026-preflight-gpu-v1"))
    p.add_argument("--include-fine",action="store_true")
    p.add_argument("--completed-job",action="append",default=[])
    p.add_argument("--self-test-only",action="store_true")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    if args.output.exists() and any(args.output.iterdir()):p.error("Fresh output directory required")
    args.output.mkdir(parents=True,exist_ok=True)
    validation=self_checks()
    (args.output/"analytical-checks.json").write_text(json.dumps(validation,indent=2)+"\n")
    print("ANALYTICAL_SELF_CHECKS_PASS",flush=True)
    if args.self_test_only:return
    if args.stage is None or not args.completed_job:p.error("Full observation requires stage and successful prerequisite jobs")
    if not all(job.isdigit() for job in args.completed_job):p.error("Integer completed-job IDs required")
    status=subprocess.check_output(["sacct","-X","-j",",".join(args.completed_job),"-n","-P","-o","JobIDRaw,State,ExitCode"],text=True)
    jobs={parts[0]:parts[1:] for line in status.splitlines() if len(parts:=line.strip().split("|"))>=3}
    for job in args.completed_job:
        if jobs.get(job,[])[:2]!=["COMPLETED","0:0"]:raise ValueError(f"Prerequisite job not successful: {job}")
    gates=[gate(args.preflight_root/"PASS512.toml")]
    if args.stage=="pilot":
        sources={name:args.preflight_root/f"pilot512_{name}" for name in ("b0","b162","axis0","axis162")}
    else:
        sources={name:args.production_root/f"ferrel2026_{name}" for name in ("b0","b162","boundary162")}
        sources.update({name:args.preflight_root/f"pilot512_{name}" for name in ("axis0","axis162")})
    if args.include_fine:
        gates.append(gate(args.preflight_root/"PASS1024.toml"));sources["fine162"]=args.preflight_root/"pilot1024_b162"
    cases={name:observe(name,folder,args.output) for name,folder in sources.items()}
    report={"status":"COMPLETE_RESOLUTION_SCREEN_NOT_CONVERGENCE","stage":args.stage,
        "compute_location":"Unity allocated CPU job","slurm_job_id":os.environ["SLURM_JOB_ID"],
        "generated_at_utc":datetime.now(timezone.utc).isoformat(),"observer_sha256":sha(__file__),
        "analytical_checks":validation,
        "successful_prerequisite_jobs":{j:jobs[j] for j in args.completed_job},"gates":gates,"cases":cases,
        "fourier_convention":"Full complex fft2 norm=forward on periodic x and even y reflection to2Ly. ky increments pi/Ly. No mode filtering. Radial shells of width min(2pi/Lx,pi/Ly).",
        "precision_convention":"Half-ULP Float32 input bounds, contractive L2 mean removal, Cauchy/Parseval L1 spectral error bounds; Float64 arithmetic allowances reported separately. Matching Float64 theta used at the first actual video plane.",
        "arithmetic_assumption":"Float64 FFT is numerically stable; gamma_(32*N_extended) is a conservative arithmetic-model allowance, not certified interval arithmetic for arbitrary FFT implementations/twiddles. Quantization and arithmetic budgets remain separate.",
        "interpretation_limits":[
            "Three horizontal planes cannot establish all-depth, vertical or bottom-boundary-layer resolution.",
            "1024²×96 is horizontal refinement only. Six-hour comparison covers startup, not developed48h convection.",
            "First12h exclusion is a convention, not a measured equilibration time.",
            "Hourly snapshots cannot resolve the reference cap buoyancy-wave period; no temporal wave/convection classification is attempted.",
            "Positive instantaneous heat cospectrum alone does not identify sustained convection; examine temporal averages, profile transport and matched controls.",
            "Tail fractions and peak scales are screens, not numerical convergence or physical-fidelity proofs. Uncertainty is global/subset bounded, not independent mode errors.",
            "Nominal regime parameters omit numerical diffusion and do not supply a calibrated Boussinesq onset threshold.",
            "No raw model fields are exported or modified; these compact diagnostics remain on Unity until publication review."]}
    (args.output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    assets=[{"file":x.name,"bytes":x.stat().st_size,"sha256":sha(x)} for x in sorted(args.output.iterdir()) if x.is_file()]
    (args.output/"manifest.json").write_text(json.dumps({"status":"READY_FOR_REVIEW","assets":assets},indent=2)+"\n")
    print("READY_FOR_REVIEW",args.output,flush=True)


if __name__=="__main__":main()
