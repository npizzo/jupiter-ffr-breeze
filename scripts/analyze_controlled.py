"""Analyze the fixed upper sink, tapered mixing, and vertical refinement campaign."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import tomllib

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from analyze_convection import cumulative_trapezoid, save_figure, time_mean


CATALOG = [
    ("controlled48_on_p020", "0.2 K: heated + cooled", "#0072B2", "-"),
    ("controlled48_off_p020", "0.2 K: cooled only", "#0072B2", "--"),
    ("controlled48_on_p002", "0.02 K: heated + cooled", "#D55E00", "-"),
    ("controlled48_off_p002", "0.02 K: cooled only", "#D55E00", "--"),
    ("controlled48_unforced", "0.02 K: unforced", "#009E73", ":"),
    ("controlled48_refined", "0.02 K: heated + cooled, Nz=128", "#CC79A7", "-."),
]


def read_case(folder, label, color, style):
    m = tomllib.loads((folder / "metadata.toml").read_text())
    d = np.genfromtxt(folder / "diagnostics.csv", delimiter=",", names=True)
    if m["status"] != "COMPLETE" or d["time_seconds"][-1] != m["stop_seconds"]:
        raise ValueError(f"Incomplete case: {folder}")
    if not all(np.isfinite(d[key]).all() for key in d.dtype.names):
        raise ValueError(f"Nonfinite diagnostic: {folder}")
    files = sorted(folder.glob("frame_*.jld2"))
    with h5py.File(files[-1]) as f:
        fields = {key: f[key][:] for key in ("theta_profile", "heat_flux_profile", "diffusive_heat_flux_profile", "zeta", "w_xz")}
        fields["t"] = float(f["t"][()])
        if "kappa_profile" in f:
            fields["kappa_profile"] = f["kappa_profile"][:]
            fields["upper_energy_source"] = f["upper_energy_source"][:]
    with h5py.File(files[0]) as f:
        initial_theta = f["theta_profile"][:]
        initial_time = float(f["t"][()])
    if initial_time != 0 or fields["t"] != m["stop_seconds"]:
        raise ValueError(f"Missing initial/final frame: {folder}")
    z = (np.arange(m["nz"]) + .5) * m["height"] / m["nz"]
    return dict(name=folder.name, m=m, d=d, t=d["time_seconds"], fields=fields,
                initial_theta=initial_theta, z=z, label=label, color=color, style=style)


def metrics(case):
    m, d, t, f, z = (case[key] for key in ("m", "d", "t", "fields", "z"))
    upper_key = "upper_theta_rate_kg_K_m2_s" if "upper_theta_rate_kg_K_m2_s" in d.dtype.names else "relaxation_theta_rate_kg_K_m2_s"
    power_key = "upper_enthalpy_W_m2" if "upper_enthalpy_W_m2" in d.dtype.names else "relaxation_enthalpy_W_m2"
    measured = d["theta_column_kg_K_m2"] - d["theta_column_kg_K_m2"][0]
    expected = cumulative_trapezoid(d["bottom_theta_flux_kg_K_m2_s"] + d[upper_key], t)
    scale = np.trapezoid(np.abs(d["bottom_theta_flux_kg_K_m2_s"]) + np.abs(d[upper_key]), t)
    # Fixed physical heights allow a consistent cap comparison across resolutions.
    levels = [82500.0, 117500.0]
    contrast = np.diff(np.interp(levels, z, f["theta_profile"]))[0]
    initial_contrast = np.diff(np.interp(levels, z, case["initial_theta"]))[0]
    case.update(measured=measured, expected=expected)
    result = {
        "nx": m["nx"], "ny": m["ny"], "nz": m["nz"],
        "bottom_heat_flux_W_m2": m["heat_flux"],
        "upper_power_min_W_m2": float(d[power_key].min()),
        "upper_power_max_W_m2": float(d[power_key].max()),
        "late_mean_upper_power_W_m2": time_mean(d[power_key], t),
        "late_mean_w_rms_m_s": time_mean(d["w_rms_center_m_s"], t),
        "final_w_rms_m_s": float(d["w_rms_center_m_s"][-1]),
        "peak_w_rms_m_s": float(d["w_rms_center_m_s"].max()),
        "peak_max_abs_w_m_s": float(d["max_abs_w_m_s"].max()),
        "late_mean_ke_J_m2": time_mean(d["kinetic_energy_J_m2"], t),
        "late_mean_vorticity_rms_over_f": time_mean(d["vorticity_rms_over_f"], t),
        "cap_contrast_retention_82p5_to_117p5_km": float(contrast / initial_contrast),
        "cap_profile_departure_rms_K": float(np.sqrt(np.mean((f["theta_profile"][z>=80000]-case["initial_theta"][z>=80000])**2))),
        "final_max_upward_convective_heat_flux_W_m2": float(f["heat_flux_profile"].max()),
        "final_min_diffusive_heat_flux_W_m2": float(f["diffusive_heat_flux_profile"].min()),
        "theta_budget_final_residual_kg_K_m2": float(measured[-1]-expected[-1]),
        "theta_budget_max_abs_residual_kg_K_m2": float(np.max(np.abs(measured-expected))),
        "theta_budget_relative_residual": float(abs(measured[-1]-expected[-1])/scale) if scale>0 else None,
        "wall_seconds": m["wall_seconds"], "iterations": m["iterations"],
        "script_sha256": m.get("script_sha256"),
    }
    if m.get("cooling_model") == "fixed_power":
        loss=m["upper_cooling_flux"]
        if not np.allclose(d[power_key],-loss,rtol=0,atol=1e-10):
            raise ValueError(f"Incorrect upper power in {case['name']}")
        if not np.allclose(d["external_enthalpy_net_W_m2"],m["heat_flux"]-loss,rtol=0,atol=1e-10):
            raise ValueError(f"Incorrect prescribed net power in {case['name']}")
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--previous", required=True)
    parser.add_argument("--output", default="gallery/controlled48")
    args=parser.parse_args()
    destination=Path(args.output); destination.mkdir(parents=True,exist_ok=True)
    cases=[read_case(Path(args.campaign)/name,label,color,style) for name,label,color,style in CATALOG]
    old=read_case(Path(args.previous)/"convection48_on_p002","Previous setup, heated 0.02 K","#777777",":")
    summary={"geometry":"traditional f-plane at 55N", "duration_hours":48,
             "cases":{c["name"]:metrics(c) for c in cases}, "previous":metrics(old), "pairs":{}}
    for tag in ("p020","p002"):
        on=next(c for c in cases if c["name"]==f"controlled48_on_{tag}")
        off=next(c for c in cases if c["name"]==f"controlled48_off_{tag}")
        keys=("nx","ny","nz","height","lx","ly","latitude","rotation_rate","f","beta","seed","perturbation",
              "convective_n2","cap_n2","cap_height","viscosity","thermal_diffusivity","cap_diffusivity",
              "diffusivity_taper_width","upper_cooling_flux","cooling_model","mixing_model","script_sha256")
        if any(on["m"][key]!=off["m"][key] for key in keys):
            raise ValueError(f"Unmatched inputs: {tag}")
        a,b=(summary["cases"][c["name"]] for c in (on,off))
        summary["pairs"][tag]={"late_w_rms_change_percent":100*(a["late_mean_w_rms_m_s"]/b["late_mean_w_rms_m_s"]-1),
                                "late_ke_change_percent":100*(a["late_mean_ke_J_m2"]/b["late_mean_ke_J_m2"]-1)}
    coarse=summary["cases"]["controlled48_on_p002"]; fine=summary["cases"]["controlled48_refined"]
    summary["vertical_refinement"]={"late_w_rms_change_percent":100*(fine["late_mean_w_rms_m_s"]/coarse["late_mean_w_rms_m_s"]-1),
                                    "late_ke_change_percent":100*(fine["late_mean_ke_J_m2"]/coarse["late_mean_ke_J_m2"]-1),
                                    "note":"Different grids have different random spatial realizations; this is a bulk sensitivity check, not pointwise convergence."}
    (destination/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    names=cases[0]["d"].dtype.names
    with (destination/"diagnostics_hourly.csv").open("w") as stream:
        writer=csv.writer(stream); writer.writerow(["case",*names])
        for case in cases:
            with gzip.open(destination/(case["name"]+"_diagnostics.csv.gz"),"wt",newline="") as full:
                full_writer=csv.writer(full); full_writer.writerow(names)
                for row in case["d"]:
                    values=[row[key] for key in names]; full_writer.writerow(values)
                    if np.isclose(row["time_seconds"]%3600,0): writer.writerow([case["name"],*values])

    fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
    panels=[("w_rms_center_m_s","Vertical RMS speed [m/s]","log"),
            ("kinetic_energy_J_m2","Column kinetic energy [J/m²]","log"),
            ("cap_mean_profile_departure_rms_K","Stable-cap mean-profile RMS change [K]","linear"),
            ("upper_enthalpy_W_m2","Prescribed upper heat input [W/m²]","linear")]
    for ax,(key,title,scale) in zip(axes.flat,panels):
        for c in cases:
            ax.plot(c["t"][1:]/3600,c["d"][key][1:],c["style"],color=c["color"],label=c["label"],linewidth=1.3)
        ax.set(title=title,xlabel="Time [Earth hours]",yscale=scale); ax.grid(alpha=.2)
    axes[0,0].plot(old["t"][1:]/3600,old["d"]["w_rms_center_m_s"][1:],":",color="gray",label=old["label"])
    axes[0,0].legend(fontsize=6)
    axes[1,1].set_ylim(-6,1)
    fig.suptitle("Revised source controls: fixed cooling, reduced cap mixing, no beta or shear")
    save_figure(fig,destination,"evolution")

    fig,axes=plt.subplots(1,3,figsize=(12,5),layout="constrained")
    selected=[cases[2],cases[3],cases[4],cases[5],old]
    for c in selected:
        for ax,key in zip(axes,("theta_profile","heat_flux_profile","diffusive_heat_flux_profile")):
            ax.plot(c["fields"][key],c["z"]/1000,c["style"],color=c["color"],label=c["label"])
    axes[0].plot(cases[2]["initial_theta"],cases[2]["z"]/1000,"k--",label="Initial profile")
    for ax,label in zip(axes,("Potential temperature [K]","Resolved convective heat flux [W/m²]","Diffusive heat flux estimate [W/m²]")):
        ax.set(xlabel=label,ylabel="Height [km]"); ax.axhline(80,color="gray",linewidth=.5); ax.grid(alpha=.2)
    axes[0].legend(fontsize=6)
    fig.suptitle("48-hour profiles; the previous uniform-diffusion setup is shown in gray")
    save_figure(fig,destination,"profiles")

    fig,axes=plt.subplots(2,3,figsize=(12,7),layout="constrained")
    for ax,c in zip(axes.flat,cases):
        ax.plot(c["t"]/3600,c["measured"],color=c["color"],label="Measured change")
        ax.plot(c["t"]/3600,c["expected"],"k--",label="Integrated sources")
        ax.set(title=c["label"],xlabel="Time [Earth hours]",ylabel="Δ∫ρ θ dz [kg K/m²]"); ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    fig.suptitle("Potential-temperature inventory; equal enthalpy powers need not cancel in this variable")
    save_figure(fig,destination,"inventory")

    fig,axes=plt.subplots(1,2,figsize=(9,5),layout="constrained")
    for c in (cases[2],cases[5]):
        axes[0].plot(c["fields"]["kappa_profile"],c["z"]/1000,c["style"],color=c["color"],label=c["label"])
        axes[1].plot(c["fields"]["upper_energy_source"],c["z"]/1000,c["style"],color=c["color"])
    axes[0].set(xlabel="Thermal diffusivity [m²/s]",xscale="log",ylabel="Height [km]")
    axes[1].set(xlabel="Volumetric energy source [W/m³]",ylabel="Height [km]")
    for ax in axes: ax.axhline(80,color="gray",linewidth=.5); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7); fig.suptitle("Prescribed profiles; upper sink integrates to −5.4 W/m²")
    save_figure(fig,destination,"prescribed_profiles")

    selected=cases[2:]
    fig,axes=plt.subplots(2,4,figsize=(14,7),layout="constrained")
    limits=[max(np.max(np.abs(c["fields"]["zeta"]))/c["m"]["f"] for c in selected),
            max(np.max(np.abs(c["fields"]["w_xz"])) for c in selected)]
    for col,c in enumerate(selected):
        for row,key in enumerate(("zeta","w_xz")):
            a=c["fields"][key]/(c["m"]["f"] if row==0 else 1)
            height=c["m"]["ly"] if row==0 else c["m"]["height"]
            im=axes[row,col].imshow(a,origin="lower",extent=(0,c["m"]["lx"]/1000,0,height/1000),
                                    aspect="auto",cmap="RdBu_r",vmin=-limits[row],vmax=limits[row],interpolation="nearest")
            axes[row,col].set(xlabel="x [km]",ylabel="y [km]" if row==0 else "z [km]")
            if row==0: axes[row,col].set_title(c["label"],fontsize=8)
            if col==3: fig.colorbar(im,ax=list(axes[row,:]),shrink=.8,label="ζ/f near 80 km" if row==0 else "w [m/s]")
    fig.suptitle("48-hour fields, common color scales; no synthetic cloud rendering")
    save_figure(fig,destination,"final_fields")
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
