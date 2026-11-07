"""Matched turbulence/thermal-forcing comparisons, with spectra and snapshots."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from analyze_controlled import read_case, metrics
from analyze_neutral import snapshot, plot_sections
from analyze_convection import save_figure, time_mean

CATALOG = [("turbulence48_heated", "Heated + cooled", "#D55E00", "-"),
           ("turbulence48_cooled", "Cooled only", "#0072B2", "--"),
           ("turbulence48_unforced", "Unforced", "#009E73", ":")]


def velocity_map(path):
    with h5py.File(path) as f:
        return {k: f[k][:] for k in ("u", "v", "w", "zeta", "horizontal_ke_profile", "vertical_ke_profile")}


def spectrum(fields):
    u, v = fields["u"], fields["v"]
    ny, nx = u.shape
    assert nx == ny
    eu = abs(np.fft.fft2(u-u.mean())/u.size)**2
    ev = abs(np.fft.fft2(v-v.mean())/v.size)**2
    modes = np.fft.fftfreq(nx)*nx
    radius = np.sqrt(modes[:, None]**2 + modes[None, :]**2)
    shells = np.rint(radius).astype(int)
    e = np.bincount(shells.ravel(), weights=((eu+ev)/2).ravel())
    assert np.isclose(e.sum(), np.mean((u-u.mean())**2+(v-v.mean())**2)/2, rtol=1e-12, atol=1e-20)
    return np.arange(len(e)), e


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", default="gallery/turbulence48")
    args=parser.parse_args()
    root, dest=Path(args.campaign),Path(args.output);dest.mkdir(parents=True,exist_ok=True)
    cases=[read_case(root/name,label,color,style) for name,label,color,style in CATALOG]
    initial_maps=[velocity_map(root/c["name"]/"frame_00000.jld2") for c in cases]
    for fields in initial_maps[1:]:
        assert all(np.array_equal(fields[k],initial_maps[0][k]) for k in fields), "Initial velocities differ"
    assert all(np.array_equal(c["initial_theta"],cases[0]["initial_theta"]) for c in cases)
    summary={"duration_hours":48,"geometry":"55N Jovian f-plane; no beta or imposed shear","cases":{}}
    finals=[]; maps=[]
    for c, initial in zip(cases,initial_maps):
        d,m,t=c["d"],c["m"],c["t"]
        assert np.isclose(m["initial_eddy_rossby"],.1,rtol=1e-12)
        assert m["initial_mass_divergence_relative"] < 1e-12
        final=velocity_map(root/c["name"]/"frame_00048.jld2");maps.append(final)
        frames=snapshot(root/c["name"]/"frame_00048.jld2");finals.append(frames)
        k,e=spectrum(final); c["spectrum"]=(k,e)
        result=metrics(c)
        result.update({"initial_column_ke_J_m2":float(d["kinetic_energy_J_m2"][0]),
            "final_column_ke_J_m2":float(d["kinetic_energy_J_m2"][-1]),
            "final_ke_retention":float(d["kinetic_energy_J_m2"][-1]/d["kinetic_energy_J_m2"][0]),
            "initial_eddy_rossby":m["initial_eddy_rossby"],
            "initial_horizontal_rms_lower_m_s":m["initial_horizontal_rms_lower_m_s"],
            "initial_vertical_rms_lower_m_s":m["initial_vertical_rms_lower_m_s"],
            "initial_mass_divergence_relative":m["initial_mass_divergence_relative"],
            "late_mean_buoyancy_power_W_m2":time_mean(d["buoyancy_power_center_W_m2"],t),
            "final_map_horizontal_ke_m2_s2":float(e.sum()),
            "final_map_fraction_ke_modes_1_2":float(e[1:3].sum()/e[1:].sum()),
            "final_map_integral_length_m":float(m["lx"]/(2*np.pi)*np.sum(e[1:]/k[1:])/e[1:].sum()),
            "spectrum_note":"Single horizontal slice near 20 km, not a full-depth energy spectrum.",
            "initial_turbulence_script_sha256":m["initial_turbulence_script_sha256"]})
        summary["cases"][c["name"]]=result
        with gzip.open(dest/(c["name"]+"_diagnostics.csv.gz"),"wt",newline="") as stream:
            w=csv.writer(stream);w.writerow(d.dtype.names);w.writerows([[row[n] for n in d.dtype.names] for row in d])
    a,b,c=[summary["cases"][x[0]] for x in CATALOG]
    summary["heating_increment"]={key+"_percent":100*(a[key]/b[key]-1) for key in ("late_mean_ke_J_m2","late_mean_w_rms_m_s","final_column_ke_J_m2")}
    summary["cooling_increment"]={key+"_percent":100*(b[key]/c[key]-1) for key in ("late_mean_ke_J_m2","late_mean_w_rms_m_s")}
    (dest/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")

    fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
    for ax,key,title in zip(axes.flat,("kinetic_energy_J_m2","w_rms_center_m_s","vorticity_rms_over_f","buoyancy_power_center_W_m2"),
        ("Column kinetic energy [J/m²]","Vertical RMS speed [m/s]","Volume RMS vertical vorticity / f","Resolved buoyancy power [W/m²]")):
        for case in cases:ax.plot(case["t"]/3600,case["d"][key],case["style"],color=case["color"],label=case["label"])
        if key in ("kinetic_energy_J_m2","w_rms_center_m_s"):ax.set_yscale("log")
        ax.set(xlabel="Time [Earth hours]",ylabel=title);ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.suptitle("Identical initial turbulence: incremental effects of thermal forcing")
    save_figure(fig,dest,"evolution")

    fig,axes=plt.subplots(1,2,figsize=(11,5),layout="constrained")
    for case,initial,final in zip(cases,initial_maps,maps):
        k,e=case["spectrum"];axes[0].loglog(k[1:],e[1:],case["style"],color=case["color"],label=case["label"])
        profile=final["horizontal_ke_profile"]+final["vertical_ke_profile"]
        axes[1].plot(profile,case["z"]/1000,case["style"],color=case["color"],label=case["label"])
    k,e=spectrum(initial_maps[0]);axes[0].loglog(k[1:],np.maximum(e[1:],1e-30),"k--",label="Initial")
    profile=initial_maps[0]["horizontal_ke_profile"]+initial_maps[0]["vertical_ke_profile"]
    axes[1].plot(profile,cases[0]["z"]/1000,"k--",label="Initial")
    axes[0].set(xlabel="Horizontal shell index kLbox/(2π)",ylabel="Shell horizontal KE [m²/s²]",ylim=(1e-12,10))
    axes[1].set(xlabel="Horizontal-face + vertical-center KE density [J/m³]",ylabel="Height [km]")
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle("48-hour scale and depth diagnostics; spectrum is one slice near 20 km")
    save_figure(fig,dest,"spectra_profiles")

    plot_sections(cases,finals,dest,"final_sections",[c["label"] for c in cases],"48-hour turbulence comparison")
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout="constrained")
    limit=max(float(np.max(abs(frame["zeta"]))/case["m"]["f"]) for frame,case in zip(maps,cases))
    for ax,case,frame in zip(axes,cases,maps):
        im=ax.imshow(frame["zeta"]/case["m"]["f"],origin="lower",cmap="RdBu_r",vmin=-limit,vmax=limit,
            extent=(0,case["m"]["lx"]/1000,0,case["m"]["ly"]/1000))
        ax.set(title=case["label"],xlabel="x [km]",ylabel="y [km]")
    fig.colorbar(im,ax=list(axes),shrink=.8,label="ζ/f near 20 km");fig.suptitle("48-hour vorticity; common color scale")
    save_figure(fig,dest,"final_vorticity")
    for case in cases:
        hours=[6,12,24,48]
        frames=[snapshot(root/case["name"]/f"frame_{h:05d}.jld2") for h in hours]
        plot_sections([case]*4,frames,dest,case["name"]+"_snapshots",[f"{h} h" for h in hours],case["label"]+": turbulence snapshots")
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__=="__main__":main()
