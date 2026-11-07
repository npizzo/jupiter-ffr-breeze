"""Analyze enhanced forcing pairs and their horizontal-resolution sensitivity."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from analyze_controlled import read_case, metrics
from analyze_neutral import snapshot, plot_sections
from analyze_turbulence import velocity_map, spectrum
from analyze_convection import save_figure, time_mean

GROUPS=[("f030", "30× flux", "#0072B2"), ("f100", "100× flux", "#D55E00"),
        ("f100_h256", "100× flux, 256² horizontal", "#CC79A7")]


def main():
    p=argparse.ArgumentParser();p.add_argument("campaign");p.add_argument("--output",default="gallery/strong48")
    args=p.parse_args();root=Path(args.campaign);dest=Path(args.output);dest.mkdir(parents=True,exist_ok=True)
    summary={"duration_hours":48,"cases":{},"pairs":{},"note":"Enhanced forcing; not Jupiter's physical heat flux. No asymptotic rapid-rotation regime is assumed."}
    pairs=[]
    for tag,label,color in GROUPS:
        factor=tag.split("_")[0];suffix="_h256" if tag.endswith("h256") else ""
        cases=[]
        for mode,style in [("heated","-"),("cooled","--")]:
            name=f"strong48_{factor}_{mode}{suffix}"
            c=read_case(root/name,label+": "+("heated + cooled" if mode=="heated" else "cooled only"),color,style)
            c["initial_map"]=velocity_map(root/name/"frame_00000.jld2")
            c["final_map"]=velocity_map(root/name/"frame_00048.jld2")
            r=metrics(c);r["late_mean_buoyancy_power_W_m2"]=time_mean(c["d"]["buoyancy_power_center_W_m2"],c["t"])
            r["final_column_ke_J_m2"]=float(c["d"]["kinetic_energy_J_m2"][-1])
            summary["cases"][name]=r;cases.append(c)
            with gzip.open(dest/(name+"_diagnostics.csv.gz"),"wt",newline="") as stream:
                w=csv.writer(stream);w.writerow(c["d"].dtype.names)
                w.writerows([[row[n] for n in c["d"].dtype.names] for row in c["d"]])
        a,b=cases
        assert all(np.array_equal(a["initial_map"][k],b["initial_map"][k]) for k in a["initial_map"])
        assert np.array_equal(a["initial_theta"],b["initial_theta"])
        keys=("nx","ny","nz","f","beta","cap_n2","convective_n2","viscosity","thermal_diffusivity","upper_cooling_flux",
              "turbulence_rossby","turbulence_seed","script_sha256","initial_turbulence_script_sha256")
        assert all(a["m"][k]==b["m"][k] for k in keys)
        ra,rb=(summary["cases"][c["name"]] for c in cases)
        summary["pairs"][tag]={key+"_change_percent":100*(ra[key]/rb[key]-1) for key in
            ("late_mean_ke_J_m2","late_mean_w_rms_m_s","final_column_ke_J_m2")}
        pairs.append(cases)
        frames=[snapshot(root/c["name"]/"frame_00048.jld2") for c in cases]
        plot_sections(cases,frames,dest,tag+"_final_sections",[c["label"] for c in cases],label+": final paired sections")
        fig,axes=plt.subplots(1,2,figsize=(9,4),layout="constrained")
        limit=max(float(np.max(abs(c["final_map"]["zeta"]))/c["m"]["f"]) for c in cases)
        for ax,c in zip(axes,cases):
            im=ax.imshow(c["final_map"]["zeta"]/c["m"]["f"],origin="lower",cmap="RdBu_r",vmin=-limit,vmax=limit,
                extent=(0,c["m"]["lx"]/1000,0,c["m"]["ly"]/1000))
            ax.set(title=c["label"],xlabel="x [km]",ylabel="y [km]")
        fig.colorbar(im,ax=list(axes),shrink=.8,label="ζ/f near 20 km")
        fig.suptitle(label+": final vorticity, common scale")
        save_figure(fig,dest,tag+"_final_vorticity")
        hours=[6,12,24,48];frames=[snapshot(root/a["name"]/f"frame_{h:05d}.jld2") for h in hours]
        assert all(frame["t"]==h*3600 for frame,h in zip(frames,hours))
        plot_sections([a]*4,frames,dest,tag+"_heated_snapshots",[f"{h} h" for h in hours],label+": heated snapshots")
        hours=[0,6,12,24,48];maps=[velocity_map(root/a["name"]/f"frame_{h:05d}.jld2") for h in hours]
        fig,axes=plt.subplots(1,5,figsize=(16,4),layout="constrained")
        limit=max(float(np.max(abs(m["zeta"]))/a["m"]["f"]) for m in maps)
        for ax,h,m in zip(axes,hours,maps):
            im=ax.imshow(m["zeta"]/a["m"]["f"],origin="lower",cmap="RdBu_r",vmin=-limit,vmax=limit,
                extent=(0,a["m"]["lx"]/1000,0,a["m"]["ly"]/1000))
            ax.set(title=f"{h} h",xlabel="x [km]",ylabel="y [km]")
        fig.colorbar(im,ax=list(axes),shrink=.8,label="ζ/f near 20 km");fig.suptitle(label+": heated vorticity, fixed scale")
        save_figure(fig,dest,tag+"_heated_vorticity")
    (dest/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout="constrained")
    for col,(cases,(_,label,_)) in enumerate(zip(pairs,GROUPS)):
        for row,(key,units) in enumerate([("kinetic_energy_J_m2","Column KE [J/m²]"),("w_rms_center_m_s","Vertical RMS speed [m/s]")]):
            for c in cases:axes[row,col].plot(c["t"]/3600,c["d"][key],c["style"],color=c["color"],label=c["label"])
            axes[row,col].set(xlabel="Time [Earth hours]",ylabel=units,yscale="log");axes[row,col].grid(alpha=.2)
        axes[0,col].set_title(label);axes[0,col].legend(fontsize=6)
    save_figure(fig,dest,"evolution")
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout="constrained")
    for ax,cases,(_,label,_) in zip(axes,pairs,GROUPS):
        for c in cases:
            k,e=spectrum(c["final_map"]);ax.loglog(k[1:],e[1:],c["style"],color=c["color"],label=c["label"])
        ax.set(title=label,xlabel="Horizontal shell index",ylabel="Shell horizontal KE [m²/s²]");ax.grid(alpha=.2);ax.legend(fontsize=6)
    fig.suptitle("Final near-20-km slice spectra; not a full-depth transfer budget")
    save_figure(fig,dest,"spectra")
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__=="__main__":main()
