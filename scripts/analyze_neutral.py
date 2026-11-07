"""Compare neutral-layer forcing controls and plot true solver snapshots."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from analyze_controlled import read_case, metrics
from analyze_convection import save_figure, time_mean

CATALOG = [
    ("neutral48_heated", "Heated + cooled", "#D55E00", "-"),
    ("neutral48_cooled", "Cooled only", "#0072B2", "--"),
    ("neutral48_unforced", "Unforced", "#009E73", ":"),
    ("neutral48_refined", "Heated + cooled, Nz=128", "#CC79A7", "-."),
]


def snapshot(path):
    with h5py.File(path) as f:
        d = {key: f[key][:] for key in ("w_xz", "theta_xz", "theta_profile", "mean_N2_profile", "zeta")}
        d["t"] = float(f["t"][()])
    # h5py sees the reversed dimensions of Julia's x-z matrices.
    assert d["theta_xz"].shape[0] == len(d["theta_profile"])
    assert d["w_xz"].shape[0] == len(d["theta_profile"])+1
    # Saved vertical velocity lives on faces; center it for matching sections.
    d["w_xz"] = (d["w_xz"][:-1]+d["w_xz"][1:])/2
    d["theta_anomaly"] = d["theta_xz"] - d["theta_profile"][:, None]
    return d


def plot_sections(cases, frames, destination, name, titles, heading):
    fig, axes = plt.subplots(2, len(frames), figsize=(14, 7), layout="constrained", squeeze=False)
    for row, key, units in [(0, "w_xz", "w [m/s]"), (1, "theta_anomaly", "θ − horizontal mean θ [K]")]:
        limit = max(float(np.max(np.abs(frame[key]))) for frame in frames)
        limit = max(limit, np.finfo(float).tiny)
        for col, (case, frame, title) in enumerate(zip(cases, frames, titles)):
            im = axes[row, col].imshow(frame[key], origin="lower", aspect="auto", cmap="RdBu_r",
                extent=(0, case["m"]["lx"]/1000, 0, case["m"]["height"]/1000),
                vmin=-limit, vmax=limit, interpolation="nearest")
            axes[row, col].axhline(80, color="gray", lw=.5)
            axes[row, col].set(xlabel="x [km]", ylabel="z [km]")
            if row == 0:
                axes[row, col].set_title(title, fontsize=9)
        fig.colorbar(im, ax=list(axes[row, :]), shrink=.8, label=units)
    fig.suptitle(heading + "\nEach row shares one color scale; actual dynamical fields, not cloud images")
    save_figure(fig, destination, name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", default="gallery/neutral48")
    args = parser.parse_args()
    root, dest = Path(args.campaign), Path(args.output)
    dest.mkdir(parents=True, exist_ok=True)
    cases = [read_case(root/name, label, color, style) for name, label, color, style in CATALOG]
    comparison_keys = set(cases[0]["m"]) - {"name", "heat_flux", "upper_cooling_flux", "wall_seconds",
        "post_warmup_wall_seconds", "configuration_sha256", "discrete_fixed_upper_power_W_m2"}
    for c in cases[1:3]:
        assert all(c["m"].get(k) == cases[0]["m"][k] for k in comparison_keys), "Unmatched control configuration"
    summary = {"duration_hours": 48, "geometry": "traditional Jovian f-plane at 55N", "cases": {}}
    finals = []
    for c in cases:
        m, d, t = c["m"], c["d"], c["t"]
        assert m["convective_n2"] == 0 and m["beta"] == 0 and m["jet_speed"] == 0
        assert np.max(abs(c["initial_theta"][c["z"] < 80000] - m["theta0"])) < 1e-9
        result = metrics(c)
        frame = snapshot(root/c["name"]/"frame_00048.jld2")
        finals.append(frame)
        lower_section = frame["w_xz"][c["z"] < 40000]
        spectrum = np.sum(abs(np.fft.rfft(lower_section-lower_section.mean(axis=1, keepdims=True), axis=1))**2, axis=0)
        # Account for the conjugate half of the real-valued section spectrum.
        spectrum[1:-1] *= 2
        modes = np.arange(len(spectrum))
        peak_mode = int(np.argmax(spectrum[1:])+1)
        result.update({
            "late_mean_buoyancy_power_W_m2": time_mean(d["buoyancy_power_center_W_m2"], t),
            "late_mean_theta_variance_K2": time_mean(d["theta_variance_K2"], t),
            "final_lower_layer_min_mean_N2_s2": float(frame["mean_N2_profile"][c["z"] < 60000].min()),
            "final_lower_layer_theta_range_K": float(np.ptp(frame["theta_profile"][c["z"] < 60000])),
            "final_section_max_abs_w_m_s": float(np.max(abs(frame["w_xz"]))),
            "final_section_max_abs_theta_anomaly_K": float(np.max(abs(frame["theta_anomaly"]))),
            "final_lower_section_peak_horizontal_wavelength_km": m["lx"]/peak_mode/1000,
            "final_lower_section_power_fraction_wavelength_below_4dx": float(spectrum[modes > m["nx"]/4].sum()/spectrum[1:].sum()),
            "section_spectrum_note": "One x-z slice below 40 km, not a full 3D spectrum or convergence test.",
        })
        summary["cases"][c["name"]] = result
        with gzip.open(dest/(c["name"]+"_diagnostics.csv.gz"), "wt", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(d.dtype.names)
            writer.writerows([[row[k] for k in d.dtype.names] for row in d])
    heated, cooled, unforced, refined = [summary["cases"][c["name"]] for c in cases]
    summary["comparison"] = {
        "heated_minus_cooled_late_w_rms_m_s": heated["late_mean_w_rms_m_s"]-cooled["late_mean_w_rms_m_s"],
        "heated_minus_unforced_late_w_rms_m_s": heated["late_mean_w_rms_m_s"]-unforced["late_mean_w_rms_m_s"],
        "refinement_late_w_rms_change_percent": 100*(refined["late_mean_w_rms_m_s"]/heated["late_mean_w_rms_m_s"]-1),
        "refinement_note": "Different spatial noise on the refined grid; a bulk sensitivity check, not a convergence demonstration.",
    }
    (dest/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for ax, key, label in zip(axes.flat,
            ("w_rms_center_m_s", "kinetic_energy_J_m2", "theta_variance_K2", "buoyancy_power_center_W_m2"),
            ("Vertical RMS speed [m/s]", "Column kinetic energy [J/m²]", "Horizontal θ variance [K²]", "Resolved buoyancy power [W/m²]")):
        for c in cases:
            ax.plot(c["t"][1:]/3600, c["d"][key][1:], c["style"], color=c["color"], label=c["label"])
        if key != "buoyancy_power_center_W_m2":
            ax.set_yscale("log")
        ax.set(xlabel="Time [Earth hours]", ylabel=label)
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Neutral initial lower layer: weak physical forcing and matched controls")
    save_figure(fig, dest, "evolution")

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), layout="constrained")
    for c, frame in zip(cases, finals):
        z = c["z"]/1000
        for ax, value in zip(axes, (frame["theta_profile"]-c["m"]["theta0"], frame["mean_N2_profile"], c["fields"]["heat_flux_profile"])):
            ax.plot(value, z, c["style"], color=c["color"], label=c["label"])
    for ax, label in zip(axes, ("Mean θ − initial lower-layer θ [K]", "Mean N² [s⁻²]", "Resolved convective heat flux [W/m²]")):
        ax.set(xlabel=label, ylabel="Height [km]", ylim=(0, 60))
        ax.axvline(0, color="gray", lw=.5)
        ax.grid(alpha=.2)
    # Show the initially neutral lower part without the much warmer stable cap.
    for ax in axes:
        lines = ax.lines[:len(cases)]
        values = np.concatenate([line.get_xdata()[c["z"] <= 60000] for line, c in zip(lines, cases)])
        lo, hi = values.min(), values.max()
        pad = max((hi-lo)*.05, 1e-20)
        ax.set_xlim(lo-pad, hi+pad)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
    axes[0].legend(fontsize=7)
    fig.suptitle("48-hour lower-layer profiles; excludes the stable cap")
    save_figure(fig, dest, "lower_profiles")

    plot_sections(cases, finals, dest, "final_sections", [c["label"] for c in cases], "48-hour comparison")
    fig, axes = plt.subplots(1, 4, figsize=(14, 4), layout="constrained")
    limit = max(float(np.max(abs(frame["zeta"]))/c["m"]["f"]) for c, frame in zip(cases, finals))
    for ax, c, frame in zip(axes, cases, finals):
        im = ax.imshow(frame["zeta"]/c["m"]["f"], origin="lower", cmap="RdBu_r",
            extent=(0, c["m"]["lx"]/1000, 0, c["m"]["ly"]/1000),
            vmin=-limit, vmax=limit, interpolation="nearest")
        ax.set(title=c["label"], xlabel="x [km]", ylabel="y [km]")
    fig.colorbar(im, ax=list(axes), shrink=.8, label="ζ/f near 80 km")
    fig.suptitle("48-hour horizontal vorticity; common physical scale across cases")
    save_figure(fig, dest, "final_vorticity")
    for c in cases:
        hours = [6, 12, 24, 48]
        frames = [snapshot(root/c["name"]/f"frame_{hour:05d}.jld2") for hour in hours]
        assert all(frame["t"] == hour*3600 for frame, hour in zip(frames, hours))
        plot_sections([c]*4, frames, dest, c["name"]+"_snapshots", [f"{hour} h" for hour in hours], c["label"]+": snapshots")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
