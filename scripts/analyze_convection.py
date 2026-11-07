"""Compare the completed, paired convection runs; publish compact gallery plots."""
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
from PIL import Image


def cumulative_trapezoid(y, t):
    return np.r_[0.0, np.cumsum((y[1:] + y[:-1]) * np.diff(t) / 2)]


def time_mean(y, t, start=86400.0):
    keep = t >= start
    return float(np.trapezoid(y[keep], t[keep]) / (t[keep][-1] - t[keep][0]))


def save_figure(fig, destination, name):
    full = destination / f"{name}_full.png"
    fig.savefig(full, dpi=130)
    plt.close(fig)
    original = Image.open(full).convert("RGB")
    width = min(840, original.width)
    preview = destination / f"{name}.png"
    while True:
        size = (width, round(original.height * width / original.width))
        original.resize(size, Image.Resampling.LANCZOS).quantize(
            colors=96, dither=Image.Dither.NONE
        ).save(preview, optimize=True)
        if preview.stat().st_size < 220_000 or width <= 500:
            break
        width = round(width * .85)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", default="gallery/convection48")
    args = parser.parse_args()
    root, destination = Path(args.campaign), Path(args.output)
    destination.mkdir(parents=True, exist_ok=True)
    cases = []
    summary = {"duration_hours": 48, "geometry": "traditional f-plane", "cases": {}, "pairs": {}}
    for tag, amplitude, color in [("p020", .2, "#0072B2"), ("p002", .02, "#D55E00")]:
        for mode in ("on", "off"):
            name = f"convection48_{mode}_{tag}"
            folder = root / name
            metadata = tomllib.loads((folder / "metadata.toml").read_text())
            if metadata["status"] != "COMPLETE":
                raise ValueError(f"Incomplete run: {name}")
            data = np.genfromtxt(folder / "diagnostics.csv", delimiter=",", names=True)
            t = data["time_seconds"]
            if t[-1] != metadata["stop_seconds"] or len(t) < 3:
                raise ValueError(f"Incomplete diagnostics: {name}")
            files = sorted(folder.glob("frame_*.jld2"))
            with h5py.File(files[-1]) as handle:
                fields = {key: handle[key][:] for key in (
                    "theta_profile", "heat_flux_profile", "diffusive_heat_flux_profile", "zeta", "w_xz"
                )}
            # JLD2 stores Julia ranges as structured scalars. Reconstruct this
            # uniform grid from its recorded extent and number of cells.
            fields["z"] = (np.arange(metadata["nz"]) + .5) * metadata["height"] / metadata["nz"]
            with h5py.File(files[0]) as handle:
                initial_theta = handle["theta_profile"][:]
            theta_change = data["theta_column_kg_K_m2"] - data["theta_column_kg_K_m2"][0]
            bottom = cumulative_trapezoid(data["bottom_theta_flux_kg_K_m2_s"], t)
            relaxation = cumulative_trapezoid(data["relaxation_theta_rate_kg_K_m2_s"], t)
            residual = theta_change - bottom - relaxation
            source_scale = abs(bottom[-1]) + np.trapezoid(np.abs(data["relaxation_theta_rate_kg_K_m2_s"]), t)
            metrics = {
                "heat_flux_W_m2": metadata["heat_flux"],
                "perturbation_parameter_K": amplitude,
                "initial_convective_rossby": np.sqrt(-metadata["convective_n2"]) / metadata["f"],
                "peak_max_abs_w_m_s": float(np.max(data["max_abs_w_m_s"])),
                "final_w_rms_m_s": float(data["w_rms_center_m_s"][-1]),
                "late_mean_w_rms_m_s": time_mean(data["w_rms_center_m_s"], t),
                "late_mean_ke_J_m2": time_mean(data["kinetic_energy_J_m2"], t),
                "final_vorticity_rms_over_f": float(data["vorticity_rms_over_f"][-1]),
                "late_mean_vorticity_rms_over_f": time_mean(data["vorticity_rms_over_f"], t),
                "late_mean_relaxation_enthalpy_W_m2": time_mean(data["relaxation_enthalpy_W_m2"], t),
                "final_max_resolved_upward_heat_flux_W_m2": float(np.max(fields["heat_flux_profile"])),
                "final_min_diffusive_heat_flux_W_m2": float(np.min(fields["diffusive_heat_flux_profile"])),
                "bottom_energy_input_J_m2": metadata["heat_flux"] * t[-1],
                "signed_relaxation_enthalpy_input_J_m2": float(np.trapezoid(data["relaxation_enthalpy_W_m2"], t)),
                "theta_column_change_kg_K_m2": float(theta_change[-1]),
                "integrated_bottom_theta_input_kg_K_m2": float(bottom[-1]),
                "integrated_relaxation_theta_input_kg_K_m2": float(relaxation[-1]),
                "theta_budget_final_residual_kg_K_m2": float(residual[-1]),
                "theta_budget_relative_residual": float(abs(residual[-1]) / max(source_scale, 1e-30)),
                "wall_seconds": metadata["wall_seconds"],
                "post_warmup_wall_seconds": metadata["post_warmup_wall_seconds"],
                "iterations": metadata["iterations"],
                "script_sha256": metadata["script_sha256"],
            }
            summary["cases"][name] = metrics
            cases.append(dict(name=name, mode=mode, tag=tag, amplitude=amplitude, color=color,
                              style="-" if mode == "on" else "--", label=f"{amplitude:g} K, heating {mode}",
                              metadata=metadata, data=data, t=t, fields=fields, initial_theta=initial_theta,
                              theta_change=theta_change, expected_theta=bottom + relaxation, residual=residual))

    for tag in ("p020", "p002"):
        on = summary["cases"][f"convection48_on_{tag}"]
        off = summary["cases"][f"convection48_off_{tag}"]
        on_case = next(c for c in cases if c["tag"] == tag and c["mode"] == "on")
        off_case = next(c for c in cases if c["tag"] == tag and c["mode"] == "off")
        for key in ("nx", "ny", "nz", "lx", "ly", "height", "latitude", "rotation_rate", "f", "beta",
                    "convective_n2", "cap_n2", "cap_height", "perturbation", "seed", "viscosity",
                    "thermal_diffusivity", "cooling_time", "sponge_time", "script_sha256"):
            if on_case["metadata"][key] != off_case["metadata"][key]:
                raise ValueError(f"Unmatched pair {tag}: {key}")
        if not np.allclose(on_case["initial_theta"], off_case["initial_theta"], rtol=0, atol=1e-12):
            raise ValueError(f"Initial mean profiles differ for {tag}")
        summary["pairs"][tag] = {
            "late_w_rms_change_percent": 100 * (on["late_mean_w_rms_m_s"] / off["late_mean_w_rms_m_s"] - 1),
            "late_ke_change_percent": 100 * (on["late_mean_ke_J_m2"] / off["late_mean_ke_J_m2"] - 1),
        }
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (destination / "diagnostics_hourly.csv").open("w") as stream:
        writer = csv.writer(stream)
        names = cases[0]["data"].dtype.names
        writer.writerow(["case", *names])
        for case in cases:
            for row in case["data"]:
                if np.isclose(row["time_seconds"] % 3600, 0):
                    writer.writerow([case["name"], *(row[key] for key in names)])
    with gzip.open(destination / "diagnostics_5min.csv.gz", "wt", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case", *names])
        for case in cases:
            for row in case["data"]:
                writer.writerow([case["name"], *(row[key] for key in names)])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    panels = [("w_rms_center_m_s", "Vertical RMS speed [m/s]", "log"),
              ("kinetic_energy_J_m2", "Column kinetic energy [J/m²]", "log"),
              ("vorticity_rms_over_f", "Volume RMS vertical vorticity / f", "linear"),
              ("relaxation_enthalpy_W_m2", "Upper-layer heat input, signed [W/m²]", "linear")]
    for ax, (key, title, scale) in zip(axes.flat, panels):
        for case in cases:
            ax.plot(case["t"][1:]/3600, case["data"][key][1:], case["style"], color=case["color"], label=case["label"])
        ax.set(title=title, xlabel="Time [Earth hours]", yscale=scale)
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    axes[1, 1].axhline(5.4, color="black", linestyle=":", label="Bottom heating when on: 5.4 W/m²")
    axes[1, 1].legend(fontsize=7)
    fig.suptitle("Controlled dry convection: fixed Jovian rotation, f-plane at 55°N")
    save_figure(fig, destination, "evolution")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for ax, case in zip(axes.flat, cases):
        ax.plot(case["t"]/3600, case["theta_change"], color=case["color"], label="Measured column change")
        ax.plot(case["t"]/3600, case["expected_theta"], "k--", label="Integrated bottom + relaxation sources")
        ax.set(title=case["label"], xlabel="Time [Earth hours]", ylabel="Δ∫ρ θ dz [kg K/m²]")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Potential-temperature inventory: source integral uses 5-minute samples")
    save_figure(fig, destination, "thermal_budget")

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), layout="constrained")
    for case in cases:
        fields = case["fields"]
        for ax, key in zip(axes, ("theta_profile", "heat_flux_profile", "diffusive_heat_flux_profile")):
            ax.plot(fields[key], fields["z"]/1000, case["style"], color=case["color"], label=case["label"])
    axes[0].plot(cases[0]["initial_theta"], cases[0]["fields"]["z"]/1000, "k:", label="Initial mean profile")
    for ax, label in zip(axes, ("Potential temperature [K]", "Resolved convective heat flux [W/m²]", "Diffusive heat flux estimate [W/m²]")):
        ax.set(xlabel=label, ylabel="Height [km]")
        ax.axhline(80, color="gray", linewidth=.5)
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle("Profiles at 48 hours; fluxes are cell-centered diagnostic estimates")
    save_figure(fig, destination, "final_profiles")

    fig, axes = plt.subplots(2, 4, figsize=(14, 7), layout="constrained")
    vort_limit = max(np.max(np.abs(c["fields"]["zeta"])) / c["metadata"]["f"] for c in cases)
    w_limit = max(np.max(np.abs(c["fields"]["w_xz"])) for c in cases)
    for column, case in enumerate(cases):
        for row, key, limit, height in [(0, "zeta", vort_limit, 1200), (1, "w_xz", w_limit, 120)]:
            field = case["fields"][key] / (case["metadata"]["f"] if row == 0 else 1)
            im = axes[row, column].imshow(field, origin="lower", extent=(0,1200,0,height), aspect="auto",
                                         cmap="RdBu_r", vmin=-limit, vmax=limit, interpolation="nearest")
            axes[row, column].set(xlabel="x [km]", ylabel="y [km]" if row == 0 else "z [km]")
            if row == 0:
                axes[row, column].set_title(case["label"], fontsize=9)
            if column == 3:
                fig.colorbar(im, ax=list(axes[row, :]), shrink=.8, label="ζ/f" if row == 0 else "w [m/s]")
    fig.suptitle("48-hour fields: common color scales across cases; no clouds modeled")
    save_figure(fig, destination, "final_fields")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
