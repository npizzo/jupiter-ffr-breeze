"""Analyze completed fixed-vortex response cases and render on Unity only.

Usage: python analyze.py CAMPAIGN_ROOT NEW_PUBLICATION_DIRECTORY
No state evolution, parameter fitting, or automatic continuation occurs here.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm
from matplotlib.ticker import FuncFormatter
from PIL import Image
import imageio_ffmpeg

F_ABS = 0.0002880957739764384
HEIGHT_M = 120000.
HEIGHTS = np.array([20000., 60000., 100000.])
SEEDS = (20261118, 20261119)
FRAME_TIMES = np.arange(21)*.5
TERMS = ("bt_power", "buoyancy_power", "viscous_power", "thermal_gradient_power",
         "thermal_diffusion_power", "thermal_relaxation_power")


def trapz(y, x):
    return float(np.sum(.5*(y[1:]+y[:-1])*np.diff(x)))


def normalized_difference(a, b, floor=1e-300):
    return float(np.max(abs(a-b))/max(float(np.max(abs(a))), float(np.max(abs(b))), floor))


def arithmetic_check(actual, expected, scale):
    error = float(np.max(abs(actual-expected)))
    assert error <= 1e-10*max(float(scale), np.finfo(float).tiny), (error, scale)


def read_cases(root):
    campaign = json.loads((root/"campaign.json").read_text())
    assert campaign["status"] == "COMPLETED" and campaign["target"] == 10.
    assert len(campaign["cases"]) == 18
    cases = {}
    for item in campaign["cases"]:
        path = root/item["name"]
        c = json.loads((path/"config.json").read_text())
        status = json.loads((path/"status.json").read_text())
        assert status["status"] == "COMPLETED" and abs(status["t_fabs"]-10.) < 1e-9
        assert c["target_time_fabs"] == 10. and not c["pilot"]
        assert c["code_sha256"] == campaign["code_sha256"]
        assert np.array_equal(c["heights_m"], HEIGHTS)
        a = np.genfromtxt(path/"scalars.csv", delimiter=",", names=True)
        cadence = float(c["scalar_interval_fabs"])
        assert cadence > 0. and abs(round(10./cadence)*cadence-10.) < 1e-9
        times = np.arange(round(10./cadence)+1)*cadence
        assert a.shape == times.shape and np.allclose(a["t_fabs"], times, rtol=0, atol=1e-9)
        assert np.min(abs(a["t_fabs"]-5.)) < 1e-9, "Late-window boundary must be an actual saved sample"
        assert all(np.isfinite(a[k]).all() for k in a.dtype.names)
        assert np.all(a["kinetic"] > 0.) and np.all(a["thermal_variance"] >= 0.)
        assert np.max(a["perturbation_rms_over_base"]) < .01
        k0 = float(a["kinetic"][0])
        kinetic_expected = -a["integral_bt_power"]+a["integral_buoyancy_power"]+a["integral_viscous_power"]
        thermal_expected = sum(a["integral_"+k] for k in TERMS[3:])
        arithmetic_check(a["kinetic_change"], a["kinetic"]-k0, np.max(a["kinetic"]))
        arithmetic_check(a["integrated_kinetic_residual"], a["kinetic_change"]-kinetic_expected,
                         k0+sum(np.max(abs(a["integral_"+k])) for k in TERMS[:3]))
        arithmetic_check(a["integrated_thermal_variance_residual"], a["thermal_variance_change"]-thermal_expected,
                         np.max(a["thermal_variance"])+sum(np.max(abs(a["integral_"+k])) for k in TERMS[3:]))
        rows = list(csv.DictReader((path/"frames.csv").open()))
        assert len(rows) == len(FRAME_TIMES)
        assert np.allclose([float(r["t_fabs"]) for r in rows], FRAME_TIMES, rtol=0, atol=1e-9)
        key = (int(c["selected_complete_vertical_modes"]), bool(c["convective_gradient_removed"]),
               int(c["seed"]), int(c["vortex_sign"]), int(c["f_over_fabs"]))
        assert key not in cases
        cases[key] = {"name": item["name"], "path": path, "config": c,
                      "data": a, "k0": k0, "frames": [path/r["file"] for r in rows]}
    expected = {(m, control, seed, sign, 1) for m in (1, 2) for control in (False, True)
                for seed in SEEDS for sign in (-1, 1)}
    expected |= {(m, False, SEEDS[0], -1, -1) for m in (1, 2)}
    assert set(cases) == expected
    k0s = [case["k0"] for case in cases.values()]
    assert max(k0s)/min(k0s)-1 < 1e-10
    return campaign, cases


def summarize(case):
    a, k0 = case["data"], case["k0"]
    late = a["t_fabs"] >= 5.-1e-9
    t, p, k = a["t_fabs"][late], a["bt_power"][late], a["kinetic"][late]
    ix = int(np.flatnonzero(late)[0])
    work = float(a["integral_bt_power"][-1]-a["integral_bt_power"][ix])
    dt = float(t[-1]-t[0])
    kinetic_scale = k0+sum(float(np.max(abs(a["integral_"+name]))) for name in TERMS[:3])
    thermal_scale = max(float(np.max(a["thermal_variance"])),
                        sum(float(np.max(abs(a["integral_"+name]))) for name in TERMS[3:]), np.finfo(float).tiny)
    residual = float(np.max(abs(a["integrated_kinetic_residual"])))
    thermal_residual = float(np.max(abs(a["integrated_thermal_variance_residual"])))
    return {"case": case["name"], "initial_kinetic": k0,
            "final_kinetic": float(a["kinetic"][-1]), "amplification_final": float(a["kinetic"][-1]/k0),
            "final_bt_power": float(a["bt_power"][-1]),
            "integrated_virtual_bt_work": float(a["integral_bt_power"][-1]),
            "integrated_virtual_bt_work_over_K0": float(a["integral_bt_power"][-1]/k0),
            "late_window_fabs_t": [5., 10.], "late_virtual_bt_work": work,
            "late_mean_bt_power": work/dt, "late_mean_bt_power_over_K0": work/(k0*dt),
            "late_mean_bt_power_over_K": trapz(p/k, t)/dt,
            "late_fraction_samples_positive_bt_power": float(np.mean(p > 0.)),
            "late_fraction_samples_growing_kinetic": float(np.mean(a["kinetic_rhs_rate"][late] > 0.)),
            "late_trapezoidal_minus_RK_work": trapz(p, t)-work,
            "kinetic_residual_final": float(a["integrated_kinetic_residual"][-1]),
            "kinetic_residual_maximum_absolute": residual,
            "kinetic_residual_max_over_K0": residual/k0,
            "kinetic_residual_max_over_work_scale": residual/kinetic_scale,
            "thermal_variance_residual_final": float(a["integrated_thermal_variance_residual"][-1]),
            "thermal_variance_residual_maximum_absolute": thermal_residual,
            "thermal_variance_residual_max_over_work_scale": thermal_residual/thermal_scale,
            "maximum_relative_divergence": float(np.max(a["divergence_relative"])),
            "maximum_perturbation_rms_over_background": float(np.max(a["perturbation_rms_over_base"])),
            "final_integrated_terms": {name: float(a["integral_"+name][-1]) for name in TERMS}}


def read_planes(path):
    with np.load(path) as d:
        assert np.array_equal(d["heights_m"], HEIGHTS)
        q = np.asarray(d["perturbation_zeta_over_fabs"], float)
        b = np.asarray(d["buoyancy_over_fabs2_H"], float)*F_ABS**2*HEIGHT_M
        assert q.shape[0] == b.shape[0] == 3 and q.shape == b.shape
        assert np.isfinite(q).all() and np.isfinite(b).all()
        return q, b


def reflection_checks(cases):
    results = {}
    for m in (1, 2):
        plus = cases[m, False, SEEDS[0], 1, 1]
        mirror = cases[m, False, SEEDS[0], -1, -1]
        columns = ("kinetic", "thermal_variance", "bt_power", "buoyancy_power", "viscous_power", "integral_bt_power")
        errors = {k: normalized_difference(plus["data"][k], mirror["data"][k], plus["k0"]*1e-12) for k in columns}
        field_errors = {"perturbation_vorticity": 0., "buoyancy": 0.}
        for p, r in zip(plus["frames"], mirror["frames"]):
            qp, bp = read_planes(p)
            qr, br = read_planes(r)
            reflect = (-np.arange(qp.shape[1])) % qp.shape[1]
            field_errors["perturbation_vorticity"] = max(field_errors["perturbation_vorticity"], normalized_difference(qr, -qp[:, reflect, :]))
            field_errors["buoyancy"] = max(field_errors["buoyancy"], normalized_difference(br, bp[:, reflect, :]))
        results[str(m)] = {"scalar_relative_differences": errors,
                           "Float32_plane_relative_differences": field_errors,
                           "numerically_consistent": max(errors.values()) < 1e-7 and max(field_errors.values()) < 2e-6}
    return results


def comparisons(cases, summaries):
    rows = []
    for m in (1, 2):
        for seed in SEEDS:
            hplus, hminus = [summaries[cases[m, False, seed, sign, 1]["name"]] for sign in (1, -1)]
            nplus, nminus = [summaries[cases[m, True, seed, sign, 1]["name"]] for sign in (1, -1)]
            row = {"selected_modes": m, "seed": seed, "window_fabs_t": [5., 10.]}
            for metric in ("late_mean_bt_power", "late_mean_bt_power_over_K0", "late_mean_bt_power_over_K"):
                hd, nd = hplus[metric]-hminus[metric], nplus[metric]-nminus[metric]
                row[metric] = {"heated_cyclone_minus_anticyclone": hd,
                               "gradient_removed_cyclone_minus_anticyclone": nd,
                               "heated_minus_gradient_removed_difference": hd-nd}
            # Numerical accounting tolerance is not sampling uncertainty.
            tol = max(10*sum(s["kinetic_residual_max_over_K0"] for s in (hplus, hminus, nplus, nminus))/5., 1e-10)
            effect = row["late_mean_bt_power_over_K0"]["heated_minus_gradient_removed_difference"]
            row["work_over_K0_accounting_tolerance"] = tol
            row["effect_sign_above_accounting_tolerance"] = int(np.sign(effect)) if abs(effect) > tol else 0
            rows.append(row)
    flags = []
    for m in (1, 2):
        signs = [r["effect_sign_above_accounting_tolerance"] for r in rows if r["selected_modes"] == m]
        if len(set(signs)) > 1:
            flags.append(f"Seed dependence: the convection-specific work difference changes sign or becomes unresolved for {m} retained mode(s).")
    for seed in SEEDS:
        signs = [r["effect_sign_above_accounting_tolerance"] for r in rows if r["seed"] == seed]
        if len(set(signs)) > 1:
            flags.append(f"Vertical truncation dependence: the convection-specific work difference changes sign or becomes unresolved for seed {seed}.")
    if all(r["effect_sign_above_accounting_tolerance"] == 0 for r in rows):
        flags.append("No convection-specific signed work difference exceeds the numerical accounting tolerance in this test.")
    return rows, flags


def plots(cases, out):
    labels = {(False, 1): "Convective gradient: cyclone", (False, -1): "Convective gradient: anticyclone",
              (True, 1): "Gradient removed: cyclone", (True, -1): "Gradient removed: anticyclone"}
    for m in (1, 2):
        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True, constrained_layout=True)
        for col, seed in enumerate(SEEDS):
            for control in (False, True):
                for sign in (1, -1):
                    c = cases[m, control, seed, sign, 1]
                    a, k0 = c["data"], c["k0"]
                    values = [a["bt_power"], a["bt_power"]/a["kinetic"], a["bt_power"]/k0]
                    for row, value in enumerate(values):
                        axes[row, col].plot(a["t_fabs"], value, color="#b62d38" if sign > 0 else "#236baf",
                                            ls="--" if control else "-", label=labels[control, sign], lw=1.5)
            axes[0, col].set_title(f"Seed {seed}")
        for row, label in enumerate(("Virtual barotropic power", "Power / current disturbance energy", "Power / initial disturbance energy")):
            axes[row, 0].set_ylabel(label)
        for ax in axes.flat:
            ax.axhline(0, color="black", lw=.6)
            ax.axvspan(5., 10., color="#888888", alpha=.09)
            ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        for ax in axes[-1]:
            ax.set_xlabel("|f| t")
        fig.suptitle(f"{m} selected convective mode(s): amplification and reinforcement are separate\n"
                     "Prescribed vortex remains fixed; positive power is prospective reinforcement", fontsize=13)
        fig.savefig(out/f"response_work_modes{m}.png", dpi=140)
        plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, constrained_layout=True)
    for col, m in enumerate((1, 2)):
        for seed in SEEDS:
            for control in (False, True):
                for sign in (1, -1):
                    c = cases[m, control, seed, sign, 1]
                    a, k0 = c["data"], c["k0"]
                    label = f"{seed}: {'gradient removed' if control else 'convective'}, {'cyclone' if sign > 0 else 'anticyclone'}"
                    axes[0, col].plot(a["t_fabs"], a["kinetic"]/k0, label=label, alpha=.8)
                    axes[1, col].plot(a["t_fabs"], a["integrated_kinetic_residual"]/k0, alpha=.8)
        axes[0, col].set_title(f"{m} selected convective mode(s)")
        axes[0, col].set_yscale("log")
        axes[0, col].legend(fontsize=6)
    axes[0, 0].set_ylabel("Disturbance kinetic energy / initial energy")
    axes[1, 0].set_ylabel("Integrated kinetic residual / initial energy")
    for ax in axes.flat:
        ax.grid(alpha=.2)
        ax.set_xlabel("|f| t")
    fig.suptitle("Measured amplification and retained RK-stage budget residual\n"
                 "Finite-step power integration is not exact quadratic-energy accounting", fontsize=13)
    fig.savefig(out/"response_amplification_and_residual.png", dpi=140)
    plt.close(fig)


def render(cases, out, field_index):
    chosen = [cases[2, False, SEEDS[0], sign, 1] for sign in (1, -1)]
    fields = [[read_planes(p)[field_index] for p in c["frames"]] for c in chosen]
    limit = max(float(np.max(abs(a))) for series in fields for a in series)
    assert limit > 0 and np.isfinite(limit)
    width = limit/100.
    norm = AsinhNorm(linear_width=width, vmin=-limit, vmax=limit)
    is_vort = field_index == 0
    name = "convective_response_vorticity" if is_vort else "convective_response_buoyancy"
    title = "Perturbation vertical vorticity" if is_vort else "Actual buoyancy perturbation — not passive tracer"
    units = "ζ′ / |f|" if is_vort else "b′ (m s⁻²)"
    length_km = chosen[0]["config"]["length_over_H"]*HEIGHT_M/1000
    extent = [-length_km/2, length_km/2]*2
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 9.1), constrained_layout=True)
    images = []
    for row in range(2):
        line = []
        for col in range(3):
            im = axes[row, col].imshow(fields[row][0][col].T, origin="lower", extent=extent,
                                       cmap="RdBu_r", norm=norm, interpolation="nearest")
            axes[row, col].set_title(f"{'Cyclone' if row == 0 else 'Anticyclone'} · z = {HEIGHTS[col]/1000:g} km")
            axes[row, col].set_xlabel("x (km)")
            axes[row, col].set_ylabel("y (km)")
            line.append(im)
        images.append(line)
    colorbar = fig.colorbar(images[0][0], ax=axes.ravel().tolist(), shrink=.85, pad=.02,
                           ticks=np.array([-1., -.1, 0., .1, 1.])*limit)
    colorbar.set_label(units+" · fixed asinh scale")
    colorbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: "0" if x == 0 else f"{x:.2e}"))
    heading = fig.suptitle("", fontsize=13)
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    plt.rcParams["animation.ffmpeg_path"] = ffmpeg
    movie = out/(name+".mp4")
    writer = FFMpegWriter(fps=4, codec="libx264", bitrate=2600,
                         extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart", "-threads", "1"])
    with writer.saving(fig, str(movie), dpi=120):
        for j, t in enumerate(FRAME_TIMES):
            for row in range(2):
                for col in range(3):
                    images[row][col].set_data(fields[row][j][col].T)
            heading.set_text(f"{title} · |f|t = {t:.1f} ({t/F_ABS/3600:.2f} hours)\n"
                             f"Two convective modes · seed {SEEDS[0]} · three reconstructed heights · fixed background vortices")
            writer.grab_frame()
        fig.savefig(out/(name+".png"), dpi=120)
    plt.close(fig)
    subprocess.run([ffmpeg, "-v", "error", "-i", str(movie), "-f", "null", "-"],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    decoded = imageio_ffmpeg.read_frames(str(movie), pix_fmt="rgb24")
    metadata = next(decoded)
    count = sum(1 for _ in decoded)
    assert count == len(FRAME_TIMES), count
    with Image.open(out/(name+".png")) as img:
        img.load()
        assert img.width > 1000 and img.height > 600 and np.std(np.asarray(img)) > 1.
    return {"file": movie.name, "poster": name+".png", "frames": count, "fps": 4,
            "duration_seconds": count/4, "decode_verified": True,
            "layout": "2 rows (cyclone, anticyclone) × 3 heights", "heights_m": HEIGHTS.tolist(),
            "field": "perturbation vertical vorticity" if is_vort else "buoyancy perturbation",
            "units": units, "passive_tracer": False,
            "scale": {"type": "asinh", "width": width, "limits": [-limit, limit],
                      "fixed_across_all_frames_cases_and_depths": True, "clipping": False},
            "decoded_size": list(metadata["size"])}


def write_derived(cases, out):
    path = out/"derived_scalars.csv"
    with path.open("w") as file:
        header = ["case", "selected_modes", "gradient_removed", "seed", "vortex_sign", "f_sign",
                  "t_fabs", "time_hours", "kinetic", "amplification", "virtual_bt_power",
                  "bt_power_over_K", "bt_power_over_K0", "integrated_virtual_bt_work",
                  "integrated_bt_work_over_K0", "buoyancy_power", "viscous_power",
                  "integrated_kinetic_residual", "integrated_kinetic_residual_over_K0",
                  "integrated_thermal_variance_residual", "perturbation_rms_over_base"]
        writer = csv.writer(file)
        writer.writerow(header)
        for key, c in sorted(cases.items()):
            m, control, seed, sign, f = key
            for a in c["data"]:
                k0 = c["k0"]
                writer.writerow([c["name"], m, int(control), seed, sign, f, a["t_fabs"], a["time_hours"],
                                 a["kinetic"], a["kinetic"]/k0, a["bt_power"], a["bt_power"]/a["kinetic"], a["bt_power"]/k0,
                                 a["integral_bt_power"], a["integral_bt_power"]/k0, a["buoyancy_power"], a["viscous_power"],
                                 a["integrated_kinetic_residual"], a["integrated_kinetic_residual"]/k0,
                                 a["integrated_thermal_variance_residual"], a["perturbation_rms_over_base"]])
    with gzip.open(out/"derived_scalars.csv.gz", "wb") as target:
        target.write(path.read_bytes())


def run(root, out):
    out.mkdir(parents=True, exist_ok=False)
    campaign, cases = read_cases(root)
    summaries = {case["name"]: summarize(case) for case in cases.values()}
    comparison, flags = comparisons(cases, summaries)
    reflection = reflection_checks(cases)
    if not all(r["numerically_consistent"] for r in reflection.values()):
        flags.append("Full-run reflected-f consistency requires numerical review.")
    max_kinetic_residual = max(s["kinetic_residual_max_over_work_scale"] for s in summaries.values())
    max_thermal_residual = max(s["thermal_variance_residual_max_over_work_scale"] for s in summaries.values())
    if max(max_kinetic_residual, max_thermal_residual) > 1e-3:
        flags.append("Integrated work residual exceeds 0.1% of the declared work scale; review before interpreting the response.")
    plots(cases, out)
    write_derived(cases, out)
    movies = [render(cases, out, i) for i in (0, 1)]
    report = {"status": "ANALYZED_MEDIA_VERIFIED", "compute_location": "Unity",
              "interpretation_status": "REVIEW_REQUIRED" if flags else "NO_QUALITATIVE_CONTRADICTION_IN_TESTED_PAIRS",
              "case_count": len(cases), "time_range_fabs_t": [0., 10.], "time_hours": [0., 10/F_ABS/3600],
              "late_window_fabs_t": [5., 10.], "seed_count": 2, "mode_counts": [1, 2],
              "source_code_sha256": campaign["code_sha256"], "basis_sha256": campaign["basis_sha256"],
              "cases": summaries, "paired_comparisons": comparison, "reflection_checks": reflection,
              "interpretation_flags": flags, "movies": movies,
              "plots": ["response_work_modes1.png", "response_work_modes2.png", "response_amplification_and_residual.png"],
              "derived_scalars": "derived_scalars.csv.gz",
              "normalization": "Dimensionless H and |f| units. PBT/K and PBT/K0 are rates per unit |f|t; buoyancy movies convert to m/s².",
              "work_semantics": "Integral PBT is prospective barotropic work from actual RK-stage accumulation. The imposed vortex never evolves.",
              "scalar_interval_fabs_by_case": {c["name"]: c["config"]["scalar_interval_fabs"] for c in cases.values()},
              "late_window_semantics": "PBT and PBT/K0 means use differences of RK-stage work integrals; PBT/K uses trapezoidal scalar samples at each case's recorded cadence.",
              "limitations": ["Linear frozen-profile Boussinesq response; no saturated convection or freely forming LSV is simulated.",
                              "Perturbation growth can occur while stresses extract energy from the prescribed vortex.",
                              "Gradient removal changes the frozen stratification; it is not instantaneous bottom-heating switch-off.",
                              "Two seeds and one horizontal resolution do not establish statistical or spatial convergence.",
                              "One/two vertical-mode comparisons test sensitivity within the selected subspaces, not general vertical convergence.",
                              "Numerical accounting tolerances are not physical sampling uncertainty.",
                              "Original atmospheric density variation is omitted by the Boussinesq approximation.",
                              "Movies reconstruct actual truncated-model fields at three heights; the scalar is buoyancy, not passive tracer.",
                              "Selected modes can be concentrated below the displayed heights. Weak upper-plane signals are retained with the same fixed scale; they do not imply a rendering failure."]}
    (out/"report.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    assets = {p.name: {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
              for p in out.iterdir() if p.is_file()}
    (out/"transport.json").write_text(json.dumps({"source_root": str(root), "assets": assets,
        "analyzer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "interpretation_status": report["interpretation_status"],
                      "flags": flags, "movies": [m["file"] for m in movies]}), flush=True)


if __name__ == "__main__":
    assert len(sys.argv) == 3, __doc__
    run(Path(sys.argv[1]), Path(sys.argv[2]))
