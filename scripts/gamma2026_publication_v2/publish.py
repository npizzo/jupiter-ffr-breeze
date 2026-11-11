#!/usr/bin/env python3
"""Analyze and render completed gamma-plane pairs on an allocated Unity CPU node.

This publication driver never runs or changes a model. Native saved arrays stay
on Unity. Output contains only movies, figures, compact diagnostics and provenance.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tomllib

import h5py
import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm
import numpy as np

DAY = 86400.0
EXPECTED_PAIRS = ((1024, 11), (1024, 22), (1024, 33), (2048, 11))
KINDS = ("polar", "flat")
LABELS = ("Polar gamma-plane", "Flat-PV trap control")
CROSSCHECK_RTOL = 5e-6
CROSSCHECK_ATOL = 1e-18
SHARED_KEYS = ("seed", "n", "gamma_m2_s", "rms_velocity_m_s", "domain_Lgamma",
               "trap_Lgamma", "initial_inverse_wavenumber_m", "annulus_fraction",
               "taper", "stepper", "aliased_fraction", "cfl", "max_dt_seconds",
               "end_days", "output_days", "movie_grid", "source_hash",
               "geophysicalflows_commit")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_frame(path, expected_n):
    with h5py.File(path, "r") as d:
        result = {key: np.asarray(d[key], dtype=np.float64) for key in
                  ("x_m", "y_m", "vorticity_s_inv", "pv_s_inv")}
        result["time_seconds"] = float(np.asarray(d["time_seconds"]))
        result["sampling_stride"] = int(np.asarray(d["sampling_stride"]))
    x, y, zeta, pv = (result[key] for key in
                       ("x_m", "y_m", "vorticity_s_inv", "pv_s_inv"))
    # h5py reverses Julia's matrix dimension order: the saved [x,y] array is [y,x].
    if zeta.shape != (len(y), len(x)) or pv.shape != zeta.shape:
        raise ValueError(f"Unexpected saved field dimensions: {path}")
    if zeta.shape != (expected_n, expected_n) or result["sampling_stride"] != 1:
        raise ValueError("Native-grid diagnostics require unreduced snapshots")
    if not all(np.isfinite(result[key]).all() for key in ("x_m", "y_m", "vorticity_s_inv", "pv_s_inv")):
        raise ValueError(f"Nonfinite saved field: {path}")
    if not math.isfinite(result["time_seconds"]):
        raise ValueError(f"Nonfinite time: {path}")
    for coordinate in (x, y):
        spacing = np.diff(coordinate)
        if not (np.all(spacing > 0) and np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-6)):
            raise ValueError("Uniform, increasing coordinates required")
    return result


def inventory(root):
    result = []
    for n, seed in EXPECTED_PAIRS:
        for kind in KINDS:
            name = f"{kind}_n{n}_seed{seed}"
            meta = root / name / "metadata.toml"
            item = {"case": name, "grid": n, "seed": seed, "kind": kind, "status": "NOT_FOUND"}
            if meta.exists():
                try:
                    m = tomllib.loads(meta.read_text())
                    item.update(status=m["status"], completed_days=m.get("completed_days"),
                                requested_days=m["end_days"])
                except (OSError, tomllib.TOMLDecodeError, KeyError):
                    item["status"] = "METADATA_IN_TRANSITION"
            result.append(item)
    return result


def case_ready(root, name):
    path = root / name / "metadata.toml"
    if not path.is_file():
        return False
    try:
        m = tomllib.loads(path.read_text())
        return m["status"] == "COMPLETE" and abs(m["completed_days"] - m["end_days"]) < 1e-7
    except (OSError, tomllib.TOMLDecodeError, KeyError):
        # A still-running case may be in its short, non-atomic metadata rewrite.
        return False


def case_analysis(root, name, out):
    folder = root / name
    metadata_path = folder / "metadata.toml"
    m = tomllib.loads(metadata_path.read_text())
    if not case_ready(root, name):
        raise ValueError(f"Case is not complete: {name}")
    if m["end_days"] != 75 or m["output_days"] != 1 or m["movie_grid"] != m["n"]:
        raise ValueError("This publication driver requires the frozen 75-day daily native-grid campaign")
    files = sorted((folder / "frames").glob("frame_*.jld2"))
    expected_frames = round(m["end_days"] / m["output_days"]) + 1
    if len(files) != expected_frames:
        raise ValueError(f"Expected {expected_frames} daily frames in {name}; got {len(files)}")
    diag = np.atleast_1d(np.genfromtxt(folder / "diagnostics.csv", names=True, delimiter=","))
    if len(diag) != len(files):
        raise ValueError("Diagnostic rows and saved times disagree")
    if not all(np.isfinite(diag[key]).all() for key in diag.dtype.names):
        raise ValueError("Nonfinite simulation diagnostic")
    s = (m["rms_velocity_m_s"] / m["gamma_m2_s"]) ** (1 / 3)
    initial = read_frame(files[0], m["n"])
    x, y = initial["x_m"], initial["y_m"]
    xx, yy = np.meshgrid(x, y)
    rr = np.hypot(xx, yy)
    trap = rr < m["trap_Lgamma"] * s
    core = rr < s
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    kx = 2 * np.pi * np.fft.fftfreq(len(x), dx)
    ky = 2 * np.pi * np.fft.fftfreq(len(y), dy)
    k2 = ky[:, None] ** 2 + kx[None, :] ** 2
    inverse_k2 = np.divide(1.0, k2, out=np.zeros_like(k2), where=k2 > 0)
    initial_pv_background = initial["pv_s_inv"] - initial["vorticity_s_inv"]
    radial_edges = np.linspace(0, m["trap_Lgamma"] * s, 51)
    radial_bin = np.searchsorted(radial_edges, rr, side="right") - 1
    radial_valid = trap & (radial_bin >= 0) & (radial_bin < 50)
    radial_count = np.bincount(radial_bin[radial_valid], minlength=50)
    rows, frames = [], []
    max_vorticity = 0.0
    for index, path in enumerate(files):
        f = initial if index == 0 else read_frame(path, m["n"])
        if not (np.array_equal(f["x_m"], x) and np.array_equal(f["y_m"], y)):
            raise ValueError("Coordinates changed during a run")
        time_days = f["time_seconds"] / DAY
        if abs(time_days - float(diag["time_days"][index])) > 1e-7 or abs(time_days - index * m["output_days"]) > 1e-7:
            raise ValueError("Saved and diagnostic times disagree")
        zeta, pv = f["vorticity_s_inv"], f["pv_s_inv"]
        maximum = float(np.max(np.abs(zeta)))
        max_vorticity = max(max_vorticity, maximum)
        # Inf deformation radius, no imposed mean flow: E = 1/2 <|grad psi|^2>.
        zhat = np.fft.fft2(zeta)
        energy = 0.5 * float(np.sum(np.abs(zhat) ** 2 * inverse_k2)) / zeta.size ** 2
        relative_z = 0.5 * float(np.mean(zeta ** 2))
        pv_z = 0.5 * float(np.mean(pv ** 2))
        observed = (energy, relative_z, pv_z)
        original = tuple(float(diag[key][index]) for key in
                         ("energy_m2_s2", "relative_enstrophy_s2", "pv_enstrophy_s2"))
        # Preserve the original comparison thresholds. Relative-vorticity movies
        # and their derived statistics are gated by energy and relative enstrophy.
        # Total-PV disagreement is a separately exposed audit failure, never a
        # successful check or an automatically attributed Float32 rounding error.
        checked = np.isclose(observed, original, rtol=CROSSCHECK_RTOL, atol=CROSSCHECK_ATOL)
        if not checked[:2].all():
            raise ValueError(f"Native vorticity/solver diagnostic disagreement in {name}, day {time_days}: {observed} vs {original}")
        row = {key: float(diag[key][index]) for key in diag.dtype.names}
        row.update(snapshot_energy_m2_s2=energy, snapshot_relative_enstrophy_s2=relative_z,
                   snapshot_pv_enstrophy_s2=pv_z,
                   snapshot_pv_enstrophy_fraction=pv_z / float(diag["pv_enstrophy_s2"][0]),
                   energy_crosscheck_passed=int(checked[0]),
                   relative_enstrophy_crosscheck_passed=int(checked[1]),
                   pv_enstrophy_crosscheck_passed=int(checked[2]),
                   relative_enstrophy_fraction=original[1] / float(diag["relative_enstrophy_s2"][0]),
                   pv_enstrophy_fraction=original[2] / float(diag["pv_enstrophy_s2"][0]),
                   snapshot_energy_relative_error=abs(energy / original[0] - 1),
                   snapshot_relative_enstrophy_relative_error=abs(relative_z / original[1] - 1),
                   snapshot_pv_enstrophy_relative_error=abs(pv_z / original[2] - 1),
                   signed_vorticity_area_integral_m2_s=float(np.sum(zeta) * dx * dy),
                   outside_trap_relative_enstrophy_fraction=float(np.sum(zeta[~trap] ** 2) / np.sum(zeta ** 2)),
                   static_background_max_abs_error_s_inv=float(np.max(np.abs(pv - zeta - initial_pv_background))))
        for label, sign in (("positive", 1), ("negative", -1)):
            weight = np.maximum(sign * zeta, 0) * trap
            total = float(np.sum(weight))
            if not total > 0:
                raise ValueError("Both vorticity signs are required for radial separation")
            cx, cy = float(np.sum(weight * xx) / total), float(np.sum(weight * yy) / total)
            row.update({f"{label}_circulation_in_trap_m2_s": total * dx * dy,
                        f"{label}_mean_radius_Lgamma": float(np.sum(weight * rr) / total / s),
                        f"{label}_rms_radius_Lgamma": float(np.sqrt(np.sum(weight * rr ** 2) / total) / s),
                        f"{label}_centroid_x_Lgamma": cx / s,
                        f"{label}_centroid_y_Lgamma": cy / s,
                        f"{label}_centroid_offset_Lgamma": math.hypot(cx, cy) / s,
                        f"{label}_core_circulation_fraction": float(np.sum(weight[core]) / total)})
        row["radial_separation_Lgamma"] = row["negative_mean_radius_Lgamma"] - row["positive_mean_radius_Lgamma"]
        initial_separation = rows[0]["radial_separation_Lgamma"] if rows else row["radial_separation_Lgamma"]
        row["radial_separation_change_Lgamma"] = row["radial_separation_Lgamma"] - initial_separation
        rows.append(row)
        frames.append({"file": path.name, "time_days": time_days, "sha256": digest(path),
                       "sampling_stride": f["sampling_stride"]})
        if index == len(files) - 1:
            radial_sum = np.bincount(radial_bin[radial_valid], weights=zeta[radial_valid], minlength=50)
            radial_mean = np.divide(radial_sum, radial_count, out=np.zeros(50), where=radial_count > 0)
    csv_path = out / f"{name}_diagnostics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    final = rows[-1]
    summary = {"case": name, "status": "COMPLETE", "grid": m["n"], "kind": m["kind"], "seed": m["seed"],
               "completed_days": final["time_days"], "saved_frames": len(files), "Lgamma_m": s,
               "simulation_wall_seconds": m["wall_seconds"], "metadata_sha256": digest(metadata_path),
               "simulation_diagnostics_sha256": digest(folder / "diagnostics.csv"),
               "source_hash": m["source_hash"], "config_hash": m["config_hash"],
               "geophysicalflows_commit": m["geophysicalflows_commit"],
               "final": final, "initial_radial_separation_Lgamma": rows[0]["radial_separation_Lgamma"],
               "final_radial_profile": {"r_over_Lgamma": ((radial_edges[:-1] + radial_edges[1:]) / (2 * s)).tolist(),
                                         "mean_vorticity_s_inv": radial_mean.tolist()},
               "snapshot_crosschecks": {"rtol": CROSSCHECK_RTOL, "atol": CROSSCHECK_ATOL,
                   "all_passed": all(r["pv_enstrophy_crosscheck_passed"] for r in rows),
                   "vorticity_energy_and_relative_enstrophy_all_passed": True,
                   "pv_enstrophy_all_passed": all(r["pv_enstrophy_crosscheck_passed"] for r in rows),
                   "pv_failed_days": [r["time_days"] for r in rows if not r["pv_enstrophy_crosscheck_passed"]],
                   "pv_audit_interpretation": "Observed total-PV discrepancy is unresolved; it is not presumed to be Float32 roundoff. Original tolerance retained. Relative-vorticity movies and derived diagnostics pass their independent energy and relative-enstrophy checks.",
                   "max_energy_relative_error": max(r["snapshot_energy_relative_error"] for r in rows),
                   "max_relative_enstrophy_relative_error": max(r["snapshot_relative_enstrophy_relative_error"] for r in rows),
                   "max_pv_enstrophy_relative_error": max(r["snapshot_pv_enstrophy_relative_error"] for r in rows)},
               "frames": frames}
    dump(out / f"{name}_summary.json", summary)
    shutil.copyfile(metadata_path, out / f"{name}_metadata.toml")
    return {"name": name, "metadata": m, "files": files, "initial": initial, "rows": rows,
            "summary": summary, "maximum": max_vorticity, "Lgamma": s}


def pair_check(cases):
    a, b = cases
    for key in SHARED_KEYS:
        if a["metadata"][key] != b["metadata"][key]:
            raise ValueError(f"Unmatched polar/control parameter: {key}")
    if not np.array_equal(a["initial"]["vorticity_s_inv"], b["initial"]["vorticity_s_inv"]):
        raise ValueError("Polar and flat control initial vorticity fields differ")
    if [r["time_days"] for r in a["rows"]] != [r["time_days"] for r in b["rows"]]:
        raise ValueError("Pair requires identical saved times")


def field_artist(ax, f, case, norm, label):
    s = case["Lgamma"]
    x, y = f["x_m"] / s, f["y_m"] / s
    dx, dy = x[1] - x[0], y[1] - y[0]
    im = ax.imshow(f["vorticity_s_inv"] / 1e-4, origin="lower", interpolation="nearest",
                   extent=(x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2),
                   cmap="RdBu_r", norm=norm)
    ax.add_patch(plt.Circle((0, 0), case["metadata"]["trap_Lgamma"], fill=False, color="0.3", lw=.6, ls="--"))
    ax.set(xlabel="x / Lγ", ylabel="y / Lγ", title=label, aspect="equal")
    return im


def verify_video(path, expected_frames, fps):
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    subprocess.run([ffmpeg, "-v", "error", "-i", str(path), "-f", "null", "-"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    frames, duration = imageio_ffmpeg.count_frames_and_secs(str(path))
    if frames != expected_frames or abs(duration - expected_frames / fps) > .05:
        raise ValueError(f"Video duration/frame count mismatch: {path}")
    if path.stat().st_size >= 24_000_000:
        raise ValueError("Movie exceeds the publication size limit; rerender with higher CRF")
    return {"frames": frames, "duration_seconds": duration, "fps": fps,
            "bytes": path.stat().st_size, "sha256": digest(path), "full_decode_passed": True}


def render_pair(cases, out, maximum, scale_cases, fps, crf):
    pair_check(cases)
    n, seed = cases[0]["metadata"]["n"], cases[0]["metadata"]["seed"]
    stem = f"gamma_flat_n{n}_seed{seed}_0_75days"
    movie = out / (stem + ".mp4")
    # A physical reference fixes the linear core of the asinh transformation.
    width = .1 * cases[0]["metadata"]["rms_velocity_m_s"] / cases[0]["metadata"]["initial_inverse_wavenumber_m"]
    limit = max(maximum, width)
    norm = AsinhNorm(linear_width=width / 1e-4, vmin=-limit / 1e-4, vmax=limit / 1e-4, clip=False)
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 7.2), layout="constrained")
    images = [field_artist(ax, c["initial"], c, norm, label) for ax, c, label in zip(axes, cases, LABELS)]
    cb = fig.colorbar(images[0], ax=axes, shrink=.82, pad=.02)
    cb.set_label("Relative vorticity ζ [10⁻⁴ s⁻¹] · fixed asinh scale")
    title = fig.suptitle(f"2D barotropic comparison · {n}² · seed {seed} · Earth day 0", fontsize=13)
    fig.supxlabel("No vertical levels or tracer are represented. Dashed circle: trap boundary. Lγ = %.0f km." % (cases[0]["Lgamma"] / 1000), fontsize=10)
    fig.canvas.draw(); fig.set_layout_engine(None)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-crf", str(crf), "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    with writer.saving(fig, str(movie), dpi=100):
        for i, row in enumerate(cases[0]["rows"]):
            for image, c in zip(images, cases):
                image.set_data(read_frame(c["files"][i], n)["vorticity_s_inv"] / 1e-4)
            title.set_text(f"2D barotropic comparison · {n}² · seed {seed} · Earth day {row['time_days']:.0f}")
            writer.grab_frame()
    fig.savefig(out / (stem + "_poster.png"), dpi=100)
    plt.close(fig)
    media = verify_video(movie, len(cases[0]["rows"]), fps)
    media.update(file=movie.name, cases=[c["name"] for c in cases], dimension=2,
                 layout="1 row × 2 matched cases; no physical depths exist",
                 every_saved_time_included=True, temporal_interpolation=False, field_clipping=False,
                 time_days=[r["time_days"] for r in cases[0]["rows"]],
                 normalization={"type": "asinh", "linear_width_s_inv": width,
                     "limits_s_inv": [-limit, limit], "fixed_in_time": True,
                     "shared_between_cases": True, "scale_reference_cases": scale_cases,
                     "definition": "Full vorticity extrema across all selected completed pairs; physical linear width = 0.1 U/Linit"})
    dump(out / (stem + "_video.json"), media)
    # Contact sheet exposes early, middle and final frames for human visual QA.
    indices = [0, len(cases[0]["files"]) // 2, len(cases[0]["files"]) - 1]
    fig, axes = plt.subplots(3, 2, figsize=(10, 14), layout="constrained")
    for row_index, frame_index in enumerate(indices):
        for col, c in enumerate(cases):
            f = read_frame(c["files"][frame_index], n)
            im = field_artist(axes[row_index, col], f, c, norm,
                              f"{LABELS[col]} · day {f['time_seconds'] / DAY:.0f}")
    fig.colorbar(im, ax=axes, shrink=.7, label="Relative vorticity ζ [10⁻⁴ s⁻¹] · fixed asinh")
    fig.suptitle(f"2D barotropic comparison · {n}² · seed {seed}")
    fig.savefig(out / (stem + "_contact.png"), dpi=100); plt.close(fig)
    return media


def diagnostic_figures(cases, out):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for case_index, c in enumerate(cases):
        rows = c["rows"]
        days = [r["time_days"] for r in rows]
        label = f"{c['metadata']['kind']} · {c['metadata']['n']}² · seed {c['metadata']['seed']}"
        for ax, key in zip(axes.flat, ("energy_fraction", "pv_enstrophy_fraction", "relative_enstrophy_fraction", "radial_separation_Lgamma")):
            ax.plot(days, [r[key] for r in rows], label=label, lw=1.2, color=f"C{case_index % 10}")
            ax.set_xlabel("Earth day"); ax.grid(alpha=.2)
        # Show the independently reconstructed PV series even when its audit
        # disagrees with the solver record; do not replace either series.
        color = axes[0, 1].lines[-1].get_color()
        axes[0, 1].plot(days, [r["snapshot_pv_enstrophy_fraction"] for r in rows],
                        color=color, ls="--", lw=.8)
    for ax, title in zip(axes.flat, ("Kinetic energy / initial", "Total-PV: solver solid / saved fields dashed", "Relative-vorticity enstrophy / initial", "Mean radius: negative − positive vorticity, / Lγ")):
        ax.set_title(title, fontsize=11)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Unforced 2D barotropic pilots · numerical budgets and descriptive sign separation")
    pv_failures = any(not c["summary"]["snapshot_crosschecks"]["pv_enstrophy_all_passed"] for c in cases)
    note = "PV audit: mismatch exceeds 5×10⁻⁶; cause unresolved. " if pv_failures else "PV audit: saved-field comparisons pass 5×10⁻⁶. "
    fig.supxlabel(note + "Relative enstrophy need not be conserved. Radial separation is not evidence of an inverse cascade.", fontsize=8)
    fig.savefig(out / "gamma_diagnostics.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for c in cases:
        profile = c["summary"]["final_radial_profile"]
        ax.plot(profile["r_over_Lgamma"], np.asarray(profile["mean_vorticity_s_inv"]) / 1e-4,
                label=f"{c['metadata']['kind']} {c['metadata']['n']}² s{c['metadata']['seed']}")
    ax.set(xlabel="Radius / Lγ", ylabel="Annular-mean relative vorticity [10⁻⁴ s⁻¹]",
           title="Day-75 radial profiles inside the trap")
    ax.axhline(0, color="k", lw=.5); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.savefig(out / "gamma_radial_profiles.png", dpi=130); plt.close(fig)


def write_report(out, cases, media, inventory_rows, requested_pairs):
    pv_failures = [c for c in cases if not c["summary"]["snapshot_crosschecks"]["pv_enstrophy_all_passed"]]
    table = ["| Case | Grid | Days | E/E₀ | Solver PV enstrophy / initial | Relative enstrophy / initial | Mean radius (negative − positive) / Lγ |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for c in cases:
        f = c["summary"]["final"]
        table.append(f"| {c['name']} | {c['metadata']['n']}² | {f['time_days']:.0f} | {f['energy_fraction']:.6f} | {f['pv_enstrophy_fraction']:.6f} | {f['relative_enstrophy_fraction']:.6f} | {f['radial_separation_Lgamma']:.5f} |")
    included = {c["name"] for c in cases}
    missing = [r for r in inventory_rows if r["case"] not in included]
    explanation = "\n".join(f"- {r['case']}: {r['status']} (not included in this publication)." for r in missing) or "All eight planned cases are included."
    text = """# Polar gamma-plane and flat-trap pilots: completed pairs

These are **two-dimensional, unforced barotropic QG** simulations, with no convection, vertical levels or passive tracer. They independently implement the reference supplementary-information initialization: a ±5% spectral annulus, quartic radial taper and common initial RMS speed. This is a 75-Earth-day early-organization test, not a mature vortex-crystal reproduction. The current solver/filter, adaptive timestep and new seeds differ from the historical calculation.

## Movies

"""
    text += "\n".join(f"- [{m['file']}]({m['file']}): synchronized polar/control comparison, all {m['frames']} saved times, two actual 2D fields." for m in media)
    text += "\n\nMovies use a symmetric asinh color map with physical vorticity labels, fixed across all times and selected cases. The range includes all selected saved extrema. No field clipping, duplicated depths or temporal interpolation is used. Each MP4 was fully decoded on Unity and its frame count and duration verified.\n\n"
    if pv_failures:
        text += "**Total-PV diagnostic audit has unresolved discrepancies.** The original relative tolerance of 5×10⁻⁶ was retained and failed at the times listed below. This publication does not certify all budget checks or make a tight total-PV conservation claim. Independently reconstructed kinetic energy and relative-vorticity enstrophy pass that same threshold at every saved time, supporting publication of the actual vorticity fields and descriptive vorticity statistics. The total-PV issue is not assumed to be Float32 roundoff.\n\n"
    else:
        text += "All saved-field energy, relative-enstrophy and total-PV comparisons passed the original 5×10⁻⁶ relative tolerance for these selected cases. This is diagnostic consistency, not proof of physical conservation.\n\n"
    text += "| Case | Maximum relative PV mismatch | Failed saved times / total | Failed Earth days |\n|---|---:|---:|---|\n"
    for c in cases:
        check = c["summary"]["snapshot_crosschecks"]
        failed = check["pv_failed_days"]
        days = ", ".join(f"{d:.0f}" for d in failed) or "none"
        text += f"| {c['name']} | {check['max_pv_enstrophy_relative_error']:.8g} | {len(failed)} / {len(c['rows'])} | {days} |\n"
    text += "\n## Numerical and descriptive diagnostics\n\n"
    text += "\n".join(table)
    pair_effects = []
    for n, seed in EXPECTED_PAIRS:
        pair = [next((c for c in cases if c["name"] == f"{kind}_n{n}_seed{seed}"), None) for kind in KINDS]
        if all(pair):
            polar_change, flat_change = (c["rows"][-1]["radial_separation_change_Lgamma"] for c in pair)
            pair_effects.append({"grid": n, "seed": seed, "polar_separation_change_Lgamma": polar_change,
                                "flat_separation_change_Lgamma": flat_change,
                                "polar_minus_flat_change_Lgamma": polar_change - flat_change})
    text += "\n\n| Grid / seed | Polar change from day 0 | Flat-control change | Polar minus flat change |\n|---|---:|---:|---:|\n"
    for effect in pair_effects:
        text += f"| {effect['grid']}² / {effect['seed']} | {effect['polar_separation_change_Lgamma']:.5f} | {effect['flat_separation_change_Lgamma']:.5f} | {effect['polar_minus_flat_change_Lgamma']:.5f} |\n"
    text += "\nAll separation changes are in Lγ. Subtracting day 0 removes the finite-realization initial sign imbalance; subtracting the paired flat-control change helps isolate the effect of the polar PV geometry. It is a descriptive matched-pair effect, not a significance test.\n"
    seed_summaries = []
    for n in sorted({e["grid"] for e in pair_effects}):
        entries = [e for e in pair_effects if e["grid"] == n]
        values = np.array([e["polar_minus_flat_change_Lgamma"] for e in entries])
        summary = {"grid": n, "seeds": [e["seed"] for e in entries], "pairs": len(entries),
                   "mean_polar_minus_flat_change_Lgamma": float(np.mean(values)),
                   "range_Lgamma": [float(np.min(values)), float(np.max(values))]}
        if len(values) > 1:
            summary["sample_standard_deviation_Lgamma"] = float(np.std(values, ddof=1))
        seed_summaries.append(summary)
    text += """

Every native Float32 snapshot was independently compared against the solver's Float64 kinetic-energy, relative-enstrophy and total-PV-enstrophy diagnostics (unchanged relative tolerance 5×10⁻⁶, absolute tolerance 10⁻¹⁸). The per-field results above distinguish passed checks from failed PV checks. Kinetic energy was reconstructed by periodic spectral inversion of saved relative vorticity at infinite deformation radius. Relative enstrophy is not an invariant when the background PV varies in space; its decline is not, by itself, numerical dissipation. The recorded total-PV diagnostic is ½⟨(ζ+η)²⟩. The PV panel shows both solver and saved-field reconstructions against the same solver initial denominator; CSVs preserve both series and their discrepancy at each time. Finite-grid filtering can remove energy and PV enstrophy, but this observation does not explain the solver/snapshot mismatch.

The radial-separation statistic is the mean radius of negative vorticity minus that of positive vorticity, each weighted by its sign's vorticity magnitude **inside r < 5 Lγ**. Positive values indicate that positive vorticity is more centrally weighted. CSVs also contain sign-weighted Cartesian centroids, RMS radii, core fractions and the initial separation. These quantities describe an Eulerian field; they do not identify or track individual vortices, and their values should be compared against the matching flat-trap control and seed variation. No inverse-cascade claim is made: this publication does not calculate nonlinear spectral transfer.

Resolution differences between the 1024² and 2048² matched-seed runs are sensitivity evidence, not a convergence proof. Day-75 morphology alone cannot establish a mature polar crystal. Comparisons with convection require a separate, dynamically matched 3D experiment.

## Campaign coverage at rendering

"""
    text += explanation + "\n\n## Provenance\n\nAll analysis, rendering and media validation ran in an allocated Unity CPU job. Raw model fields remain on Unity. Per-case metadata, native diagnostic CSVs, frame hashes and video normalization records accompany this report. The JSON manifest records the exact renderer and source revisions.\n"
    (out / "README.md").write_text(text)
    refinements = []
    for kind in KINDS:
        pair = [next((c for c in cases if c["name"] == f"{kind}_n{n}_seed11"), None) for n in (1024, 2048)]
        if all(pair):
            keys = ("energy_fraction", "pv_enstrophy_fraction", "relative_enstrophy_fraction", "radial_separation_Lgamma")
            refinements.append({"kind": kind, "seed": 11,
                                "difference_2048_minus_1024": {key: pair[1]["rows"][-1][key] - pair[0]["rows"][-1][key] for key in keys},
                                "interpretation": "Matched-seed grid sensitivity only; no convergence claim"})
    return {"status": "COMPLETE_WITH_PV_DIAGNOSTIC_DISCREPANCY" if pv_failures else "COMPLETE",
            "compute_location": "Unity allocated CPU job", "slurm_job_id": os.environ["SLURM_JOB_ID"],
            "output_dimension": 2, "raw_fields_downloaded": False, "requested_pairs": requested_pairs,
            "cases_included": [c["name"] for c in cases], "campaign_inventory_at_rendering": inventory_rows,
            "cases_not_included": missing, "media": media, "grid_sensitivity": refinements,
            "paired_radial_separation_changes": pair_effects, "seed_summaries": seed_summaries,
            "snapshot_crosscheck_passed": not pv_failures,
            "vorticity_energy_and_relative_enstrophy_checks_passed": True,
            "pv_snapshot_crosscheck_passed": not pv_failures,
            "pv_audit_failures": {c["name"]: c["summary"]["snapshot_crosschecks"] for c in pv_failures},
            "publication_scope": "Verified actual vorticity movies and descriptive vorticity diagnostics; total-PV consistency assessed separately, with unresolved failures retained",
            "renderer_sha256": digest(__file__),
            "limitations": ["No convection, physical depth samples or tracer in this model",
                            "Any failed total-PV snapshot/solver comparison is unresolved and disallows a tight total-PV conservation claim",
                            "75 days is an early-organization pilot, not a mature-crystal test",
                            "Sign-weighted field centroids are not tracked vortices",
                            "No nonlinear spectral-transfer diagnostic; no inverse-cascade inference",
                            "One refined seed pair cannot establish ensemble convergence"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", nargs="+", default=[f"{n}:{seed}" for n, seed in EXPECTED_PAIRS])
    parser.add_argument("--require-all", action="store_true", help="Fail unless every selected pair is complete")
    parser.add_argument("--fps", type=int, default=6)
    parser.add_argument("--crf", type=int, default=28)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Scientific publication computation must run in an allocated Unity job")
    if args.output.exists():
        raise ValueError("Use a fresh publication directory; existing assets are immutable")
    if args.fps < 1 or not 18 <= args.crf <= 39:
        raise ValueError("Invalid movie encoding parameters")
    requested = [tuple(map(int, pair.split(":"))) for pair in args.pairs]
    if len(set(requested)) != len(requested) or any(pair not in EXPECTED_PAIRS for pair in requested):
        raise ValueError("Unknown or duplicate pair specification")
    available = [pair for pair in requested if all(case_ready(args.root, f"{kind}_n{pair[0]}_seed{pair[1]}") for kind in KINDS)]
    if not available or (args.require_all and available != requested):
        raise ValueError("Required complete polar/flat pairs are not available")
    args.output.mkdir(parents=True)
    plt.rcParams.update({"animation.ffmpeg_path": imageio_ffmpeg.get_ffmpeg_exe(), "font.size": 10})
    inventory_rows = inventory(args.root)
    pairs, cases = [], []
    for n, seed in available:
        pair = [case_analysis(args.root, f"{kind}_n{n}_seed{seed}", args.output) for kind in KINDS]
        pair_check(pair); pairs.append(pair); cases.extend(pair)
    limit = max(c["maximum"] for c in cases)
    scale_cases = [c["name"] for c in cases]
    media = [render_pair(pair, args.output, limit, scale_cases, args.fps, args.crf) for pair in pairs]
    diagnostic_figures(cases, args.output)
    report = write_report(args.output, cases, media, inventory_rows, args.pairs)
    report["rendered_at_utc"] = datetime.now(timezone.utc).isoformat()
    report["publication_runtime"] = {"numpy": np.__version__, "matplotlib": matplotlib.__version__,
                                      "h5py": h5py.__version__, "imageio_ffmpeg": imageio_ffmpeg.__version__}
    report["assets"] = {p.name: {"bytes": p.stat().st_size, "sha256": digest(p)} for p in sorted(args.output.iterdir()) if p.is_file()}
    dump(args.output / "publication.json", report)
    print("READY_FOR_PUBLICATION", args.output, flush=True)


if __name__ == "__main__":
    main()
