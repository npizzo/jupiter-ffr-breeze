"""Completed Ferrel screening reports and actual three-depth movies; Unity only.

This observer never modifies simulation outputs. Native sampled step closure is
reported separately from approximate, centered eddy-stress interpretation.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tomllib
import shutil
from datetime import datetime, timezone

import h5py
import imageio_ffmpeg
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm
from PIL import Image

plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
TERMS = ("advection", "coriolis", "diffusion", "external forcing")
LABELS = {"b0": "Unheated", "b162": "Heated · 162 W/m²",
          "boundary162": "Heated · boundary sensitivity",
          "axis0": "Axisymmetric · unheated", "axis162": "Axisymmetric · heated",
          "fine162": "Heated · 1024² refinement"}
CAVEATS = [
    "Dry anelastic beta-channel transient mechanism screening; no steady or deep Jovian-cell claim.",
    "The first 12 Earth hours are excluded as a conservative startup window, not a demonstrated equilibration timescale. Six-hour pilots and controls lie entirely within that window.",
    "Native momentum closure samples complete SSPRK3 steps; it is not a continuous long-run integrated momentum budget.",
    "Centered covariance stresses, force work and mean-flow transfer are approximate diagnostics, not exact WENO fluxes or a closed kinetic-energy budget.",
    "Positive covariance times mean shear feeds zonal-mean kinetic energy. This is not by itself proof of an inverse cascade.",
    "A single seed, short duration, cap relaxation, sponge and global thermal feedback limit physical attribution.",
    "The boundary case changes both cap relaxation and vertical-velocity damping; it cannot isolate their individual effects.",
    "No residual/TEM circulation is inferred in the nearly neutral layer; mean and eddy tracer fluxes are reported directly.",
    "Psi is overturning transport per unit zonal length in kg/(m s), not a planet-integrated kg/s streamfunction.",
    "Temporal variability and coherence are descriptive; short windows and a single seed do not establish statistical significance.",
]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scalar(d, key):
    value = float(d[key][()])
    if not np.isfinite(value):
        raise ValueError(f"Nonfinite scalar {key}")
    return value


def array(d, key, shape=None):
    a = np.asarray(d[key][()], dtype=float)
    if shape is not None and a.shape != shape:
        raise ValueError(f"{key}: {a.shape} != expected {shape}")
    if not np.isfinite(a).all():
        raise ValueError(f"Nonfinite {key}")
    return a


def time_mean(values, times, start):
    selected = times >= start - 1e-7
    t = times[selected]
    a = np.asarray(values)[selected]
    if len(t) < 2 or t[-1] <= t[0]:
        raise ValueError("At least two saved times required in averaging window")
    # np.trapz works in both current and older Unity NumPy environments.
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return trap(a, t, axis=0) / (t[-1] - t[0]), [float(t[0]), float(t[-1])]


def gate(path):
    g = tomllib.loads(path.read_text())
    if g["status"] != "PASS":
        raise ValueError(f"Validation gate has not passed: {path}")
    root = Path.cwd()
    for name, expected in g["source_hashes"].items():
        if digest(root / "scripts/ferrel2026_v3" / name) != expected:
            raise ValueError(f"Scientific source differs from gate: {name}")
    for name, expected in g["configuration_hashes"].items():
        if digest(root / "configs/ferrel2026" / f"{name}.toml") != expected:
            raise ValueError(f"Configuration differs from gate: {name}")
    return {"file": str(path), "sha256": digest(path), "status": "PASS"}


def collect(name, folder, out):
    m = tomllib.loads((folder / "metadata.toml").read_text())
    if m["status"] != "COMPLETE":
        raise ValueError(f"Refusing incomplete output: {folder}")
    if m["meridional_boundary"] != "channel" or m["coriolis_model"] != "beta-plane":
        raise ValueError("This report is only for the Ferrel beta-channel campaign")
    if not np.isfinite(m["f"]) or m["f"]==0:
        raise ValueError("Finite nonzero f is required for vorticity normalization")
    nx, ny, nz = (m[k] for k in ("nx", "ny", "nz"))
    dy, dz = m["ly"] / ny, m["height"] / nz
    files = sorted(folder.glob("frame_*.jld2"))
    times, rows, frames = [], [], []
    limits = np.zeros(3)
    tracer_range = [0.0, 1.0]
    maps = {k: [] for k in ("mass_streamfunction", "zonal_u", "native_rhou_rhs",
                            "centered_native_advection_remainder", "eddy_to_mean_ke_transfer_y",
                            "eddy_to_mean_ke_transfer_z", "mean_rho_vtracer", "mean_rho_wtracer",
                            "zonal_rho_vtracer", "zonal_rho_wtracer")}
    shapes = {k: (nz, ny) for k in maps}
    shapes.update(mass_streamfunction=(nz + 1, ny + 1), native_rhou_rhs=(4, nz, ny))
    for file in files:
        with h5py.File(file) as d:
            t = scalar(d, "t")
            heights = array(d, "video_heights", (3,))
            if not np.array_equal(heights, m["video_heights_actual_m"]):
                raise ValueError("Saved video heights changed")
            zeta = array(d, "vorticity_slices", (3, ny, nx)) / m["f"]
            tracer = array(d, "tracer_slices", (3, ny, nx))
            limits = np.maximum(limits, np.max(np.abs(zeta), axis=(1, 2)))
            tracer_range = [min(tracer_range[0], float(tracer.min())),
                            max(tracer_range[1], float(tracer.max()))]
            values = {k: array(d, k, shapes[k]) for k in maps}
            rho = array(d, "density", (nz,))
            rv = array(d, "native_rhov", (nz, ny + 1))
            rw = array(d, "native_rhow", (nz + 1, ny))
            psi = values["mass_streamfunction"]
            if any(d[key].dtype != np.dtype("float64") for key in ("native_rhov", "native_rhow", "mass_streamfunction")):
                raise ValueError("Native mass-flux roundoff check requires saved Float64 fields")
            expected = np.vstack((np.zeros((1, ny + 1)), -np.cumsum(rv, axis=0) * dz))
            error = float(np.max(np.abs(psi - expected)))
            # Check the saved native definition without mixing its units with velocity.
            accumulation = dz * float(np.max(np.sum(np.abs(rv), axis=0)))
            bound = 8 * (nz + 2) * np.finfo(float).eps * max(accumulation, float(np.max(np.abs(psi))), np.finfo(float).tiny)
            if error > bound or np.max(np.abs(psi[0])) != 0:
                raise ValueError("Saved mass streamfunction does not integrate native rho*v")
            div = np.diff(rv, axis=1) / dy + np.diff(rw, axis=0) / dz
            reconstruction = np.diff(psi, axis=1) / dy - rw
            scale = max(float(np.max(np.abs(rv))), float(np.max(np.abs(rw))), 1e-10)
            if np.max(np.abs(div)) > 1e-6 * scale / dz + 1e-12:
                raise ValueError("Native mass continuity check failed")
            if np.max(np.abs(reconstruction)) > 1e-6 * scale + 1e-9:
                raise ValueError("Native mass circulation reconstruction failed")
            for axis in ("v", "w"):
                total = array(d, f"total_rho_{axis}tracer", (nz, ny))
                rebuilt = values[f"mean_rho_{axis}tracer"] + values[f"zonal_rho_{axis}tracer"]
                if not np.allclose(total, rebuilt, rtol=1e-10, atol=1e-12):
                    raise ValueError("Mean/eddy tracer flux decomposition failed")
            mean_v = array(d, "zonal_v", (nz, ny))
            mean_w = array(d, "zonal_w", (nz, ny))
            weighted_rms = lambda a: float(np.sqrt(np.sum(rho[:, None] * a*a) / (ny * rho.sum())))
            row = {"time_seconds": t, "earth_day": t / 86400,
                   "max_abs_mass_streamfunction_kg_m_s": float(np.max(np.abs(psi))),
                   "density_weighted_zonal_mean_v_rms_m_s": weighted_rms(mean_v),
                   "density_weighted_zonal_mean_w_rms_m_s": weighted_rms(mean_w),
                   "mass_divergence_max_kg_m3_s": float(np.max(np.abs(div))),
                   "mass_flux_reconstruction_max_kg_m2_s": float(np.max(np.abs(reconstruction))),
                   "psi_native_integration_max_error_kg_m_s": error,
                   "native_advection_rms_kg_m2_s2": float(np.sqrt(np.mean(values["native_rhou_rhs"][0]**2))),
                   "centered_advection_remainder_rms_kg_m2_s2": float(np.sqrt(np.mean(values["centered_native_advection_remainder"]**2)))}
            for key in ("eddy_to_mean_ke_transfer_W_m2", "centered_eddy_force_work_W_m2",
                        "mean_ke_J_m2", "kinetic_energy", "tracer_column", "tracer_min", "tracer_max",
                        "cap_min", "cap_mean", "cap_enthalpy", "budget_correction", "buoyancy_power"):
                row[key] = scalar(d, key)
            for key in ("mean_rho_vtracer", "mean_rho_wtracer", "zonal_rho_vtracer", "zonal_rho_wtracer"):
                row[f"{key}_rms_kg_m2_s"] = float(np.sqrt(np.mean(values[key]**2)))
            if name.startswith("axis"):
                row["axisymmetry_max_velocity_departure_m_s"] = max(
                    float(np.max(np.abs(array(d, key) - np.mean(array(d, key), axis=1, keepdims=True))))
                    for key in ("u", "v", "w"))
                if row["axisymmetry_max_velocity_departure_m_s"] > 1e-8:
                    raise ValueError("Axisymmetric control developed an x-dependent velocity")
                if abs(row["eddy_to_mean_ke_transfer_W_m2"]) > 1e-8:
                    raise ValueError("Axisymmetric control has nonzero eddy transfer")
            for key, value in values.items():
                maps[key].append(value)
        times.append(t); rows.append(row)
        video_hash=hashlib.sha256()
        for values in (zeta,tracer):video_hash.update(np.asarray(values,dtype="<f8",order="C").tobytes())
        frames.append({"file": file.name, "time_seconds": t, "bytes": file.stat().st_size,
                       "video_fields_sha256":video_hash.hexdigest(),
                       "video_hash_encoding":"C-order little-endian Float64: concatenated zeta/f0 then tracer, each shape(3,ny,nx)"})
    times = np.asarray(times)
    expected = np.arange(0, m["stop_seconds"] + m["output_seconds"] / 2, m["output_seconds"])
    if times.shape != expected.shape or not np.allclose(times, expected, rtol=0, atol=1e-6):
        raise ValueError(f"Missing or nonchronological frames in {folder}")
    if len(times) < 2:
        raise ValueError("Need at least two saved frames")
    checks = np.genfromtxt(folder / "native_budget/checks.csv", delimiter=",", names=True)
    checks = np.atleast_1d(checks)
    if len(checks) != m["native_budget_checks"] or not len(checks):
        raise ValueError("Native sampled step count disagrees with metadata")
    if not all(np.isfinite(checks[k]).all() for k in checks.dtype.names):
        raise ValueError("Nonfinite native budget report")
    start = 43200.0 if times[-1] > 43200 else 0.0
    means = {}
    for k, values in maps.items():
        means[k], window = time_mean(values, times, start)
    integral_transfer = float(time_mean([r["eddy_to_mean_ke_transfer_W_m2"] for r in rows], times, start)[0])
    psi_second_moment = time_mean(np.asarray(maps["mass_streamfunction"])**2, times, start)[0]
    psi_std = np.sqrt(np.maximum(0, psi_second_moment-means["mass_streamfunction"]**2))
    psi_abs_mean = time_mean(np.abs(maps["mass_streamfunction"]), times, start)[0]
    coherence = float(np.sum(np.abs(means["mass_streamfunction"]))/max(float(np.sum(psi_abs_mean)),1e-30))
    summary = {"case": name, "label": LABELS[name], "source": str(folder),
               "status": "COMPLETE", "grid": [nx, ny, nz], "heat_flux_W_m2": m["heat_flux"],
               "stop_seconds": m["stop_seconds"], "frames": len(files), "window_seconds": window,
               "window_scope": "12 h startup exclusion applied; equilibration unproven" if start else "entirely within excluded startup window",
               "mean_centered_eddy_to_mean_ke_transfer_W_m2": integral_transfer,
               "mean_mass_streamfunction_max_abs_kg_m_s": float(np.max(np.abs(means["mass_streamfunction"]))),
               "mass_streamfunction_temporal_std_spatial_rms_kg_m_s": float(np.sqrt(np.mean(psi_std**2))),
               "mass_streamfunction_abs_mean_over_mean_abs": coherence,
               "coherence_definition": "sum(abs(time-mean Psi))/sum(time-mean abs(Psi)); descriptive only",
               "native_sampled_steps": len(checks),
               "native_step_max_absolute_residual_kg_m2_s": float(checks["absolute_residual"].max()),
               "native_step_max_relative_to_leading_increment": float(checks["relative_to_leading_increment"].max()),
               "native_step_max_roundoff_units": float(checks["roundoff_units"].max()),
               "max_mass_divergence_kg_m3_s": max(r["mass_divergence_max_kg_m3_s"] for r in rows),
               "max_mass_flux_reconstruction_error_kg_m2_s": max(r["mass_flux_reconstruction_max_kg_m2_s"] for r in rows),
               "metadata_sha256": digest(folder / "metadata.toml"),
               "native_checks_sha256": digest(folder / "native_budget/checks.csv"),
               "scientific_driver_sha256": m["script_sha256"],
               "configuration_sha256": m["configuration_sha256"],
               "heights_actual_m": m["video_heights_actual_m"],
               "wall_seconds": m["wall_seconds"], "frame_manifest": frames,
               "first": rows[0], "last": rows[-1]}
    summary["direct_tracer_flux_rms_of_time_mean_kg_m2_s"] = {
        key:float(np.sqrt(np.mean(means[key]**2))) for key in
        ("mean_rho_vtracer", "mean_rho_wtracer", "zonal_rho_vtracer", "zonal_rho_wtracer")}
    if name.startswith("axis"):
        summary["axisymmetric_full_window_checks"] = {
            "max_velocity_departure_m_s":max(r["axisymmetry_max_velocity_departure_m_s"] for r in rows),
            "max_abs_centered_eddy_transfer_W_m2":max(abs(r["eddy_to_mean_ke_transfer_W_m2"]) for r in rows)}
    with (out / f"{name}_diagnostics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    print("ANALYZED", name, len(files), "frames", flush=True)
    return {"name": name, "m": m, "files": files, "times": times, "limits": limits,
            "tracer_range": tracer_range, "summary": summary, "means": means, "rows": rows}


def movie(cases, destination, limits, crange, fps):
    times = cases[0]["times"]
    # Preserve all frames while allowing a seven-frame six-hour control to be
    # watched for seven seconds rather than flashing by at production cadence.
    fps = min(fps, max(1, len(times)//8))
    reference = cases[0]["m"]
    for case in cases[1:]:
        if times.shape != case["times"].shape or not np.allclose(times, case["times"], rtol=0, atol=1e-6):
            raise ValueError("Synchronized movie requires identical times")
        for key in ("lx", "ly", "height", "f", "video_heights_actual_m"):
            if case["m"][key] != reference[key]:
                raise ValueError(f"Matched movie mismatch: {key}")
    fig, axes = plt.subplots(2*len(cases), 3, figsize=(15, 7.2*len(cases)), layout="constrained", squeeze=False)
    images = []
    for index, case in enumerate(cases):
        m = case["m"]; panel = []
        with h5py.File(case["files"][0]) as data:
            for row, key in enumerate(("vorticity_slices", "tracer_slices")):
                values = data[key][:] / m["f"] if row == 0 else data[key][:]
                for col, height in enumerate(m["video_heights_actual_m"]):
                    ax = axes[2*index+row, col]
                    norm = AsinhNorm(linear_width=.1 if row == 0 else .01,
                                     vmin=-limits[col] if row == 0 else crange[0],
                                     vmax=limits[col] if row == 0 else crange[1], clip=False)
                    im = ax.imshow(values[col], origin="lower", interpolation="nearest", norm=norm,
                                   cmap="RdBu_r" if row == 0 else "viridis",
                                   extent=(0, m["lx"]/1000, -m["ly"]/2000, m["ly"]/2000))
                    quantity = "ζ/f₀ · asinh" if row == 0 else "tracer C · asinh"
                    ax.set(title=f'{LABELS[case["name"]]} · {quantity}\nz = {height/1000:.2f} km above bottom', xlabel="x [km]", ylabel="y [km]")
                    # Four x ticks and five y ticks stay legible in each depth panel.
                    ax.set_xticks(np.linspace(0, m["lx"]/1000, 4))
                    ax.set_yticks(np.linspace(-m["ly"]/2000, m["ly"]/2000, 5))
                    ax.title.set_fontsize(9)
                    cb = fig.colorbar(im, ax=ax, shrink=.85)
                    ticks = ([-5,-3,-1,-.3,-.1,0,.1,.3,1,3,5] if row == 0 else [0,.01,.03,.1,.3,1])
                    ticks = [v for v in ticks if norm.vmin <= v <= norm.vmax]
                    cb.set_ticks(ticks); cb.set_ticklabels([f"{v:g}" for v in ticks])
                    panel.append((key, col, im))
        images.append(panel)
    title = fig.suptitle("Ferrel mechanism screening · Earth day 0.000", fontsize=12)
    fig.canvas.draw(); fig.set_layout_engine(None)
    writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-crf", "29", "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    with writer.saving(fig, str(destination), dpi=90):
        for i, t in enumerate(times):
            for case, panel in zip(cases, images):
                with h5py.File(case["files"][i]) as data:
                    for key, level, im in panel:
                        field = data[key][level]
                        im.set_data(field / case["m"]["f"] if key == "vorticity_slices" else field)
            title.set_text(f'Ferrel mechanism screening · Earth day {t/86400:.3f} · {"startup window" if t <= 43200 else "transient evolution"}')
            writer.grab_frame()
    fig.savefig(destination.with_suffix(".png"), dpi=90); plt.close(fig)
    with Image.open(destination.with_suffix(".png")) as poster:
        poster_size=list(poster.size);poster.verify()
    count, duration = imageio_ffmpeg.count_frames_and_secs(str(destination))
    if count != len(times) or abs(duration - count/fps) > .05:
        raise ValueError("Encoded movie frame count/duration mismatch")
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-v", "error", "-i", str(destination), "-f", "null", "-"], check=True)
    probe_path=shutil.which("ffprobe")
    if probe_path:
        probe=json.loads(subprocess.check_output([probe_path,"-v","error","-select_streams","v:0","-show_entries","stream=codec_name,pix_fmt,width,height,avg_frame_rate","-of","json",str(destination)],text=True))["streams"][0]
        if probe["codec_name"]!="h264" or probe["pix_fmt"]!="yuv420p":
            raise ValueError("Unexpected browser video encoding")
        if [probe["width"],probe["height"]]!=poster_size:
            raise ValueError("Movie/poster dimensions differ")
    else:
        reader=imageio_ffmpeg.read_frames(str(destination))
        try: probe=next(reader)
        finally: reader.close()
        if probe.get("codec")!="h264" or "yuv420p" not in probe.get("pix_fmt",""):
            raise ValueError("Unexpected browser video encoding")
        if list(probe["size"])!=poster_size:
            raise ValueError("Movie/poster dimensions differ")
    if destination.stat().st_size > 24_000_000:
        raise ValueError("Finished movie exceeds 24 MB publication limit; re-render at explicit lower quality on Unity")
    provenance = {"file": destination.name, "poster": destination.with_suffix(".png").name,
                  "cases": [c["name"] for c in cases], "layout": f"{2*len(cases)} rows × 3 actual depth columns",
                  "row_order": [f'{c["name"]}: {field}' for c in cases for field in ("vorticity", "tracer")],
                  "frames": count, "fps": fps, "duration_seconds": duration,
                  "times_seconds": times.tolist(), "heights_actual_m": reference["video_heights_actual_m"],
                  "vorticity_asinh_width": .1, "vorticity_limits_by_depth": [[-float(v),float(v)] for v in limits],
                  "tracer_asinh_width": .01, "tracer_limits": crange,
                  "fixed_across_time": True, "shared_between_cases": True,
                  "no_field_clipping_or_interpolation": True, "sha256": digest(destination),
                  "bytes": destination.stat().st_size, "decode_validation": "PASS"}
    provenance["stream_probe"]=probe
    destination.with_suffix(".json").write_text(json.dumps(provenance, indent=2)+"\n")
    print("RENDERED_AND_DECODED", destination.name, count, "frames", destination.stat().st_size, "bytes", flush=True)
    return provenance


def figures(cases, out):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    metrics = [("mean_ke_J_m2", "Zonal-mean zonal KE [J/m²]"),
               ("max_abs_mass_streamfunction_kg_m_s", "Max |native Ψ| [kg/(m s)]"),
               ("eddy_to_mean_ke_transfer_W_m2", "Centered eddy→mean KE transfer [W/m²]"),
               ("centered_advection_remainder_rms_kg_m2_s2", "Centered/native advection remainder RMS [kg/(m² s²)]")]
    for ax, (metric, label) in zip(axes.flat, metrics):
        for case in cases:
            ax.plot(case["times"]/86400, [r[metric] for r in case["rows"]], label=LABELS[case["name"]])
        ax.axvspan(0, .5, alpha=.08, color="black"); ax.axhline(0, color="black", linewidth=.5)
        ax.set(xlabel="Earth day", ylabel=label); ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    fig.suptitle("Transient Ferrel screening · shaded first 12 h excluded as startup")
    fig.savefig(out / "diagnostics.png", dpi=130); plt.close(fig)
    for case in cases:
        m, means = case["m"], case["means"]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout="constrained")
        fields = [(means["mass_streamfunction"], "Native mass streamfunction Ψ [kg/(m s)]"),
                  (means["zonal_u"], "Zonal mean u [m/s]"),
                  (means["eddy_to_mean_ke_transfer_y"]+means["eddy_to_mean_ke_transfer_z"], "Centered eddy→mean KE transfer [W/m³]")]
        for ax, (data, label) in zip(axes, fields):
            limit = max(float(np.max(np.abs(data))), 1e-30)
            yf=np.linspace(-m["ly"]/2000,m["ly"]/2000,m["ny"]+1)
            zf=np.linspace(0,m["height"]/1000,m["nz"]+1)
            if data.shape==(m["nz"]+1,m["ny"]+1):
                im=ax.pcolormesh(yf,zf,data,shading="gouraud",cmap="RdBu_r",vmin=-limit,vmax=limit)
            else:
                im=ax.pcolormesh(yf,zf,data,shading="flat",cmap="RdBu_r",vmin=-limit,vmax=limit)
            ax.set(title=label, xlabel="y [km]", ylabel="height above bottom [km]")
            ax.title.set_fontsize(9); fig.colorbar(im, ax=ax, shrink=.8)
        window=np.asarray(case["summary"]["window_seconds"])/3600
        fig.suptitle(f'{LABELS[case["name"]]} · time mean {window[0]:g}–{window[1]:g} h · {case["summary"]["window_scope"]}')
        fig.savefig(out/f'{case["name"]}_circulation.png',dpi=130);plt.close(fig)
        fig,axes=plt.subplots(2,3,figsize=(14,8),layout="constrained")
        fieldlist=[(means["native_rhou_rhs"][i],f"Native {term} RHS") for i,term in enumerate(TERMS)]
        fieldlist.append((means["centered_native_advection_remainder"],"Centered/native advection remainder"))
        for ax,(data,label) in zip(axes.flat,fieldlist):
            limit=max(float(np.max(np.abs(data))),1e-30)
            im=ax.pcolormesh(yf,zf,data,shading="flat",cmap="RdBu_r",vmin=-limit,vmax=limit)
            ax.set(title=label,xlabel="y [km]",ylabel="height [km]");ax.title.set_fontsize(9)
            fig.colorbar(im,ax=ax,shrink=.8,label="kg/(m² s²)")
        axes.flat[-1].axis("off")
        fig.suptitle(f'{LABELS[case["name"]]} · mean sampled snapshot RHS {window[0]:g}–{window[1]:g} h\nNative step-closure checks are separate; this figure is not a closed long-run integral')
        fig.savefig(out/f'{case["name"]}_momentum.png',dpi=110);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
        for row,axis in enumerate(("v","w")):
            fields=[means[f"mean_rho_{axis}tracer"],means[f"zonal_rho_{axis}tracer"]]
            limit=max(float(np.max(np.abs(fields))),1e-30)
            for col,(data,label) in enumerate(zip(fields,("Eulerian mean","Eddy covariance"))):
                ax=axes[row,col]
                im=ax.pcolormesh(yf,zf,data,shading="flat",cmap="RdBu_r",vmin=-limit,vmax=limit)
                ax.set(title=f"{label} rho*{axis}*C",xlabel="y [km]",ylabel="height [km]")
                fig.colorbar(im,ax=ax,shrink=.8,label="kg/(m² s)")
        fig.suptitle(f'{LABELS[case["name"]]} · centered tracer fluxes {window[0]:g}–{window[1]:g} h\nDirect mean/eddy decomposition; not native WENO face transport')
        fig.savefig(out/f'{case["name"]}_tracer_transport.png',dpi=110);plt.close(fig)


def comparisons(cases):
    result={}; by_name={case["name"]:case for case in cases}
    for left,right in (("b0","b162"),("b162","boundary162"),("axis0","axis162"),("axis0","b0"),("axis162","b162"),("b162","fine162")):
        if left not in by_name or right not in by_name: continue
        a,b=by_name[left],by_name[right]
        matches=[(i,j) for i,t in enumerate(a["times"]) for j,s in enumerate(b["times"]) if abs(t-s)<=1e-6]
        if len(matches)<2: continue
        ia,ib=map(np.asarray,zip(*matches));shared=a["times"][ia]
        metrics={}
        for key in ("mean_ke_J_m2","eddy_to_mean_ke_transfer_W_m2","max_abs_mass_streamfunction_kg_m_s","density_weighted_zonal_mean_v_rms_m_s","density_weighted_zonal_mean_w_rms_m_s"):
            delta=np.asarray([b["rows"][i][key] for i in ib])-np.asarray([a["rows"][i][key] for i in ia])
            start=43200.0 if shared[-1]>43200 else 0.0
            mean,window=time_mean(delta,shared,start)
            metrics[key]={"mean_difference":float(mean),"final_difference":float(delta[-1])}
        result[f"{right}_minus_{left}"]={"window_seconds":window,"common_times":len(shared),"metrics":metrics,
            "scope":"Matched diagnostics at coincident saved times; not statistical significance or fieldwise convergence."}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",choices=("pilot","production"),required=True)
    p.add_argument("--production-root",type=Path,default=Path("data/ferrel2026-gpu-v1"))
    p.add_argument("--preflight-root",type=Path,default=Path("data/ferrel2026-preflight-gpu-v1"))
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--include-fine",action="store_true")
    p.add_argument("--fps",type=int,default=6)
    p.add_argument("--completed-job",action="append",required=True,help="Successful prerequisite Slurm job; repeat for independent refinement")
    args=p.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        p.error("Run analysis/rendering inside an allocated Unity job")
    if args.fps<1: p.error("fps must be positive")
    if not all(job.isdigit() for job in args.completed_job): p.error("Completed job IDs must be integers")
    accounting=subprocess.check_output(["sacct","-X","-j",",".join(args.completed_job),"-n","-P","-o","JobIDRaw,State,ExitCode"],text=True)
    states={row[0]:row[1:] for line in accounting.splitlines() if (row:=line.strip().split("|")) and len(row)>=3}
    for job in args.completed_job:
        if job not in states or states[job][:2]!=["COMPLETED","0:0"]:
            raise ValueError(f"Prerequisite job {job} has not completed successfully: {states.get(job)}")
    if args.output.exists() and any(args.output.iterdir()):
        p.error("Use a fresh output directory; refusing to overwrite publication assets")
    args.output.mkdir(parents=True,exist_ok=True)
    gates=[gate(args.preflight_root/"PASS512.toml")]
    if args.stage=="pilot":
        sources={name:args.preflight_root/f"pilot512_{name}" for name in ("b0","b162","axis0","axis162")}
    else:
        sources={name:args.production_root/f"ferrel2026_{name}" for name in ("b0","b162","boundary162")}
        sources.update({name:args.preflight_root/f"pilot512_{name}" for name in ("axis0","axis162")})
    if args.include_fine:
        gates.append(gate(args.preflight_root/"PASS1024.toml"))
        sources["fine162"]=args.preflight_root/"pilot1024_b162"
    cases=[collect(name,folder,args.output) for name,folder in sources.items()]
    limits=np.maximum(np.max([c["limits"] for c in cases],axis=0),1e-12)
    crange=[min(c["tracer_range"][0] for c in cases),max(c["tracer_range"][1] for c in cases)]
    media=[]
    for case in cases:
        media.append(movie([case],args.output/f'{case["name"]}_vorticity_tracer.mp4',limits,crange,args.fps))
    by_name={case["name"]:case for case in cases}
    for left,right,title in (("b0","b162","heating_comparison"),("b162","boundary162","boundary_comparison"),("axis0","axis162","axisymmetric_comparison")):
        if left in by_name and right in by_name:
            media.append(movie([by_name[left],by_name[right]],args.output/f"{title}.mp4",limits,crange,args.fps))
    figures(cases,args.output)
    report={"status":"COMPLETE","stage":args.stage,"compute_location":"Unity allocated CPU job",
            "generated_at_utc":datetime.now(timezone.utc).isoformat(),
            "slurm_job_id":os.environ["SLURM_JOB_ID"],"renderer_sha256":digest(__file__),
            "successful_prerequisite_jobs":{job:states[job] for job in args.completed_job},
            "raw_fields_downloaded":False,"gates":gates,"caveats":CAVEATS,
            "cases":{c["name"]:c["summary"] for c in cases},"comparisons":comparisons(cases),"media":media,
            "normalization":{"vorticity":"asinh ζ/f0, width0.1; per-depth full extrema shared over all cases and times",
                             "tracer":"asinh C, width0.01; physical bounds expanded to observed extrema; common across all cases, depths and times"}}
    (args.output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    text=["# Ferrel transient mechanism screening", "",f"Completed {args.stage} publication assets; generated and decoded on Unity.","",
          "[Watch on the GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/).", "",
          "| Case | Grid | Duration (Earth h) | Native sampled closure steps |", "|---|---|---:|---:|"]
    for c in cases:
        s=c["summary"];text.append(f'| {s["label"]} | {" × ".join(map(str,s["grid"]))} | {s["stop_seconds"]/3600:g} | {s["native_sampled_steps"]} |')
    text += ["","These outputs do not yet establish statistically steady Ferrel cells. The first 12 hours are conservatively excluded as startup, without claiming a measured equilibration time; six-hour controls/refinement probe this early window only.","",
             "Movies use actual saved vorticity/tracer at three heights. Time-independent asinh color scales are shared across cases, with a common tracer range and separate full-range vorticity limits for each depth. No automatic per-frame rescaling, field clipping or inferred depths.","",
             "See `report.json` for numeric diagnostics, exact averaging windows, limitations, input hashes and media checksums. Native step closure and approximate centered covariance interpretation are reported separately.",""]
    for caveat in CAVEATS:text.append("- "+caveat)
    (args.output/"README.md").write_text("\n".join(text)+"\n")
    assets=[{"file":f.name,"bytes":f.stat().st_size,"sha256":digest(f)} for f in sorted(args.output.iterdir()) if f.is_file()]
    (args.output/"publication-manifest.json").write_text(json.dumps({"status":"READY_FOR_REVIEW","stage":args.stage,"assets":assets},indent=2)+"\n")
    print("READY_FOR_REVIEW",args.output,flush=True)


if __name__=="__main__":main()
