"""Analyze and publish the completed finite-mode campaign; execute on Unity only.

Usage: python analyze.py CAMPAIGN_ROOT RENDER_OUTPUT

Scientific integration is not performed here. Only completed saved diagnostics
and 64 x 64 modal reconstructions are read. The output is a self-contained gallery
section; this script never edits the repository landing page or uploads files.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np


BUDGETS = ("kinetic", "thermal_variance", "potential_moment", "barotropic_kinetic")
DISPLAY_NAMES = {
    "plus-seed18": "Convective seed A, f = +1",
    "plus-seed19": "Convective seed B, f = +1",
    "minus-reflected18": "Reflected seed A, f = -1",
    "zero-gradient18": "Zero imposed gradient, seed A",
    "refined-plus18": "Nested extension, seed A, f = +1",
}


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def series(case, key):
    values = np.asarray([row[key] for row in case["rows"]], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite diagnostic {case['name']}: {key}")
    return values


def cumulative_trapezoid(values, times):
    increments = np.diff(times) * (values[1:] + values[:-1]) / 2
    return np.concatenate(([0.0], np.cumsum(increments)))


def weighted_mean(values, times):
    if len(times) < 2 or times[-1] == times[0]:
        return float(values[-1])
    return float(cumulative_trapezoid(values, times)[-1] / (times[-1] - times[0]))


def load_completed(root):
    root = Path(root)
    config = read_json(root / "config.json")
    if not (root / "validation" / "PASS.json").is_file():
        raise RuntimeError("No passing release-validation record; publication blocked")
    cases = []
    for spec in config["cases"]:
        directory = root / "runs" / spec["name"]
        required = [directory / name for name in
                    ("COMPLETE.json", "diagnostics.json", "planes.npz", "metadata.json")]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise RuntimeError("Incomplete case; publication blocked: " + ", ".join(missing))
        complete = read_json(directory / "COMPLETE.json")
        metadata = read_json(directory / "metadata.json")
        rows = read_json(directory / "diagnostics.json")
        if not rows:
            raise ValueError(f"No diagnostics for {spec['name']}")
        with np.load(directory / "planes.npz", allow_pickle=False) as archive:
            planes = {key: np.asarray(archive[key]) for key in
                      ("times", "heights", "x", "y", "vorticity", "buoyancy_perturbation")}
        times = np.asarray([row["time"] for row in rows], dtype=float)
        frames = planes["times"]
        if not (np.isfinite(times).all() and np.isfinite(frames).all()
                and np.all(np.diff(times) > 0) and np.all(np.diff(frames) > 0)):
            raise ValueError(f"Invalid time coordinate for {spec['name']}")
        if (abs(times[0]) > 1e-10 or abs(frames[0]) > 1e-10
                or not np.isclose(times[-1], config["duration"])
                or not np.isclose(frames[-1], config["duration"])):
            raise ValueError(f"Saved output does not cover full declared interval: {spec['name']}")
        if len(planes["heights"]) != 3:
            raise ValueError("Three actual saved heights are required")
        expected = (len(frames), 3, config["nxy"], config["nxy"])
        for key in ("vorticity", "buoyancy_perturbation"):
            if planes[key].shape != expected or not np.isfinite(planes[key]).all():
                raise ValueError(f"Invalid saved field {spec['name']}: {key}")
        if int(complete["saved_frames"]) != len(frames):
            raise ValueError(f"Frame count disagrees with completion record: {spec['name']}")
        if metadata["case"] != spec:
            raise ValueError(f"Case metadata mismatch: {spec['name']}")
        cases.append({"name": spec["name"], "spec": spec, "directory": directory,
                      "metadata": metadata, "complete": complete,
                      "rows": rows, "time": times, "planes": planes})
    return config, cases


def color_scales(cases):
    scales = {}
    for key, cmap, label in (
        ("vorticity", "RdBu_r", "Vertical vorticity zeta / |f0|"),
        ("buoyancy_perturbation", "BrBG", "Buoyancy perturbation theta / (|f0|^2 H)"),
    ):
        peak = max(float(np.max(np.abs(case["planes"][key]))) for case in cases)
        # A zero field needs a nonzero plotting interval but stays exactly neutral.
        limit = peak if peak > 0 else 1.0
        scales[key] = {"label": label, "colormap": cmap,
                       "physical_min": -limit, "physical_max": limit,
                       "actual_global_absolute_max": peak,
                       "linear_width": limit / 1000,
                       "mapping": "asinh(value / linear_width) / asinh(physical_max / linear_width)",
                       "fixed_over": "all cases, all times, all three heights",
                       "clipping": False,
                       "zero_field_fallback": peak == 0}
    return scales


def calculate_metrics(config, cases):
    result = {"campaign": config["campaign"], "config": config,
              "meaning": config["meaning"],
              "control_caveat": config["control_caveat"],
              "statistical_window": "final quarter of each completed interval",
              "third_moment_integral_method": "trapezoidal quadrature of saved diagnostic samples; not an RK-stage closure test",
              "cases": {}}
    for case in cases:
        t = case["time"]
        mask = t >= t[0] + .75 * (t[-1] - t[0])
        sign = case["spec"]["f_sign"]
        aligned = sign * series(case, "barotropic_vorticity_skewness")
        final = case["rows"][-1]
        defects = {key: float(np.max(np.abs(series(case, key + "_integrated_budget_relative_defect"))))
                   for key in BUDGETS}
        if max(defects.values()) > 5e-5:
            raise RuntimeError(f"Budget release limit exceeded in {case['name']}")
        rates = {}
        for short, key in (
            ("convective", "barotropic_third_moment_convective"),
            ("projected_barotropic_self_advection", "barotropic_third_moment_self_advective"),
            ("viscous", "barotropic_third_moment_viscous"),
        ):
            values = series(case, key)
            rates[short] = {"signed_integral": float(cumulative_trapezoid(values, t)[-1]),
                            "f_aligned_integral": float(sign * cumulative_trapezoid(values, t)[-1]),
                            "absolute_integral": float(cumulative_trapezoid(np.abs(values), t)[-1]),
                            "final_quarter_mean": weighted_mean(values[mask], t[mask])}
        third_rhs = series(case, "barotropic_vorticity_third_moment_tendency")
        summed = sum(series(case, key) for key in (
            "barotropic_third_moment_convective", "barotropic_third_moment_self_advective",
            "barotropic_third_moment_viscous", "barotropic_third_moment_coriolis",
            "barotropic_third_moment_buoyancy"))
        third_abs_residual = float(np.max(np.abs(third_rhs - summed)))
        rms = series(case, "barotropic_vorticity_rms")
        usable = mask & (rms > max(1e-12, 1e-8 * float(np.max(rms))))
        signed_mean = weighted_mean(aligned[usable], t[usable]) if np.any(usable) else None
        metrics = {
            "title": DISPLAY_NAMES.get(case["name"], case["name"]),
            "duration": float(t[-1]), "frames": int(len(case["planes"]["times"])),
            "saved_heights": [float(z) for z in case["planes"]["heights"]],
            "K": case["spec"]["K"], "Nv": case["spec"]["Nv"],
            "f": case["metadata"]["model"]["f"],
            "G": case["metadata"]["model"]["G"],
            "final_kinetic": float(final["kinetic"]),
            "final_barotropic_kinetic": float(final["barotropic_kinetic"]),
            "final_barotropic_energy_fraction": float(final["barotropic_energy_fraction"]),
            "final_lowest_barotropic_shell_kinetic": float(final["lowest_barotropic_shell_kinetic"]),
            "final_lowest_horizontal_shell_kinetic": float(final["lowest_horizontal_shell_kinetic"]),
            "final_vorticity_skewness": float(final["barotropic_vorticity_skewness"]),
            "final_quarter_f_aligned_skewness_mean": signed_mean,
            "final_quarter_f_aligned_skewness_min": float(np.min(aligned[usable])) if np.any(usable) else None,
            "final_quarter_f_aligned_skewness_max": float(np.max(aligned[usable])) if np.any(usable) else None,
            "final_quarter_barotropic_power_mean": weighted_mean(series(case, "barotropic_power")[mask], t[mask]),
            "final_quarter_barotropic_energy_mean": weighted_mean(series(case, "barotropic_kinetic")[mask], t[mask]),
            "final_quarter_outer_band_energy_fraction_mean": weighted_mean(series(case, "outer_horizontal_band_energy_fraction")[mask], t[mask]),
            "max_integrated_budget_relative_defect": defects,
            "max_relative_divergence": float(np.max(series(case, "relative_divergence"))),
            "max_advective_cfl": float(np.max(series(case, "advective_cfl_sum"))),
            "max_potential_exchange_absolute_residual": float(np.max(np.abs(series(case, "potential_buoyancy_exchange_residual")))),
            "max_stress_exchange_absolute_residual": float(np.max(np.abs(series(case, "barotropic_power") - series(case, "barotropic_stress_power")))),
            "max_third_moment_rhs_absolute_residual": third_abs_residual,
            "third_moment_contributions": rates,
            "integration_elapsed_seconds": float(case["complete"]["elapsed_seconds"]),
            "source_sha256": case["metadata"]["source_sha256"],
            "input_sha256": {name: checksum(case["directory"] / name) for name in
                             ("diagnostics.json", "metadata.json", "COMPLETE.json", "planes.npz")},
        }
        result["cases"][case["name"]] = metrics
    by_name = {case["name"]: case for case in cases}
    if {"plus-seed18", "minus-reflected18"} <= set(by_name):
        plus, minus = (by_name[key] for key in ("plus-seed18", "minus-reflected18"))
        if not np.array_equal(plus["planes"]["times"], minus["planes"]["times"]):
            raise ValueError("Reflected pair does not have matching frame times")
        reflection = {}
        for key, parity in (("vorticity", -1), ("buoyancy_perturbation", 1)):
            a, b = plus["planes"][key], minus["planes"][key]
            reflected = parity * np.take(a, (-np.arange(a.shape[-2])) % a.shape[-2], axis=-2)
            residual = b - reflected
            reflection[key] = {
                "relative_l2_over_saved_sequence": float(np.linalg.norm(residual.ravel()) / max(np.linalg.norm(a.ravel()), 1e-300)),
                "max_absolute_residual": float(np.max(np.abs(residual))),
            }
        result["reflected_pair_saved_fields"] = reflection
    return result


def make_scientific_plot(cases, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
    colors = plt.get_cmap("tab10")(np.arange(len(cases)))
    for case, color in zip(cases, colors):
        t = case["time"]
        name = case["name"]
        sign = case["spec"]["f_sign"]
        axes[0, 0].semilogy(t, np.maximum(series(case, "barotropic_kinetic"), 1e-30), color=color, label=name)
        aligned = sign * series(case, "barotropic_vorticity_skewness")
        rms = series(case, "barotropic_vorticity_rms")
        aligned[rms <= max(1e-12, 1e-8 * float(np.max(rms)))] = np.nan
        axes[0, 1].plot(t, aligned, color=color)
        axes[1, 0].semilogy(t, np.maximum(series(case, "lowest_barotropic_shell_kinetic"), 1e-30), color=color)
        for key, style in (("barotropic_third_moment_convective", "-"),
                           ("barotropic_third_moment_self_advective", "--")):
            cumulative = sign * cumulative_trapezoid(series(case, key), t)
            axes[1, 1].plot(t, cumulative, style, color=color)
        axes[2, 0].plot(t, series(case, "barotropic_power"), color=color)
        defects = np.maximum.reduce([np.abs(series(case, key + "_integrated_budget_relative_defect")) for key in BUDGETS])
        axes[2, 1].semilogy(t, np.maximum(defects, 1e-18), color=color)
    axes[0, 0].set_title("Absolute barotropic kinetic energy")
    axes[0, 0].set_ylabel("E_BT")
    axes[0, 0].legend(fontsize=7, loc="best")
    axes[0, 1].set_title("Rotation-aligned barotropic skewness")
    axes[0, 1].set_ylabel("sign(f) <Z^3> / <Z^2>^(3/2)")
    axes[1, 0].set_title("Absolute lowest-shell barotropic energy")
    axes[1, 0].set_ylabel("E_BT in 0 < mx^2 + my^2 <= 2")
    axes[1, 1].set_title("Cubic-moment production: convective vs projection")
    axes[1, 1].set_ylabel("sign(f) x cumulative production")
    axes[1, 1].plot([], [], "k-", label="Convective stress")
    axes[1, 1].plot([], [], "k--", label="Projected BT self-advection")
    axes[1, 1].legend(fontsize=8)
    axes[2, 0].set_title("Signed convective transfer into barotropic flow")
    axes[2, 0].set_ylabel("P_BT (positive supplies barotropic energy)")
    axes[2, 1].set_title("Largest integrated budget defect")
    axes[2, 1].set_ylabel("max relative defect: K, Q, M, E_BT")
    axes[2, 1].axhline(5e-5, color="black", linestyle=":", linewidth=1)
    for ax in axes.flat:
        ax.set_xlabel("Time |f0| t")
        ax.grid(alpha=.2)
    for ax in (axes[0, 1], axes[1, 1], axes[2, 0]):
        ax.axhline(0, color="0.6", linewidth=.7)
    axes[1, 1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axes[2, 0].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    fig.suptitle("Finite-mode rotating-convection reference: all five completed cases", fontsize=14)
    fig.savefig(Path(output) / "diagnostics.png", dpi=140)
    fig.savefig(Path(output) / "diagnostics.pdf")
    plt.close(fig)


def mapped_rgb(values, scale, cmap):
    limit, width = scale["physical_max"], scale["linear_width"]
    normalized = .5 + .5 * np.arcsinh(values / width) / np.arcsinh(limit / width)
    if float(np.min(normalized)) < -1e-12 or float(np.max(normalized)) > 1 + 1e-12:
        raise ValueError("Field exceeds fixed publication color scale")
    return np.asarray(cmap(np.clip(normalized, 0, 1), bytes=True)[..., :3], dtype=np.uint8)


def verify_movie(path, expected_frames, fps, dimensions):
    command = ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
               "-show_entries", "stream=width,height,nb_read_frames,codec_name,duration",
               "-of", "json", str(path)]
    probe = subprocess.run(command, check=True, text=True, capture_output=True)
    data = json.loads(probe.stdout)
    stream = data["streams"][0]
    if (int(stream["nb_read_frames"]) != expected_frames
            or (int(stream["width"]), int(stream["height"])) != dimensions):
        raise RuntimeError(f"Encoded movie dimensions or frame count failed: {path}")
    expected_duration = expected_frames / fps
    if abs(float(stream["duration"]) - expected_duration) > 1 / fps + 1e-6:
        raise RuntimeError(f"Encoded movie duration failed: {path}")
    decoded = subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-threads", "1",
                              "-i", str(path), "-f", "null", "-"],
                             check=True, text=True, capture_output=True)
    if decoded.stderr.strip():
        raise RuntimeError(f"Movie decoder reported errors: {decoded.stderr}")
    return {"ffprobe": stream, "full_decode_passed": True,
            "sha256": checksum(path), "bytes": path.stat().st_size}


def render_movie(case, scales, output, fps):
    from PIL import Image, ImageDraw, ImageFont
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import colormaps, font_manager

    font_path = font_manager.findfont("DejaVu Sans")
    fonts = {size: ImageFont.truetype(font_path, size) for size in (13, 15, 17, 20, 23)}
    dimensions = (1100, 840)
    panel = 256
    columns = [68, 365, 662]
    tops = [155, 505]
    base = Image.new("RGB", dimensions, "white")
    draw = ImageDraw.Draw(base)
    title = DISPLAY_NAMES.get(case["name"], case["name"])
    model = case["metadata"]["model"]
    draw.text((30, 12), title, fill="#172336", font=fonts[23])
    draw.text((30, 46), "Nonlinear Galerkin Boussinesq reference; fixed-temperature, free-slip walls; no cap", fill="#334155", font=fonts[17])
    params = (f"K_h={model['K']}, N_v={model['Nv']}    L/H={model['horizontal_period']:.5f}    "
              f"G={model['G']:.6g}    nu=kappa={model['nu']:.1e}    f={model['f']:+.0f}")
    draw.text((30, 72), params, fill="#334155", font=fonts[15])
    fields = ("vorticity", "buoyancy_perturbation")
    cmaps = {key: colormaps[scales[key]["colormap"]] for key in fields}
    for row, key in enumerate(fields):
        y = tops[row]
        label = "Vertical vorticity zeta / |f0|" if row == 0 else "Buoyancy perturbation theta / (|f0|^2 H), NOT a tracer"
        draw.text((columns[0], y - 49), label, fill="#172336", font=fonts[17])
        for col, height in enumerate(case["planes"]["heights"]):
            x = columns[col]
            draw.text((x + 164, y - 22), f"z/H = {height:g}", fill="#172336", font=fonts[15])
            draw.rectangle((x - 1, y - 1, x + panel, y + panel), outline="#334155")
            draw.text((x, y + panel + 3), "0", fill="#475569", font=fonts[13])
            draw.text((x + 111, y + panel + 3), "x/H", fill="#475569", font=fonts[13])
            draw.text((x + panel - 37, y + panel + 3), f"{model['horizontal_period']:.3f}", fill="#475569", font=fonts[13])
        draw.text((16, y + 115), "y/H", fill="#475569", font=fonts[13])
        draw.text((13, y), f"{model['horizontal_period']:.3f}", fill="#475569", font=fonts[13])
        draw.text((42, y + panel - 15), "0", fill="#475569", font=fonts[13])
        bar_x, bar_w = 954, 23
        transformed = np.linspace(1, -1, panel)
        scale = scales[key]
        values = scale["linear_width"] * np.sinh(transformed * np.arcsinh(scale["physical_max"] / scale["linear_width"]))
        bar = mapped_rgb(values[:, None], scale, cmaps[key])
        bar = np.repeat(bar, bar_w, axis=1)
        base.paste(Image.fromarray(bar), (bar_x, y))
        draw.rectangle((bar_x - 1, y - 1, bar_x + bar_w, y + panel), outline="#334155")
        lim = scale["physical_max"]
        for value in (-lim, -.01 * lim, 0, .01 * lim, lim):
            unit = np.arcsinh(value / scale["linear_width"]) / np.arcsinh(lim / scale["linear_width"])
            yy = y + round((1 - unit) * (panel - 1) / 2)
            draw.line((bar_x + bar_w, yy, bar_x + bar_w + 5, yy), fill="#334155")
            draw.text((bar_x + bar_w + 8, yy - 7), f"{value:.2g}", fill="#334155", font=fonts[13])
        draw.text((bar_x - 3, y + panel + 4), "fixed asinh", fill="#334155", font=fonts[13])
    footer = (f"Exact retained-mode reconstruction sampled at {case['planes']['vorticity'].shape[-1]} x "
              f"{case['planes']['vorticity'].shape[-1]}; display enlarged only. Common scales across every case/time/depth.")
    draw.text((30, 798), footer, fill="#475569", font=fonts[13])
    draw.text((30, 817), "Positive zeta is cyclonic for f > 0; negative zeta is cyclonic for f < 0. Dimensionless time is not Earth days.", fill="#475569", font=fonts[13])
    destination = Path(output) / f"{case['name']}.mp4"
    log = Path(output) / f"{case['name']}-encode.log"
    command = ["ffmpeg", "-y", "-v", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{dimensions[0]}x{dimensions[1]}", "-r", str(fps), "-i", "-",
               "-an", "-c:v", "libx264", "-preset", "ultrafast", "-crf", "24",
               "-threads", "1", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(destination)]
    with log.open("wb") as logstream:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=logstream)
        try:
            for frame_index, timestamp in enumerate(case["planes"]["times"]):
                frame = base.copy()
                for row, key in enumerate(fields):
                    for col in range(3):
                        # Saved coordinates are [x,y]; image rows run from high y to low y.
                        values = case["planes"][key][frame_index, col].T[::-1, :]
                        rgb = mapped_rgb(values, scales[key], cmaps[key])
                        image = Image.fromarray(rgb).resize((panel, panel), Image.Resampling.NEAREST)
                        frame.paste(image, (columns[col], tops[row]))
                textdraw = ImageDraw.Draw(frame)
                textdraw.text((820, 14), f"|f0| t = {timestamp:6.1f}", fill="#172336", font=fonts[20])
                if frame_index == len(case["planes"]["times"]) - 1:
                    frame.save(Path(output) / f"{case['name']}-poster.png")
                process.stdin.write(frame.tobytes())
            process.stdin.close()
            if process.wait() != 0:
                raise RuntimeError(f"Movie encoding failed; inspect {log}")
        except BaseException:
            if process.poll() is None:
                process.kill()
                process.wait()
            raise
    verification = verify_movie(destination, len(case["planes"]["times"]), fps, dimensions)
    verification.update(file=destination.name, poster=f"{case['name']}-poster.png", fps=fps,
                        field_pixels=[64, 64], displayed_panel_pixels=[panel, panel],
                        spatial_display_interpolation="nearest; no additional dynamical resolution")
    return verification


def fmt(value):
    return "undefined (negligible barotropic signal)" if value is None else f"{value:.4g}"


def narrative(metrics):
    names = metrics["cases"]
    paragraphs = [
        f"All {len(names)} declared cases completed through |f0|t={metrics['config']['duration']:g}. "
        "These are finite-mode, fixed-temperature rotating Boussinesq experiments, not the capped Breeze atmosphere.",
        "The movies show actual reconstructed vertical vorticity above buoyancy perturbation at three saved heights. "
        "There is no passive tracer. Each field has one common symmetric asinh scale across all cases, times and depths.",
    ]
    if "plus-seed18" in names and "plus-seed19" in names:
        a, b = (names[key] for key in ("plus-seed18", "plus-seed19"))
        paragraphs.append("Over the final quarter, rotation-aligned barotropic skewness averaged "
                          f"{fmt(a['final_quarter_f_aligned_skewness_mean'])} for seed A and "
                          f"{fmt(b['final_quarter_f_aligned_skewness_mean'])} for seed B. "
                          "Positive aligned skewness means a cyclonic statistical bias; it does not establish a single persistent vortex.")
    if "refined-plus18" in names:
        refined = names["refined-plus18"]
        paragraphs.append("The nested extension's final-quarter aligned skewness was "
                          f"{fmt(refined['final_quarter_f_aligned_skewness_mean'])}; its mean power into the barotropic flow was "
                          f"{fmt(refined['final_quarter_barotropic_power_mean'])}. "
                          "This comparison is a truncation-sensitivity check, not a converged DNS result.")
    worst = max(max(item["max_integrated_budget_relative_defect"].values()) for item in names.values())
    paragraphs.append(f"The largest recorded integrated relative budget defect was {worst:.3g} "
                      "across kinetic energy, thermal variance, potential moment and barotropic kinetic energy. "
                      "The cubic-moment plot separates convective production from projected barotropic self-advection; "
                      "the latter need not vanish in a finite Galerkin model.")
    if "reflected_pair_saved_fields" in metrics:
        mirror = metrics["reflected_pair_saved_fields"]
        paragraphs.append("The completed reflected pair's saved-field relative L2 mismatch was "
                          f"{mirror['vorticity']['relative_l2_over_saved_sequence']:.3g} for vorticity and "
                          f"{mirror['buoyancy_perturbation']['relative_l2_over_saved_sequence']:.3g} for buoyancy. "
                          "Reflection reverses both f and vorticity, while preserving scalar buoyancy.")
    paragraphs.extend([
        "The zero-imposed-gradient control retains the same initial buoyancy perturbation and its force. "
        "It has a finite thermal reservoir and is not pure velocity-only decay.",
        "The narrow box contains only a few critical convective wavelengths. Two seeds and one extension cannot "
        "establish a robust Jovian sign-selection mechanism or a broad inverse cascade. Absolute low-shell energy "
        "and signed stress power are provided so morphology and energy fractions are not mistaken for transfer evidence.",
    ])
    return paragraphs


def write_publication(output, metrics, scales):
    output = Path(output)
    paragraphs = narrative(metrics)
    md = ["# Nonlinear Galerkin convection: completed reference experiment", ""]
    for text in paragraphs:
        md.extend([text, ""])
    md.extend(["| Case | Final E_BT | Final low-shell E_BT | Final-quarter aligned skewness | Mean P_BT |",
               "|---|---:|---:|---:|---:|"])
    table_rows = []
    cards = []
    for name, item in metrics["cases"].items():
        numbers = [fmt(item[key]) for key in ("final_barotropic_kinetic", "final_lowest_barotropic_shell_kinetic",
                                              "final_quarter_f_aligned_skewness_mean", "final_quarter_barotropic_power_mean")]
        md.append(f"| {item['title']} | " + " | ".join(numbers) + " |")
        table_rows.append("<tr><th>" + html.escape(item["title"]) + "</th>" + "".join(f"<td>{html.escape(x)}</td>" for x in numbers) + "</tr>")
        movie = item["movie"]
        cards.append(f'''<section class="movie" id="{html.escape(name)}"><h2>{html.escape(item['title'])}</h2>
<p>Cutoff K<sub>h</sub>={item['K']}, N<sub>v</sub>={item['Nv']}; G={item['G']:.6g};
three actual heights z/H={', '.join(f'{z:g}' for z in item['saved_heights'])}.</p>
<video controls playsinline preload="metadata" poster="{html.escape(movie['poster'])}">
<source src="{html.escape(movie['file'])}" type="video/mp4">Your browser cannot display this video.</video>
<p class="small">Top: vertical vorticity. Bottom: buoyancy perturbation, not passive tracer.
{item['frames']} frames, {movie['fps']:g} frames/s. <a href="{html.escape(movie['file'])}">Open movie</a>.</p></section>''')
    md.extend(["", "![Diagnostic overview](diagnostics.png)", "", "## Movies", ""])
    for name, item in metrics["cases"].items():
        md.append(f"- [{item['title']}]({item['movie']['file']})")
    md.extend(["", "Fixed color scales and movie verification are recorded in `publication-metrics.json`. "
               "Cumulative third-moment production uses saved-sample trapezoidal quadrature, not an RK-stage budget test.", ""])
    (output / "README.md").write_text("\n".join(md))
    scale_text = " ".join(f"{key}: full physical range [{scale['physical_min']:.4g}, {scale['physical_max']:.4g}], "
                          f"asinh linear width {scale['linear_width']:.4g}." for key, scale in scales.items())
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Nonlinear Galerkin convection · Jupiter flow experiments</title>
<style>body{{font:16px/1.55 system-ui,sans-serif;color:#192b3f;background:#f5f7fa;margin:0}}main{{max-width:1160px;margin:auto;padding:28px 22px 60px}}h1{{font-size:2rem;line-height:1.2}}h2{{font-size:1.25rem}}a{{color:#165caa}}.notice,.movie,.overview{{background:white;border:1px solid #d7e0e9;border-radius:10px;padding:20px;margin:24px 0}}.small{{font-size:.88rem;color:#516172}}video,img{{display:block;width:100%;height:auto;border-radius:5px}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}td,th{{text-align:right;padding:9px;border-bottom:1px solid #d7e0e9}}th:first-child{{text-align:left}}.tablewrap{{overflow-x:auto}}nav{{margin:15px 0}}nav a{{margin-right:14px}}code{{font-size:.9em}}</style></head>
<body><main><a href="https://npizzo.github.io/jupiter-ffr-breeze/">← Gallery landing page</a>
<h1>Fully coupled nonlinear Galerkin convection</h1>
<p>Completed finite-mode reference tests · dimensionless time 0–{metrics['config']['duration']:g}</p>
<div class="notice">{''.join('<p>' + html.escape(text) + '</p>' for text in paragraphs)}</div>
<nav>{''.join(f'<a href="#{html.escape(name)}">{html.escape(item["title"])}</a>' for name, item in metrics['cases'].items())}</nav>
<section class="overview"><h2>Transfer, sign and numerical checks</h2>
<img src="diagnostics.png" alt="Five-case barotropic energy, aligned skewness, low-shell energy, cubic-moment production, stress power and integrated budget defects">
<p class="small">Solid cumulative-production curves: convective stress. Dashed: projected barotropic self-advection.
Colors match cases in the legend. Cumulative third-moment terms use diagnostic-sample quadrature; they are not RK-stage closure tests.</p>
<div class="tablewrap"><table><thead><tr><th>Case</th><th>Final E<sub>BT</sub></th><th>Final low-shell E<sub>BT</sub></th><th>Late aligned skewness</th><th>Late mean P<sub>BT</sub></th></tr></thead><tbody>{''.join(table_rows)}</tbody></table></div></section>
{''.join(cards)}
<section class="notice"><h2>Scales and provenance</h2><p>{html.escape(scale_text)}</p>
<p>The same mapping is fixed across all movies, times and depths; no framewise normalization or field clipping.
Saved 64×64 samples are exact retained-mode reconstructions, enlarged for display. They do not imply a 64³ dynamical calculation.
Each movie passed frame-count, duration and full-decode checks on Unity.</p>
<p><a href="publication-metrics.json">Metrics, configuration, checksums and media verification</a> ·
<a href="diagnostics.pdf">Diagnostic figure PDF</a> ·
<a href="https://github.com/npizzo/jupiter-ffr-breeze/tree/main/experiments/galerkin-convection-v1">Equations and experiment plan</a> ·
<a href="https://github.com/npizzo/jupiter-ffr-breeze/tree/main/scripts/galerkin_convection_v1">Implementation</a> ·
<a href="https://arxiv.org/abs/2505.07678">Mode-selection reference</a></p></section>
</main></body></html>'''
    (output / "index.html").write_text(document)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("render_output", type=Path)
    parser.add_argument("--fps", type=float, default=12.0)
    args = parser.parse_args()
    if args.fps <= 0 or not math.isfinite(args.fps):
        raise ValueError("fps must be positive and finite")
    for executable in ("ffmpeg", "ffprobe"):
        if not shutil.which(executable):
            raise RuntimeError(f"Required Unity publication tool missing: {executable}")
    # Reuse permits a failed render to be retried, but completion is written last.
    args.render_output.mkdir(parents=True, exist_ok=True)
    marker = args.render_output / "PUBLICATION_READY.json"
    if marker.exists():
        raise RuntimeError("Publication already verified; use a new output directory to revise it")
    started = time.monotonic()
    config, cases = load_completed(args.campaign_root)
    metrics = calculate_metrics(config, cases)
    scales = color_scales(cases)
    metrics["color_scales"] = scales
    metrics["analysis_source_sha256"] = checksum(__file__)
    metrics["validation_record_sha256"] = checksum(args.campaign_root / "validation" / "PASS.json")
    make_scientific_plot(cases, args.render_output)
    for case in cases:
        print(f"Rendering {case['name']} from completed saved output", flush=True)
        metrics["cases"][case["name"]]["movie"] = render_movie(case, scales, args.render_output, args.fps)
    metrics["analysis_render_elapsed_seconds"] = time.monotonic() - started
    write_json(args.render_output / "publication-metrics.json", metrics)
    write_publication(args.render_output, metrics, scales)
    write_json(marker, {"status": "verified_finished_assets", "cases": [case["name"] for case in cases],
                        "metrics_sha256": checksum(args.render_output / "publication-metrics.json"),
                        "index_sha256": checksum(args.render_output / "index.html"),
                        "elapsed_seconds": time.monotonic() - started,
                        "uploaded": False})
    print(json.dumps({"publication_ready": str(args.render_output),
                      "elapsed_seconds": time.monotonic() - started}), flush=True)


if __name__ == "__main__":
    main()
