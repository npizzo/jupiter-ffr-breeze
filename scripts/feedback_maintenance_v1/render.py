"""Render the completed two-case maintenance intervention on Unity.

Usage: python render.py CAMPAIGN_ROOT NEW_PUBLICATION_DIRECTORY
Requires CAMPAIGN_ROOT/{full,no_low_quadratic} completed through t=450.
"""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import sys
import tomllib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg


CASES = ("full", "no_low_quadratic")
LABELS = {"full": "Full quadratic + cubic feedback",
          "no_low_quadratic": "Quadratic output removed only at 0 < |k| ≤ 2"}
TERMS = ("advection", "linear", "quadratic", "cubic", "viscous")
N = 512
TIMES = np.arange(400., 450.5, .5)
STEM = "maintenance_controls_512_400_450"


def field(path):
    a = np.fromfile(path, dtype=np.float32)
    assert a.size == N * N and np.isfinite(a).all(), str(path)
    return a.reshape((N, N), order="F")


def run(root, out):
    out.mkdir(parents=True, exist_ok=False)
    frames, scalars, configs = {}, {}, {}
    for case in CASES:
        p = root / case
        c = tomllib.loads((p / "config.toml").read_text())
        s = tomllib.loads((p / "status.toml").read_text())
        assert s["status"] == "COMPLETED" and s["t"] == s["target"] == 450.
        assert c["case"] == case and c["start_time"] == 400. and c["target_time"] == 450.
        assert c["grid"] == [N, N] and np.allclose(c["domain"], [2*np.pi, 2*np.pi])
        assert c["a"] == .0002 and c["b"] == .002 and c["c"] == .0002
        assert c["f"] == 1. and c["beta"] == 0 and c["viscosity"] == .0001
        assert c["forcing_band"] == [12., 16.] and not c["allow_zonal_forcing"]
        assert c["velocity_drag"] == c["hyperviscosity"] == 0
        assert c["cubic_unchanged"] and c["remove_quadratic_output_k_le_2"] == (case == "no_low_quadratic")
        rows = list(csv.DictReader((p / "frames.csv").open()))
        assert np.array_equal([float(x["t"]) for x in rows], TIMES)
        frames[case] = [p / x["file"] for x in rows]
        a = np.genfromtxt(p / "scalars.csv", delimiter=",", names=True)
        assert np.array_equal(a["t"], np.arange(400., 450.25, .25))
        for name in a.dtype.names:
            assert np.isfinite(a[name]).all(), (case, name)
        # The intervention removes this source only; no diagnostic assumption
        # that the remaining advection or cubic power is small is imposed.
        assert np.max(np.abs(a["linear_lowk_energy_rate"])) < 1e-12
        if case == "no_low_quadratic":
            assert np.max(np.abs(a["quadratic_lowk_energy_rate"])) < 1e-12
        assert np.allclose(a["lowk_energy_change"], a["energy_k_le_2"]-a["energy_k_le_2"][0], rtol=1e-12, atol=1e-13)
        integrated = sum(a[f"integral_{k}_lowk_work"] for k in TERMS)
        assert np.allclose(a["sampled_lowk_budget_residual"], a["lowk_energy_change"]-integrated, rtol=1e-11, atol=1e-13)
        configs[case], scalars[case] = c, a
    assert configs["full"]["source_checkpoint_sha256"] == configs["no_low_quadratic"]["source_checkpoint_sha256"]
    assert configs["full"]["code_sha256"] == configs["no_low_quadratic"]["code_sha256"]
    assert np.array_equal(field(frames["full"][0]), field(frames["no_low_quadratic"][0]))

    # Recover the native inverse-Laplacian energy and source work independently
    # from three saved Float32 fields, using the known C-grid derivative symbol.
    k = np.fft.fftfreq(N, d=1/N)
    k2 = k[:, None]**2 + k[None, :]**2
    dx = 2*np.pi/N
    lam = (4*np.sin(k[:, None]*dx/2)**2 + 4*np.sin(k[None, :]*dx/2)**2) / dx**2
    invlam = np.divide(1., lam, out=np.zeros_like(lam), where=lam != 0)
    low = (k2 > 0) & (k2 <= 4)
    reconstructed_checks = {}
    for case in CASES:
        error = 0.
        for i in (0, 50, 100):
            q = field(frames[case][i]).astype(np.float64)
            qh = np.fft.fft2(q);a = scalars[case][2*i]
            energy = .5*np.sum(np.abs(qh[low])**2 * invlam[low]) / N**4
            assert np.isclose(energy, a["energy_k_le_2"], rtol=2e-5, atol=2e-9)
            for term, source in (("quadratic", .002*q*q), ("cubic", -.0002*q*q*q)):
                sh = np.fft.fft2(source-source.mean())
                if case == "no_low_quadratic" and term == "quadratic":
                    sh[low] = 0.
                value = np.sum((qh.conj()*sh).real[low]*invlam[low])/N**4
                recorded = a[f"{term}_lowk_energy_rate"]
                assert np.isclose(value, recorded, rtol=2e-5, atol=2e-9), (case, i, term)
                error = max(error, float(abs(value-recorded)))
        reconstructed_checks[case] = {"sampled_frame_indices": [0, 50, 100], "maximum_source_work_absolute_difference": error}

    limit = max(float(np.abs(field(p)).max()) for paths in frames.values() for p in paths)
    assert limit > 0
    norm = AsinhNorm(linear_width=.1, vmin=-limit, vmax=limit)
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    plt.rcParams["animation.ffmpeg_path"] = ff
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.7), layout="constrained")
    images = []
    for ax, case in zip(axes, CASES):
        im = ax.imshow(field(frames[case][0]).T, origin="lower", extent=(0, 2*np.pi, 0, 2*np.pi),
                       cmap="RdBu_r", norm=norm, interpolation="nearest")
        ax.set(title=LABELS[case], xlabel="x", ylabel="y")
        images.append(im)
    fig.colorbar(images[0], ax=list(axes), shrink=.8, label="Vertical vorticity ζ · shared fixed asinh scale")
    title = fig.suptitle("")
    fig.supxlabel("512² · identical mature t=400 initial state · f=+1 · horizontal 2D, no tracer or depths\n"
                  "Maintenance test, not spontaneous formation · positive ζ is cyclonic", fontsize=10)
    writer = FFMpegWriter(fps=6, codec="libx264", extra_args=["-b:v", "2200k", "-maxrate", "2600k",
                           "-bufsize", "5200k", "-threads", "1", "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    movie = out / (STEM + ".mp4")
    with writer.saving(fig, str(movie), dpi=100):
        for i, t in enumerate(TIMES):
            for im, case in zip(images, CASES):
                im.set_data(field(frames[case][i]).T)
            title.set_text(f"What maintains the large vortices? · t = {t:g} · elapsed since intervention = {t-400:g}")
            writer.grab_frame()
            if i == len(TIMES)-1:
                fig.savefig(out / (STEM + ".png"), dpi=100)
    plt.close(fig)
    count, duration = imageio_ffmpeg.count_frames_and_secs(str(movie))
    assert count == 101 and abs(duration-101/6) < .1
    subprocess.run([ff, "-v", "error", "-xerror", "-i", str(movie), "-f", "null", "-"], check=True)
    assert movie.stat().st_size < 25*1024**2
    print("MEDIA_VERIFIED", movie.name, flush=True)

    fig, axs = plt.subplots(3, 2, figsize=(14, 12), layout="constrained")
    for case in CASES:
        a = scalars[case]
        axs[0, 0].plot(a["t"], a["energy_k_le_2"], label=LABELS[case])
        axs[0, 1].plot(a["t"], a["skewness"], label=LABELS[case])
    axs[0, 0].set(title="Absolute energy in 0 < |k| ≤ 2", ylabel="Energy")
    axs[0, 1].set(title="Full-field vorticity skewness", ylabel="Skewness")
    axs[0, 0].legend(fontsize=8)
    colors = {"advection": "#3577b5", "linear": "#555555", "quadratic": "#d57922", "cubic": "#ad4b88", "viscous": "#398259"}
    for j, case in enumerate(CASES):
        a = scalars[case]
        for k in TERMS:
            axs[1, j].plot(a["t"], a[f"{k}_lowk_energy_rate"], label=k, color=colors[k])
            axs[2, j].plot(a["t"], a[f"integral_{k}_lowk_work"], label=k, color=colors[k])
        axs[1, j].set(title=LABELS[case], ylabel="Native instantaneous low-mode power")
        axs[2, j].plot(a["t"], a["lowk_energy_change"], color="black", lw=2, label="Measured energy change")
        axs[2, j].plot(a["t"], a["sampled_lowk_budget_residual"], color="#df2626", ls="--", lw=2, label="Unassigned sampling residual")
        axs[2, j].set(title="Cumulative sampled work vs measured change", ylabel="Energy change / integrated work")
        axs[1, j].legend(fontsize=8);axs[2, j].legend(fontsize=8)
    for ax in axs.flat:
        ax.axhline(0, color="gray", lw=.5);ax.grid(alpha=.2);ax.set(xlabel="Dimensionless model time")
    fig.suptitle("Matched maintenance intervention · native work sampled every 0.25 time units\n"
                 "Trapezoidal cumulative work is not an exact every-stage numerical budget", fontsize=14)
    fig.savefig(out / "maintenance_lowmode_work.png", dpi=150);plt.close(fig)

    summary = {}
    for case in CASES:
        a = scalars[case]
        integral_values = {k: float(a[f"integral_{k}_lowk_work"][-1]) for k in TERMS}
        sum_abs_integrals = sum(abs(v) for v in integral_values.values())
        final_residual = float(a["sampled_lowk_budget_residual"][-1])
        # Report residual size explicitly without inventing a physical pass/fail
        # threshold. The reader can see cancellation and sampling limitations.
        summary[case] = {
            "initial": {k: float(a[k][0]) for k in a.dtype.names},
            "final": {k: float(a[k][-1]) for k in a.dtype.names},
            "low_energy_change": float(a["lowk_energy_change"][-1]),
            "low_energy_final_over_initial": float(a["energy_k_le_2"][-1] / a["energy_k_le_2"][0]),
            "sampled_cumulative_lowmode_work": integral_values,
            "sampled_budget_residual_final": final_residual,
            "sampled_budget_residual_maximum_absolute": float(np.abs(a["sampled_lowk_budget_residual"]).max()),
            "final_residual_over_sum_absolute_final_term_integrals": final_residual / max(sum_abs_integrals, 1e-30),
            "final_residual_over_initial_low_energy": final_residual / a["energy_k_le_2"][0],
            "maximum_absolute_quadratic_lowmode_power": float(np.abs(a["quadratic_lowk_energy_rate"]).max()),
        }
        (out / (case + "_scalars.csv")).write_bytes((root / case / "scalars.csv").read_bytes())
    public = {
        "status": "COMPLETED_MEDIA_VERIFIED", "compute_location": "Unity", "grid": [N, N],
        "domain": [float(2*np.pi), float(2*np.pi)], "time_range": [400, 450], "elapsed_range": [0, 50],
        "time_units": "dimensionless", "matched_initial_fields": True,
        "parameters": {"f": 1., "beta": 0., "a": .0002, "b": .002, "c": .0002, "viscosity": .0001,
                       "forcing_band": [12, 16], "allow_zonal_forcing": False, "drag": 0., "hyperviscosity": 0.},
        "intervention": "Remove only the quadratic source output in 0<|k|≤2; cubic source and all other terms remain unchanged",
        "movie": {"file": movie.name, "frames": count, "frames_per_second": 6, "duration_seconds": duration,
                  "full_decode_passed": True, "layout": "1 row × 2 synchronized horizontal 2D vorticity panels",
                  "frame_interval": .5, "tracer": False, "depth_layers": False},
        "color_scale": {"type": "symmetric asinh", "width": .1, "limits": [-limit, limit],
                        "fixed_across_cases_and_time": True, "clipping": False},
        "diagnostic_interval": .25, "summary": summary, "independent_saved_field_checks": reconstructed_checks,
        "limitations": [
            "This tests maintenance from a mature shared state, not spontaneous formation or independence from initial conditions.",
            "The imposed polynomial closure is phenomenological, not a derived or validated convection closure.",
            "Direct low-mode stress drive can be physically legitimate; dependence on it would identify this model's mechanism, not falsify convective LSV formation.",
            "Native instantaneous tendency work is sampled every 0.25 model time units. Cumulative trapezoidal work is not exact every-stage accounting; the unassigned residual is explicitly reported and plotted.",
            "Skewness and low-mode energy do not alone establish a single LSV or an inverse cascade. One initial state and one resolution do not establish robustness."
        ],
    }
    (out / "report.json").write_text(json.dumps(public, indent=2, allow_nan=False) + "\n")
    assets = {p.name: {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
              for p in out.iterdir() if p.is_file()}
    (out / "transport.json").write_text(json.dumps({"source_root": str(root), "configs": configs,
        "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "assets": assets}, indent=2) + "\n")
    print("PUBLICATION_READY", flush=True)


if __name__ == "__main__":
    assert len(sys.argv) == 3, __doc__
    run(Path(sys.argv[1]), Path(sys.argv[2]))
