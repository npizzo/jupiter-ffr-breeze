"""Render one case or a matched pair with a fixed, shared vorticity scale."""
import argparse
import json
from pathlib import Path
import tomllib

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, writers
import numpy as np


def read_frame(path, coriolis):
    with h5py.File(path) as stream:
        t = float(stream["t"][()])
        zeta = stream["zeta"][:] / coriolis
    if not np.isfinite(zeta).all() or not np.isfinite(t):
        raise ValueError(f"Nonfinite frame: {path.name}")
    return t, zeta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cases", nargs="+", type=Path)
    p.add_argument("--output", required=True, type=Path, help="MP4 file; also writes a poster and JSON record")
    p.add_argument("--fps", type=int, default=6, help="Saved frames per video second")
    p.add_argument("--allow-sparse", action="store_true", help="Explicitly labeled preview of incomplete/sparse data")
    args = p.parse_args()
    if not 1 <= len(args.cases) <= 2 or args.fps < 1 or args.output.suffix != ".mp4":
        p.error("Supply one or two cases, positive fps, and an .mp4 output")
    if not writers.is_available("ffmpeg"):
        raise RuntimeError("Install system ffmpeg with the libx264 encoder to render MP4 videos")
    cases = []
    limit = 0.0
    for folder in args.cases:
        m = tomllib.loads((folder / "metadata.toml").read_text())
        if not np.isfinite(m["f"]) or m["f"] == 0:
            raise ValueError("A finite, nonzero f is required")
        if not args.allow_sparse and m["status"] != "COMPLETE":
            raise ValueError(f"Case not complete: {folder.name}")
        files = sorted(folder.glob("frame_*.jld2"))
        times = []
        for file in files:
            t, field = read_frame(file, m["f"])
            if field.shape != (m["ny"], m["nx"]):
                raise ValueError(f"Unexpected map shape: {file.name}")
            times.append(t)
            limit = max(limit, float(np.max(np.abs(field))))
        times = np.asarray(times)
        if len(times) < 2 or not np.all(np.diff(times) > 0):
            raise ValueError("Need at least two frames with strictly increasing times")
        if not args.allow_sparse:
            expected = np.arange(0, m["stop_seconds"] + .5 * m["output_seconds"], m["output_seconds"])
            if times.shape != expected.shape or not np.allclose(times, expected, rtol=0, atol=1e-6):
                raise ValueError("Missing frames or inconsistent time spacing; retrieve every saved frame")
        cases.append((m, files, times))
    if len(cases) == 2:
        keys = ("nx", "ny", "nz", "lx", "ly", "height", "diagnostic_height_m", "f", "upper_cooling_flux")
        if not np.array_equal(cases[0][2], cases[1][2]) or any(cases[0][0][k] != cases[1][0][k] for k in keys):
            raise ValueError("Paired video requires matching grids, map heights, f, cooling and frame times")
    limit = max(limit, 1e-12)
    fig, axes = plt.subplots(1, len(cases), figsize=(6 * len(cases), 5), squeeze=False, layout="constrained")
    artists = []
    for ax, (m, files, _) in zip(axes.flat, cases):
        _, field = read_frame(files[0], m["f"])
        k = np.clip(round(m.get("diagnostic_height", m["cap_height"]) / m["height"] * m["nz"]), 1, m["nz"])
        z_km = (k - .5) * m["height"] / m["nz"] / 1000
        im = ax.imshow(field, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit,
                       extent=(0, m["lx"] / 1000, 0, m["ly"] / 1000), interpolation="nearest")
        label = "Heated + cooled" if m["heat_flux"] > 0 else "Cooling only"
        ax.set(title=f"{label}: bottom {m['heat_flux']:g} W/m²\nz = {z_km:.2f} km", xlabel="x [km]", ylabel="y [km]")
        artists.append(im)
    fig.colorbar(artists[0], ax=list(axes.flat), shrink=.8, label="Vertical vorticity ζ/f (fixed scale)")
    prefix = "Sparse preview | " if args.allow_sparse else ""
    title = fig.suptitle("")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=args.fps, codec="libx264", extra_args=["-crf", "24", "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    with writer.saving(fig, str(args.output), dpi=100):
        for index, t in enumerate(cases[0][2]):
            for artist, (m, files, _) in zip(artists, cases):
                _, field = read_frame(files[index], m["f"])
                artist.set_data(field)
            title.set_text(f"{prefix}Vorticity evolution | {t/3600:.1f} Earth hours")
            writer.grab_frame()
    fig.savefig(args.output.with_suffix(".png"), dpi=100)
    plt.close(fig)
    record = {"cases": [m["name"] for m, _, _ in cases], "frames": len(cases[0][2]),
              "times_seconds": cases[0][2].tolist(), "frames_per_second": args.fps,
              "color_limits_zeta_over_f": [-limit, limit], "normalization": "one shared maximum over every frame and case; no clipping",
              "sparse_preview": args.allow_sparse, "poster": "final frame",
              "note": "Horizontal dry-model vorticity; not cloud imagery. No temporal interpolation."}
    args.output.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Rendered {record['frames']} frames: {args.output}")


if __name__ == "__main__":
    main()
