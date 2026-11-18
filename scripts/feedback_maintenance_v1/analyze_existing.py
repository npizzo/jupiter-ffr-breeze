"""Analyze completed polynomial-feedback runs on Unity, without advancing them.

Usage: python analyze_existing.py CAMPAIGN_ROOT NEW_OUTPUT_DIRECTORY
Only publication plots and scalar reports are written. Raw fields stay on Unity.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import tomllib

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


N = 512
L = 2 * np.pi
DX = L / N
NU = 1e-4
CASES = ("cubic_control", "biased_plus", "biased_minus")
K = np.fft.fftfreq(N, d=1 / N)
K2 = K[:, None] ** 2 + K[None, :] ** 2
LAM = (4 * np.sin(K[:, None] * DX / 2) ** 2
       + 4 * np.sin(K[None, :] * DX / 2) ** 2) / DX ** 2
INVLAM = np.divide(1., LAM, out=np.zeros_like(LAM), where=LAM != 0)
LOW = (K2 > 0) & (K2 <= 4)
SHELL = (K2 >= 12 ** 2) & (K2 <= 16 ** 2)
SHELL[0, :] = False
NORM = N ** 4


def velocity(q):
    """Native C-grid inverse curl with zero mean velocity."""
    psi = np.fft.ifft2(-np.fft.fft2(q) * INVLAM).real
    return -(np.roll(psi, -1, axis=1) - psi) / DX, (np.roll(psi, -1, axis=0) - psi) / DX


def curl(u, v):
    return (v - np.roll(v, 1, axis=0) - u + np.roll(u, 1, axis=1)) / DX


def centered_advection(u, v):
    """Second-order conservative C-grid momentum-flux reconstruction.

    This has independent algebraic checks, but no comparison against the running
    Oceananigans native tendency; report it as a reconstruction, not an exact
    native budget. Pressure contributes no work to the reconstructed solenoidal
    low-mode velocity in the periodic domain.
    """
    ux = .5 * (u + np.roll(u, -1, axis=0))
    vy = .5 * (v + np.roll(v, -1, axis=1))
    cross = .25 * (u + np.roll(u, 1, axis=1)) * (v + np.roll(v, 1, axis=0))
    au = -(ux * ux - np.roll(ux * ux, 1, axis=0)
           + np.roll(cross, -1, axis=1) - cross) / DX
    av = -(np.roll(cross, -1, axis=0) - cross
           + vy * vy - np.roll(vy * vy, 1, axis=1)) / DX
    return au, av


def work(qh, source, mask=None):
    w = (qh.conj() * np.fft.fft2(source)).real * INVLAM / NORM
    return float(w.sum() if mask is None else w[mask].sum())


def components(q, sign, threshold):
    """Four-neighbor connectivity including BOTH periodic seams."""
    mask = sign * q > threshold
    labs, count = ndimage.label(mask)
    parent = np.arange(count + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def unite(a, b):
        if a and b:
            a, b = find(a), find(b)
            if a != b:
                parent[b] = a

    for a, b in zip(labs[0], labs[-1]):
        unite(int(a), int(b))
    for a, b in zip(labs[:, 0], labs[:, -1]):
        unite(int(a), int(b))
    roots = np.array([find(i) for i in range(count + 1)])
    labels = roots[labs]
    area = np.bincount(labels.ravel(), minlength=count + 1)
    circulation = np.bincount(labels.ravel(), weights=(sign * q).ravel(), minlength=count + 1)
    ids = np.flatnonzero(area >= 4)
    ids = ids[ids != 0]
    ids = ids[np.argsort(circulation[ids])[::-1]]
    positive_total = np.maximum(sign * q, 0).sum()
    return {
        "count_at_least_four_cells": int(len(ids)),
        "count_each_carrying_at_least_five_percent_positive_circulation": int(
            np.sum(circulation[ids] >= .05 * positive_total)),
        "largest_area_fraction": float(area[ids].max() / q.size) if len(ids) else 0.,
        "largest_core_circulation_fraction": float(circulation[ids[0]] / positive_total) if len(ids) else 0.,
        "leading_core_area_fractions": [float(area[i] / q.size) for i in ids[:6]],
        "leading_core_circulation_fractions": [float(circulation[i] / positive_total) for i in ids[:6]],
    }


def self_checks():
    x = np.arange(N)[:, None] * DX
    y = np.arange(N)[None, :] * DX
    q = np.sin(2*x + 3*y) + .6*np.cos(4*x-y) + .3*np.sin(x+2*y)
    u, v = velocity(q)
    au, av = centered_advection(u, v)
    assert np.max(np.abs(curl(u, v) - q)) < 1e-10
    assert abs(np.mean(u*au + v*av)) < 1e-13
    qh = np.fft.fft2(q)
    assert abs(work(qh, curl(au, av)) - np.mean(u*au + v*av)) < 1e-13
    assert abs(np.mean(au)) + abs(np.mean(av)) < 1e-13
    # A core crossing both periodic seams must be ONE connected component.
    test = np.zeros((16, 16));test[[0, -1], :] = 0
    test[np.ix_([0, 1, 14, 15], [0, 1, 14, 15])] = 3
    c = components(test, 1, 2)
    assert c["count_at_least_four_cells"] == 1
    assert c["largest_area_fraction"] == 16/256
    return {"inverse_curl": "PASS", "centered_advection_energy_conservation": "PASS",
            "curl_vs_momentum_work": "PASS", "periodic_connectivity": "PASS"}


def inventory(root, case):
    frames, scalars, provenance = {}, {}, []
    for end in (80, 160, 240, 320, 400):
        p = root / case / f"to_{end:03d}"
        config = tomllib.loads((p / "config.toml").read_text())
        state = tomllib.loads((p / "status.toml").read_text())
        assert state["status"] == "COMPLETED" and state["t"] == end
        assert config["grid"] == [N, N] and config["viscosity"] == NU
        assert config["a"] == .0002 and config["c"] == .0002 and config["beta"] == 0
        assert config["b"] == (0 if case == "cubic_control" else .002)
        assert config["forcing_band"] == [12., 16.] and not config["allow_zonal_forcing"]
        assert config["f"] == (-1 if case == "biased_minus" else 1)
        for row in csv.DictReader((p / "frames.csv").open()):
            t = float(row["t"])
            path = p / row["file"]
            if t in frames:
                assert path.read_bytes() == frames[t].read_bytes()
            frames[t] = path
        for row in csv.DictReader((p / "scalars.csv").open()):
            row = {k: float(v) for k, v in row.items()}
            t = row["t"]
            if t in scalars:
                for k in ("energy", "enstrophy", "energy_k_le_2", "quadratic_lowk_energy_rate", "cubic_lowk_energy_rate"):
                    assert np.isclose(row[k], scalars[t][k], rtol=1e-10, atol=1e-12)
            scalars[t] = row
        provenance.append({"segment_end": end, "config_sha256": hashlib.sha256((p / "config.toml").read_bytes()).hexdigest(),
                           "scalar_sha256": hashlib.sha256((p / "scalars.csv").read_bytes()).hexdigest()})
    assert sorted(frames) == list(np.arange(0, 400.5, .5)) == sorted(scalars)
    return frames, scalars, config, provenance


def load_field(p):
    q = np.fromfile(p, dtype=np.float32)
    assert q.size == N*N and np.isfinite(q).all()
    return q.reshape((N, N), order="F").astype(np.float64)


def analyze(root, out):
    out.mkdir(parents=True, exist_ok=False)
    validation = self_checks()
    results, histories = {}, {}
    for case in CASES:
        frames, scalars, config, provenance = inventory(root, case)
        initial = load_field(frames[0.])
        threshold = 2 * np.sqrt(np.mean(initial ** 2))
        rows, morphology = [], []
        max_errors = dict(energy=0., enstrophy=0., low_energy=0., quadratic_work=0., cubic_work=0.)
        for t in np.arange(0, 405, 5.):
            q = load_field(frames[t]);qh = np.fft.fft2(q);s = scalars[t]
            u, v = velocity(q)
            e_modes = .5 * np.abs(qh) ** 2 * INVLAM / NORM
            energy, lowe = float(e_modes.sum()), float(e_modes[LOW].sum())
            q2 = config["b"]*config["f"]*q*q;q2 -= q2.mean()
            q3 = -config["c"]*q*q*q;q3 -= q3.mean()
            i2, i3 = work(qh, q2, LOW), work(qh, q3, LOW)
            for name, a, b in (("energy", energy, s["energy"]), ("enstrophy", .5*np.mean(q*q), s["enstrophy"]),
                               ("low_energy", lowe, s["energy_k_le_2"]), ("quadratic_work", i2, s["quadratic_lowk_energy_rate"]),
                               ("cubic_work", i3, s["cubic_lowk_energy_rate"])):
                error = abs(a-b) / max(abs(b), 1e-8)
                max_errors[name] = max(max_errors[name], float(error))
                assert np.isclose(a, b, rtol=2e-5, atol=2e-9), (case, t, name, a, b)
            au, av = centered_advection(u, v)
            adv = work(qh, curl(au, av), LOW)
            visc = float(-NU * np.sum(np.abs(qh[LOW])**2) / NORM)
            lin = float(config["a"]*np.sum(np.abs(qh[LOW])**2 * K2[LOW]*SHELL[LOW]*INVLAM[LOW]) / NORM)
            assert lin == 0.
            kmean = float(np.sum(np.sqrt(K2)*e_modes)/max(energy, 1e-30))
            lam_energy = float(2*np.pi*np.sum(np.divide(e_modes, np.sqrt(K2), out=np.zeros_like(e_modes), where=K2>0))/max(energy, 1e-30))
            urms = float(np.sqrt(2*energy))
            rate = adv + i2 + i3 + visc + lin
            m = components(q, np.sign(config["f"]), threshold)
            row = dict(t=float(t), energy=energy, energy_k_le_2=lowe, low_energy_fraction=lowe/energy,
                       quadratic_lowk_power=i2, cubic_lowk_power=i3, reconstructed_advection_lowk_power=adv,
                       viscosity_lowk_power=visc, linear_lowk_power=lin, reconstructed_total_lowk_power=rate,
                       reconstructed_advection_total_power=float(np.mean(u*au + v*av)),
                       energy_weighted_wavelength=lam_energy, energy_weighted_mean_wavenumber=kmean,
                       urms=urms, domain_crossing_time=L/max(urms, 1e-30),
                       integral_scale_crossing_time=lam_energy/max(urms, 1e-30),
                       skewness=s["skewness"], largest_core_area_fraction=m["largest_area_fraction"],
                       major_core_count=m["count_each_carrying_at_least_five_percent_positive_circulation"])
            rows.append(row);morphology.append(dict(t=float(t), **m))
        # All original scalar output is used for time-integrating source terms;
        # the advective reconstruction uses the deliberately sparse five-unit sample.
        st = np.array(sorted(scalars)); sr = [scalars[t] for t in st]
        summary = {}
        for start in (0., 300., 350.):
            sampled = [r for r in rows if r["t"] >= start]
            t = np.array([r["t"] for r in sampled]);sel=st>=start
            integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
            integral = {k: float(integrate([r[k] for r in sampled], t)) for k in
                        ("reconstructed_advection_lowk_power", "viscosity_lowk_power")}
            integral["quadratic_lowk_power"] = float(integrate([r["quadratic_lowk_energy_rate"] for r in sr if r["t"]>=start], st[sel]))
            integral["cubic_lowk_power"] = float(integrate([r["cubic_lowk_energy_rate"] for r in sr if r["t"]>=start], st[sel]))
            change = scalars[400.]["energy_k_le_2"]-scalars[start]["energy_k_le_2"]
            summary[f"{start:g}_400"] = {"low_energy_change": change, "sampled_term_integrals": integral,
                                       "unassigned_difference_including_sampling_and_reconstruction": change-sum(integral.values())}
        history_file = out / f"{case}_sampled.csv"
        with history_file.open("w") as fp:
            w=csv.DictWriter(fp, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        results[case] = {"final": rows[-1], "windows": summary, "scalar_reconstruction_max_relative_errors": max_errors,
                         "core_threshold": threshold, "core_threshold_rule": "2 times original initial RMS vorticity, fixed in time",
                         "periodic_core_samples": morphology, "provenance": provenance}
        histories[case]=rows
        print("ANALYZED", case, "t400", json.dumps(rows[-1]), flush=True)
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for case,rows in histories.items():
        tt=[r["t"] for r in rows]
        for ax,key,label in zip(axs.flat,("energy_k_le_2","energy_weighted_wavelength","largest_core_area_fraction","major_core_count"),
                                ("Absolute energy in |k|≤2","Energy-weighted wavelength","Largest threshold-core area / domain","Cores holding ≥5% of signed-positive circulation")):
            ax.plot(tt,[r[key] for r in rows],label=case);ax.set(xlabel="Dimensionless time",title=label);ax.grid(alpha=.25)
    axs[0,0].legend();fig.suptitle("Existing polynomial model: large-scale organization and periodic core counts")
    fig.savefig(out/"existing_organization.png",dpi=150);plt.close(fig)
    rows=histories["biased_plus"];fig,axs=plt.subplots(2,1,figsize=(12,8),layout="constrained")
    for key,label in (("reconstructed_advection_lowk_power","Centered momentum-advection reconstruction"),
                      ("quadratic_lowk_power","Quadratic source"),("cubic_lowk_power","Cubic source"),
                      ("viscosity_lowk_power","Molecular viscosity")):
        for ax in axs:
            ax.plot([r["t"] for r in rows],[r[key] for r in rows],label=label)
    axs[0].set(title="Low-mode energy supply: biased +f");axs[1].set(xlim=(300,400),title="Late-time window")
    for ax in axs:ax.axhline(0,color="black",lw=.7);ax.grid(alpha=.2);ax.set(xlabel="Dimensionless time",ylabel="Power into |k|≤2")
    axs[0].legend();fig.suptitle("Saved-state reconstruction, NOT a closed all-stage numerical budget")
    fig.savefig(out/"existing_lowmode_power.png",dpi=150);plt.close(fig)
    public={"status":"ANALYSIS_COMPLETE","compute_location":"Unity","cases":results,"validation":validation,
            "frame_sample_interval":5.,"scalar_interval":.5,"time_range":[0,400],"low_mode_definition":"0<|k|≤2",
            "limitations":["Native staggered inverse-curl energy and nonlinear work reproduce saved scalar data within Float32 field precision.",
             "Advection is an independently energy-conserving centered C-grid momentum reconstruction; it has not been compared to the native Oceananigans tendency.",
             "Term integrals use trapezoidal time sampling, with coarser sampling for advection and viscosity. The unassigned difference is retained, not labeled physical forcing.",
             "Connected cores use one declared fixed vorticity threshold, four-neighbor periodic connectivity and a four-cell minimum; counts are descriptive, not an objective LSV theorem.",
             "One initial realization per sign; reflected minus case is not independent. No convergence or inverse-cascade claim follows from morphology."]}
    (out/"report.json").write_text(json.dumps(public,indent=2,allow_nan=False)+"\n")
    assets={p.name:{"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
    (out/"PASS.json").write_text(json.dumps({"status":"PASS","assets":assets,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+"\n")
    print("ANALYSIS_READY",flush=True)


if __name__ == "__main__":
    assert len(sys.argv) == 3, __doc__
    analyze(Path(sys.argv[1]),Path(sys.argv[2]))
