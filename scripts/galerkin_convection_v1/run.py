"""Bounded reference-layer Galerkin campaign. Execute scientific work on Unity."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from model import GalerkinModel


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def prepare(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    modes = []
    for nx in range(-4, 5):
        for ny in range(-4, 5):
            if nx == ny == 0:
                continue
            for nz in (1, 2):
                kh2 = 20.**2*(nx*nx+ny*ny)
                kz2 = (np.pi*nz)**2
                critical = (kz2+1e-10*(kh2+kz2)**3)/kh2
                modes.append((float(critical), nx, ny, nz))
    critical, nx, ny, nz = min(modes)
    config = {
        "campaign": "galerkin-convection-v1", "kfund": 20.,
        "nu": 1e-5, "kappa": 1e-5, "G": 2*critical, "f": 1.,
        "dt": .02, "duration": 200., "diagnostic_interval": .5,
        "snapshot_interval": 1., "heights": [.2, .5, .8], "nxy": 64,
        "velocity_rms": 1e-4, "temperature_rms": 1e-5,
        "discrete_onset_G": critical, "one_critical_mode": [nx, ny, nz],
        "cpu_budget_seconds": 1800,
        "cases": [
            {"name": "plus-seed18", "K": 4, "Nv": 2, "seed": 20261118, "f_sign": 1, "gradient": True},
            {"name": "minus-reflected18", "K": 4, "Nv": 2, "seed": 20261118, "f_sign": -1, "gradient": True},
            {"name": "plus-seed19", "K": 4, "Nv": 2, "seed": 20261119, "f_sign": 1, "gradient": True},
            {"name": "zero-gradient18", "K": 4, "Nv": 2, "seed": 20261118, "f_sign": 1, "gradient": False},
            {"name": "refined-plus18", "K": 6, "Nv": 3, "seed": 20261118, "f_sign": 1, "gradient": True},
        ],
        "meaning": "Finite-mode nonlinear rotating-convection reference, not the capped Breeze atmosphere or a resolved Jupiter simulation.",
        "control_caveat": "Zero-gradient control retains matched initial buoyancy disturbances, which remain dynamically active.",
        "initial_condition": "Broad real divergence-free convective disturbance, zero barotropic and horizontally averaged modes; refined run embeds the identical base coefficients.",
    }
    write_json(output/"config.json", config)
    print(json.dumps(config), flush=True)


def make_model(config, case):
    return GalerkinModel(K=case["K"], Nv=case["Nv"], kfund=config["kfund"],
                         nu=config["nu"], kappa=config["kappa"],
                         G=config["G"] if case["gradient"] else 0.,
                         f=config["f"]*case["f_sign"])


def embed(source, old, new):
    result = np.zeros((4,)+new.shape, complex)
    for i, nx in enumerate(old.ix):
        for j, ny in enumerate(old.iy):
            for k, nz in enumerate(old.iz):
                if old.mask_u[i,j,k] or old.mask_theta[i,j,k]:
                    result[:, nx % new.shape[0], ny % new.shape[1], nz % new.shape[2]] = source[:,i,j,k]
    return new.project(result)


def run(config_path, case_name, root):
    config = json.loads(Path(config_path).read_text())
    case = next(c for c in config["cases"] if c["name"] == case_name)
    model = make_model(config, case)
    base_case = dict(case, K=4, Nv=2, f_sign=1, gradient=True)
    base = make_model(config, base_case)
    state = base.initial_state(seed=case["seed"], velocity_rms=config["velocity_rms"],
                               temperature_rms=config["temperature_rms"], barotropic=False)
    state = embed(state, base, model)
    if case["f_sign"] < 0:
        state = model.reflect(state)
    output = Path(root)/case_name
    output.mkdir(parents=True, exist_ok=False)
    source_dir = Path(__file__).parent
    metadata = {"config": config, "case": case, "model": model.metadata(),
                "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in source_dir.glob("*.py")}}
    write_json(output/"metadata.json", metadata)
    np.savez_compressed(output/"initial.npz", state=state)
    initial = model.diagnostics(state)
    accum, absolute_work = {}, {}
    diagnostic_rows, plane_rows, times = [], [], []
    dt = config["dt"]
    steps = round(config["duration"]/dt)
    diagnostic_stride = round(config["diagnostic_interval"]/dt)
    snapshot_stride = round(config["snapshot_interval"]/dt)
    started = time.monotonic()
    def preserve_failure(reason):
        np.savez_compressed(output/"failure-state.npz", state=state, time=t)
        write_json(output/"diagnostics.json", diagnostic_rows)
        if plane_rows:
            np.savez_compressed(output/"partial-planes.npz", times=times,
                                heights=config["heights"], x=plane_rows[0]["x"], y=plane_rows[0]["y"],
                                **{key: np.stack([r[key] for r in plane_rows]) for key in
                                   ("vorticity", "buoyancy_perturbation", "barotropic_vorticity")})
        write_json(output/"FAILED.json", {"reason": reason, "time": t})
        raise RuntimeError(reason)
    for step in range(steps+1):
        t = step*dt
        if step % diagnostic_stride == 0 or step == steps:
            row = model.diagnostics(state)
            row.update(time=t, step=step, integrated_work=dict(accum),
                       accumulated_absolute_work=dict(absolute_work))
            for quantity in ("kinetic", "thermal_variance", "potential_moment", "barotropic_kinetic"):
                defect = row[quantity]-initial[quantity]-accum.get(quantity, 0.)
                scale = max(abs(initial[quantity]), abs(row[quantity]), absolute_work.get(quantity, 0.), 1e-25)
                row[quantity+"_integrated_budget_relative_defect"] = defect/scale
            # A conservative per-coordinate advective CFL on the quadrature grid.
            fields = model.physical(state)
            rates = [np.max(abs(fields[d]))*model.shape[d]/([2*np.pi/config["kfund"]]*2+[2.])[d]
                     for d in range(3)]
            row["advective_cfl_sum"] = float(dt*sum(rates))
            if not all(np.isfinite(v) for v in row.values() if isinstance(v, (float, int))):
                preserve_failure("Nonfinite diagnostic; production stopped")
            if row["advective_cfl_sum"] > .6:
                preserve_failure("Declared CFL limit exceeded; production stopped")
            if any(abs(row[q+"_integrated_budget_relative_defect"]) > 5e-5 for q in
                   ("kinetic", "thermal_variance", "potential_moment", "barotropic_kinetic")):
                preserve_failure("Integrated budget defect exceeded 5e-5; production stopped")
            diagnostic_rows.append(row)
            if step % (10*diagnostic_stride) == 0:
                print(json.dumps({"case": case_name, "t": t, "K": row["kinetic"],
                                  "Kbt": row["barotropic_kinetic"],
                                  "skew_bt": row["barotropic_vorticity_skewness"],
                                  "elapsed_seconds": time.monotonic()-started}), flush=True)
                write_json(output/"progress.json", diagnostic_rows[-1])
        if step % snapshot_stride == 0 or step == steps:
            times.append(t)
            plane_rows.append(model.planes(state, heights=config["heights"], nxy=config["nxy"]))
        if step == steps:
            break
        state, work = model.rk4_step_with_work(state, dt)
        if not np.isfinite(state).all():
            preserve_failure("Nonfinite state after RK step; production stopped")
        for key, value in work.items():
            accum[key] = accum.get(key, 0.)+value
            absolute_work[key] = absolute_work.get(key, 0.)+abs(value)
        if time.monotonic()-started > config["cpu_budget_seconds"]:
            preserve_failure("Campaign wall-time guard exceeded")
    np.savez_compressed(output/"final.npz", state=state)
    planes = {key: np.stack([r[key] for r in plane_rows])
              for key in ("vorticity", "buoyancy_perturbation", "barotropic_vorticity")}
    np.savez_compressed(output/"planes.npz", times=times, heights=config["heights"],
                        x=plane_rows[0]["x"], y=plane_rows[0]["y"], **planes)
    write_json(output/"diagnostics.json", diagnostic_rows)
    write_json(output/"COMPLETE.json", {"initial": initial, "final": diagnostic_rows[-1],
                "elapsed_seconds": time.monotonic()-started, "saved_frames": len(times)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("output")
    p = sub.add_parser("run")
    p.add_argument("config")
    p.add_argument("case")
    p.add_argument("root")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.output)
    else:
        run(args.config, args.case, args.root)
