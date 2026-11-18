"""Bounded 18-case linear response campaign; all execution on Unity."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from model import F_ABS
from response import Response, load_basis, with_control

SEEDS = (20261118, 20261119)


def code_hashes():
    root = Path(__file__).parent
    names = ("model.py", "optimized.py", "response.py", "run.py", "check_response.py")
    return {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in names}


def run_case(response, out, seed, target, dt, reflected=False, metadata=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    state = response.seed(seed)
    if reflected:
        state = response.reflect_state(state)
    initial = response.energy(state)
    integrals = {k: 0. for k in ("bt_power", "buoyancy_power", "viscous_power", "thermal_gradient_power", "thermal_diffusion_power", "thermal_relaxation_power")}
    config = {"model": "Frozen-profile Boussinesq linear response to a prescribed vortex",
              "compute_location": "Unity", "n": response.n, "length_over_H": response.length,
              "radius_over_H": response.radius, "f_over_fabs": response.f,
              "vortex_sign": response.vortex_sign, "core_abs_zeta_over_fabs": .25,
              "convective_gradient_removed": response.system.control,
              "seed": seed, "reflected_seed": reflected, "target_time_fabs": target,
              "dt_ceiling_fabs": dt, "initial_perturbation_rms_over_base": 1e-6,
              "linear_amplitude_guard": .01, "vertical_velocity_functions": response.nw,
              "horizontal_velocity_functions": response.nh, "thermal_functions": response.system.nb+1,
              "frame_interval_fabs": .5, "scalar_interval_fabs": .2,
              "heights_m": [20000, 60000, 100000], "code_sha256": code_hashes(),
              "interpretation": "Virtual stress feedback; background vortex and thermal profile are held fixed.",
              "tracer": "No passive tracer is evolved; saved scalar is buoyancy perturbation.",
              "budget": "RK4-stage integrated powers with independently measured quadratic-energy residual.",
              **(metadata or {})}
    (out/"config.json").write_text(json.dumps(config, indent=2)+"\n")
    np.savez_compressed(out/"background.npz", U=response.U, V=response.V, zeta=response.Z,
                        x=response.x, y=response.y)
    (out/"frames.csv").write_text("t_fabs,file\n")
    scalar_path = out/"scalars.csv"
    t = 0.
    next_scalar, next_frame = .2, .5
    frame_count = 0
    scalar_count = 0
    last_scalar = last_frame = -1.
    timings = {"setup_wall_seconds": time.monotonic()-started,
               "scalar_output_wall_seconds": 0., "frame_output_wall_seconds": 0.,
               "checkpoint_wall_seconds": 0.}

    def save_scalar():
        nonlocal scalar_count, last_scalar
        if abs(t-last_scalar) < 1e-10:
            return
        sample_started = time.monotonic()
        d = response.diagnostics(state)
        assert np.isfinite(list(d.values())).all()
        assert d["perturbation_rms_over_base"] < .01, "Linear amplitude guard reached"
        dk = d["kinetic"]-initial["kinetic"]
        db = d["thermal_variance"]-initial["thermal_variance"]
        ik = -integrals["bt_power"]+integrals["buoyancy_power"]+integrals["viscous_power"]
        ib = sum(integrals[k] for k in ("thermal_gradient_power", "thermal_diffusion_power", "thermal_relaxation_power"))
        row = {"t_fabs": t, "time_hours": t/F_ABS/3600, **d,
               **{"integral_"+k: v for k, v in integrals.items()},
               "kinetic_change": dk, "thermal_variance_change": db,
               "integrated_kinetic_residual": dk-ik,
               "integrated_thermal_variance_residual": db-ib,
               "bt_power_over_perturbation_energy": d["bt_power"]/max(d["kinetic"], np.finfo(float).tiny)}
        with scalar_path.open("a") as file:
            writer = csv.DictWriter(file, fieldnames=list(row))
            if scalar_count == 0:
                writer.writeheader()
            writer.writerow(row)
        scalar_count += 1
        last_scalar = t
        timings["scalar_output_wall_seconds"] += time.monotonic()-sample_started

    def save_frame():
        nonlocal frame_count, last_frame
        if abs(t-last_frame) < 1e-10:
            return
        frame_started = time.monotonic()
        name = f"planes_{frame_count:04d}.npz"
        np.savez_compressed(out/name, t_fabs=t, **response.plane_fields(state))
        with (out/"frames.csv").open("a") as file:
            file.write(f"{t:.12g},{name}\n")
        frame_count += 1
        last_frame = t
        timings["frame_output_wall_seconds"] += time.monotonic()-frame_started

    save_scalar()
    save_frame()
    steps = 0
    while t < target-1e-10:
        step = min(dt, target-t, next_scalar-t, next_frame-t)
        assert step > 1e-12
        state, work = response.rk4(state, step)
        for key in integrals:
            integrals[key] += work[key]
        t += step
        steps += 1
        if abs(t-next_scalar) < 1e-9:
            save_scalar()
            next_scalar = round(next_scalar+.2, 10)
        if abs(t-next_frame) < 1e-9:
            save_frame()
            next_frame = round(next_frame+.5, 10)
        if time.monotonic()-started > 2400:
            raise RuntimeError("Per-case emergency wall limit exceeded")
    save_scalar()
    save_frame()
    checkpoint_started = time.monotonic()
    np.savez_compressed(out/"checkpoint.npz", state=state, t_fabs=t)
    timings["checkpoint_wall_seconds"] = time.monotonic()-checkpoint_started
    status = {"status": "COMPLETED", "t_fabs": t, "steps": steps,
              "frames": frame_count, "scalar_samples": scalar_count,
              "elapsed_seconds": time.monotonic()-started, "timings": timings}
    (out/"status.json").write_text(json.dumps(status, indent=2)+"\n")
    print(json.dumps({"case": out.name, **status}), flush=True)
    return status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("basis")
    ap.add_argument("out")
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--t-stop", type=float)
    ap.add_argument("--dt", type=float)
    ap.add_argument("--validation")
    args = ap.parse_args()
    systems, geometry = load_basis(args.basis)
    target = args.t_stop if args.t_stop is not None else (.4 if args.pilot else 10.)
    assert 0 < target <= 10.
    if args.pilot:
        response = Response(systems[1], **geometry)
        dt = args.dt or response.stable_step()
        run_case(response, args.out, SEEDS[0], target, dt,
                 metadata={"selected_complete_vertical_modes": 1, "pilot": True})
        return
    if not args.validation:
        raise RuntimeError("Production requires independent response validation")
    gate = json.loads(Path(args.validation).read_text())
    assert gate["status"] == "PASS"
    assert gate["code_sha256"] == code_hashes()
    assert gate["basis_sha256"] == hashlib.sha256(Path(args.basis).read_bytes()).hexdigest()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    runs = []
    for nm in (1, 2):
        for control in (False, True):
            system = with_control(systems[nm]) if control else systems[nm]
            for seed in SEEDS:
                for sign in (1, -1):
                    runs.append((nm, control, seed, sign, 1., False, system))
    for nm in (1, 2):
        runs.append((nm, False, SEEDS[0], -1, -1., True, systems[nm]))
    tested_dt = float(gate["checks"]["throughput"]["production_dt"])
    dt = args.dt or tested_dt
    assert abs(dt-tested_dt) < 1e-12, "Production must use the accuracy- and cost-tested timestep"
    assert dt <= .1
    record = {"status": "RUNNING", "cases": [], "dt_ceiling": dt, "target": target,
              "basis_sha256": gate["basis_sha256"], "code_sha256": code_hashes()}
    (out/"campaign.json").write_text(json.dumps(record, indent=2)+"\n")
    for nm, control, seed, sign, f, reflected, system in runs:
        name = f"m{nm}_{'neutral' if control else 'heated'}_seed{seed}_{'plus' if sign>0 else 'minus'}_f{'plus' if f>0 else 'minus'}"
        response = Response(system, f=f, vortex_sign=sign, **geometry)
        status = run_case(response, out/name, seed, target, dt, reflected=reflected,
                          metadata={"selected_complete_vertical_modes": nm, "pilot": False})
        record["cases"].append({"name": name, **status})
        (out/"campaign.json").write_text(json.dumps(record, indent=2)+"\n")
        if time.monotonic()-started > 2400:
            raise RuntimeError("Campaign internal wall limit reached; completed cases retained")
    record.update(status="COMPLETED", elapsed_seconds=time.monotonic()-started)
    (out/"campaign.json").write_text(json.dumps(record, indent=2)+"\n")


if __name__ == "__main__":
    main()
