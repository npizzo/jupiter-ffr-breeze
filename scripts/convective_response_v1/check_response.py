"""Independent finite-window response release checks. Execute on Unity only.

The variable-sign reference N2 does not define a positive b^2/N2 energy.
Kinetic energy and buoyancy variance are therefore checked separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def inner(a, matrix, b):
    return float(np.einsum("aij,ab,bij->", a, matrix, b, optimize=True)/a.shape[-1]**2)


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-30))


def independent_powers(response, state):
    """Physical-space work; independent of response diagnostic implementations."""
    s = response.system
    nh, nw, nv = s.nh, s.nw, s.nvelocity
    q = response.physical(state)
    u, v, w, b = q[:nh], q[nh:2*nh], q[2*nh:nv], q[nv:]
    qx = response.physical(1j*response.kx*state)
    qy = response.physical(1j*response.ky*state)
    ux, vx, wx, bx = qx[:nh], qx[nh:2*nh], qx[2*nh:nv], qx[nv:]
    uy, vy, wy, by = qy[:nh], qy[nh:2*nh], qy[2*nh:nv], qy[nv:]
    rxx = np.einsum("aij,ab,bij->ij", u, s.H, u, optimize=True)
    rxy = np.einsum("aij,ab,bij->ij", u, s.H, v, optimize=True)
    ryy = np.einsum("aij,ab,bij->ij", v, s.H, v, optimize=True)
    stress_u = response.physical(-1j*response.kx*response.spectral(rxx)
                                 -1j*response.ky*response.spectral(rxy))
    stress_v = response.physical(-1j*response.kx*response.spectral(rxy)
                                 -1j*response.ky*response.spectral(ryy))
    bt_divergence = float(np.mean(response.U*stress_u+response.V*stress_v))
    bt_strain = float(np.mean(rxx*response.Ux+rxy*(response.Uy+response.Vx)+ryy*response.Vy))
    viscous = -s.nu*(inner(ux, s.H, ux)+inner(uy, s.H, uy)
                     +inner(vx, s.H, vx)+inner(vy, s.H, vy)
                     +inner(wx, s.W, wx)+inner(wy, s.W, wy)
                     +inner(u, s.Vh, u)+inner(v, s.Vh, v)+inner(w, s.Vw, w))
    return {"bt_power": bt_divergence, "bt_strain_power": bt_strain,
            "buoyancy_power": inner(w, s.C, b), "viscous_power": viscous,
            "thermal_gradient_power": -inner(b, s.B, w),
            "thermal_diffusion_power": -(inner(bx, s.Kx, bx)+inner(by, s.Kx, by)+inner(b, s.Kz, b)),
            "thermal_relaxation_power": -inner(b, s.Gamma, b)}


def check_energy(response, state):
    s = response.system
    q = response.physical(state)
    result = response.rhs(state, with_work=True)
    derivative, reported = result
    dq = response.physical(derivative)
    nv = s.nvelocity
    dkinetic = inner(q[:nv], s.M, dq[:nv])
    dthermal = float(np.mean(np.sum(q[nv:]*dq[nv:], axis=0)))
    p = independent_powers(response, state)
    assert p["viscous_power"] <= 0.
    assert p["thermal_diffusion_power"] <= 0.
    assert p["thermal_relaxation_power"] <= 0.
    kinetic_terms = [-p["bt_power"], p["buoyancy_power"], p["viscous_power"]]
    thermal_terms = [p["thermal_gradient_power"], p["thermal_diffusion_power"], p["thermal_relaxation_power"]]
    kinetic_error = abs(dkinetic-sum(kinetic_terms))/max(sum(map(abs, kinetic_terms)), 1e-30)
    thermal_error = abs(dthermal-sum(thermal_terms))/max(sum(map(abs, thermal_terms)), 1e-30)
    stress_error = abs(p["bt_power"]-p["bt_strain_power"])/max(abs(p["bt_power"]), abs(p["bt_strain_power"]), 1e-20)
    report_error = max(abs(reported[k]-p[k])/max(abs(p[k]), abs(reported[k]), 1e-20)
                       for k in reported if k in p)
    e = response.energy(state)
    kinetic = .5*inner(q[:nv], s.M, q[:nv])
    thermal = .5*float(np.mean(np.sum(q[nv:]**2, axis=0)))
    assert abs(e["kinetic"]-kinetic) < 2e-11*max(kinetic, 1e-20)
    assert abs(e["thermal_variance"]-thermal) < 2e-11*max(thermal, 1e-20)
    checks = {"kinetic_balance_relative": kinetic_error,
              "buoyancy_variance_balance_relative": thermal_error,
              "stress_divergence_strain_relative": stress_error,
              "native_report_independent_power_relative": report_error}
    assert max(checks.values()) < 2e-9, checks
    return checks


def reflection(response, state):
    physical = response.physical(state)
    reflected = physical[:, (-np.arange(response.n)) % response.n, :].copy()
    reflected[:response.system.nh] *= -1
    return response.spectral(reflected)


def check_projection(response, seed=199):
    rng = np.random.default_rng(seed)
    s = response.system
    a = response.spectral(rng.normal(size=(s.nvelocity+s.nb+1, response.n, response.n)))
    y = response.project(a)
    d = (1j*response.kx*y[:s.nh]+1j*response.ky*y[s.nh:2*s.nh]
         +np.einsum("ab,bij->aij", s.J, y[2*s.nh:s.nvelocity], optimize=True))
    continuity = float(np.linalg.norm(d)/max(np.linalg.norm(y[:s.nvelocity])*np.sqrt(np.max(response.k2)), 1e-30))
    idempotence = relative(response.project(y), y)
    hermitian = relative(response.spectral(response.physical(y)), y)
    # Metric orthogonality includes removed Fourier modes as well as pressure.
    q, removed = response.physical(y), response.physical(a-y)
    orthogonality = abs(inner(q[:s.nvelocity], s.M, removed[:s.nvelocity]))
    orthogonality /= max(np.sqrt(inner(q[:s.nvelocity], s.M, q[:s.nvelocity])
                                 *inner(removed[:s.nvelocity], s.M, removed[:s.nvelocity])), 1e-30)
    checks = {"continuity_relative": continuity, "projection_idempotence": idempotence,
              "real_field_hermitian_relative": hermitian, "projection_work_orthogonality": orthogonality}
    assert max(checks.values()) < 2e-10, checks
    return y, checks


def centrifugal_margin(response):
    # Actual filtered periodic velocity; the axisymmetric criterion is screened
    # only in the prescribed central vortex, not at the periodic box corners.
    r2 = response.X**2+response.Ygrid**2
    twice_omega = np.divide(2*(-response.Ygrid*response.U+response.X*response.V),
                            r2, out=response.Z.copy(), where=r2 > 1e-28)
    core = r2 <= (2*response.radius)**2
    phi = (response.f+twice_omega)*(response.f+response.Z)/response.f**2
    margin = float(np.min(phi[core]))
    assert margin > .4, margin
    return margin


def same_physical_initial_state(response1, response2, seed):
    y1, y2 = response1.seed(seed), response2.seed(seed)
    blocks = []
    for response, y in ((response1, y1), (response2, y2)):
        s = response.system
        blocks.append(np.concatenate((
            np.einsum("ab,bij->aij", s.h_coeff, y[:s.nh], optimize=True),
            np.einsum("ab,bij->aij", s.h_coeff, y[s.nh:2*s.nh], optimize=True),
            np.einsum("ab,bij->aij", s.g_coeff, y[2*s.nh:s.nvelocity], optimize=True),
            y[s.nvelocity:]), axis=0))
    error = relative(blocks[0], blocks[1])
    assert error < 2e-10, error
    return error


def advance(response, initial, dt, steps):
    y = initial.copy()
    integrated = {}
    for _ in range(steps):
        y, work = response.rk4(y, dt)
        for key, value in work.items():
            integrated[key] = integrated.get(key, 0.)+value
    return y, integrated


def timestep_checks(response, dt):
    y = response.seed(31415)
    coarse, wc = advance(response, y, dt, 10)
    fine, wf = advance(response, y, dt/2, 20)
    finer, _ = advance(response, y, dt/4, 40)
    difference = relative(coarse, fine)
    fine_difference = relative(fine, finer)
    assert difference < 1e-3, difference
    assert fine_difference < max(.15*difference, 1e-11), (difference, fine_difference)
    def residual(last, work):
        delta = response.energy(last)["kinetic"]-response.energy(y)["kinetic"]
        expected = -work["bt_power"]+work["buoyancy_power"]+work["viscous_power"]
        return float(delta-expected)
    rc, rf = residual(coarse, wc), residual(fine, wf)
    # RK quadrature of instantaneous powers is not an exact quadratic invariant.
    assert abs(rf) < max(.2*abs(rc), 1e-12*response.energy(y)["kinetic"]), (rc, rf)
    return {"coarse_dt": dt, "coarse_half_dt_relative_state_difference": difference,
            "half_quarter_dt_relative_state_difference": fine_difference,
            "coarse_dt_kinetic_work_residual": rc, "half_dt_kinetic_work_residual": rf}


def hashes():
    names = ("model.py", "optimized.py", "response.py", "run.py", "check_response.py")
    return {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in names}


def check_real_complex_product(response):
    rng = np.random.default_rng(22198)
    matrix = rng.normal(size=(7, 5))
    state = rng.normal(size=(5, 13, 9))+1j*rng.normal(size=(5, 13, 9))
    errors = []
    for field in (state, state[:, ::2, :], state[..., ::2]):
        actual = response.multiply(matrix, field)
        expected = (matrix@field.reshape(5, -1)).reshape((7,)+field.shape[1:])
        errors.append(relative(actual, expected))
    error = max(errors)
    assert error < 2e-14, error
    return error


def count_event_aligned_steps(dt, target=10., scalar_interval=.2, frame_interval=.5):
    t, next_scalar, next_frame, steps = 0., scalar_interval, frame_interval, 0
    while t < target-1e-10:
        step = min(dt, target-t, next_scalar-t, next_frame-t)
        assert step > 1e-12
        t += step
        steps += 1
        if abs(t-next_scalar) < 1e-9:
            next_scalar = round(next_scalar+scalar_interval, 10)
        if abs(t-next_frame) < 1e-9:
            next_frame = round(next_frame+frame_interval, 10)
    return steps


def main(basis, out):
    from optimized import AugmentedSystem
    from response import Response, load_basis

    started = time.monotonic()
    basis, out = Path(basis).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    systems, metadata = load_basis(basis)
    length = float(metadata["length"])
    radius = float(metadata["radius"])
    kh_seed = float(metadata["kh_seed"])
    responses = {m: Response(s, n=128, length=length, radius=radius,
                              vortex_sign=1, f=1, core_rossby=.25, kh_seed=kh_seed)
                 for m, s in systems.items()}
    checks = {}
    checks["real_matrix_complex_field_relative"] = check_real_complex_product(responses[1])
    stable_steps = []
    for m, response in responses.items():
        random, projection = check_projection(response)
        checks[f"modes_{m}_projection"] = projection
        checks[f"modes_{m}_energy"] = check_energy(response, random)
        checks[f"modes_{m}_cyclone_centrifugal_margin"] = centrifugal_margin(response)
        anti = Response(response.system, n=128, length=length, radius=radius,
                        vortex_sign=-1, f=1, core_rossby=.25, kh_seed=kh_seed)
        checks[f"modes_{m}_anticyclone_centrifugal_margin"] = centrifugal_margin(anti)
        control = AugmentedSystem(response.system.g_coeff, response.system.h_coeff,
                                  response.system.nb, response.system.n2_si, control=True)
        no_gradient = Response(control, n=128, length=length, radius=radius,
                               vortex_sign=1, f=1, core_rossby=.25, kh_seed=kh_seed)
        stable_steps.extend((response.stable_step(), no_gradient.stable_step()))
        for seed in (31415, 27182):
            y = response.seed(seed)
            assert np.array_equal(y, anti.seed(seed)), "Vortex sign changed disturbance seed"
            assert np.array_equal(y, no_gradient.seed(seed)), "Gradient removal changed disturbance seed"
        linear_error = relative(response.rhs(2.3*random), 2.3*response.rhs(random))
        assert linear_error < 2e-10, linear_error
        checks[f"modes_{m}_linearity_relative"] = linear_error
        _, powers = response.rhs(random, with_work=True)
        _, doubled_powers = response.rhs(2*random, with_work=True)
        power_scale_error = max(abs(doubled_powers[k]-4*powers[k])
                                /max(abs(doubled_powers[k]), abs(4*powers[k]), 1e-20)
                                for k in powers)
        assert power_scale_error < 2e-10, power_scale_error
        checks[f"modes_{m}_quadratic_power_scaling_relative"] = power_scale_error
        reflected = Response(response.system, n=128, length=length, radius=radius,
                             vortex_sign=-1, f=-1, core_rossby=.25, kh_seed=kh_seed)
        mirror_error = relative(reflection(response, response.rhs(random)),
                                reflected.rhs(reflection(response, random)))
        assert mirror_error < 2e-10, mirror_error
        checks[f"modes_{m}_reflected_f_rhs_relative"] = mirror_error
        _, reflected_powers = reflected.rhs(reflection(response, random), with_work=True)
        mirror_power_error = max(abs(reflected_powers[k]-powers[k])
                                 /max(abs(reflected_powers[k]), abs(powers[k]), 1e-20)
                                 for k in powers)
        assert mirror_power_error < 2e-9, mirror_power_error
        checks[f"modes_{m}_reflected_f_power_relative"] = mirror_power_error
        print("ALGEBRA_PASS", m, flush=True)
    checks["nested_initial_physical_relative"] = {
        str(seed): same_physical_initial_state(responses[1], responses[2], seed)
        for seed in (31415, 27182)}
    response = responses[2]
    production_dt = min(stable_steps)
    checks["timestep"] = timestep_checks(response, production_dt)

    # Exact actual-driver invocation, with the configured cadence and final frame.
    pilot = out/"driver-pilot"
    pilot_start = time.monotonic()
    command = [sys.executable, str(Path(__file__).with_name("run.py")), str(basis),
               str(pilot), "--pilot", "--dt", str(production_dt)]
    with (out/"driver-pilot.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    pilot_wall = time.monotonic()-pilot_start
    import csv
    scalar_rows = list(csv.DictReader((pilot/"scalars.csv").open()))
    frame_rows = list(csv.DictReader((pilot/"frames.csv").open()))
    scalar_times = [float(row["t_fabs"]) for row in scalar_rows]
    frame_times = [float(row["t_fabs"]) for row in frame_rows]
    assert np.allclose(scalar_times, [0., .2, .4], rtol=0, atol=2e-12), scalar_times
    assert np.allclose(frame_times, [0., .4], rtol=0, atol=2e-12), frame_times
    status = json.loads((pilot/"status.json").read_text())
    assert status["status"] == "COMPLETED", status
    assert abs(status["t_fabs"]-.4) < 2e-12, status
    assert (pilot/"checkpoint.npz").stat().st_size > 0
    checks["actual_driver_pilot"] = {"status": "PASS", "target": .4,
                                     "scalar_times": scalar_times, "frame_times": frame_times,
                                     "elapsed_seconds": pilot_wall}
    # Warmed full steps, including stage work; benchmark the larger retained space.
    y = response.seed(31415)
    y, _ = response.rk4(y, production_dt)
    benchmark_start = time.monotonic()
    for _ in range(20):
        y, _ = response.rk4(y, production_dt)
    seconds_per_step = (time.monotonic()-benchmark_start)/20
    diagnostic_start = time.monotonic()
    response.diagnostics(y)
    seconds_per_diagnostic = time.monotonic()-diagnostic_start
    # Account exactly for .2 scalar and .5 frame event-aligned timesteps,
    # real scalar/frame/checkpoint output and repeated setup. The larger space
    # sets the timestep benchmark; scaling output by its coefficient count is
    # conservative for the fixed three-height fields.
    steps_per_case = count_event_aligned_steps(production_dt)
    assert status["steps"] == count_event_aligned_steps(production_dt, target=.4)
    timing = status["timings"]
    scalar_seconds = max(timing["scalar_output_wall_seconds"]/status["scalar_samples"],
                         seconds_per_diagnostic)
    size_ratio = responses[2].nvar/responses[1].nvar
    frame_seconds = timing["frame_output_wall_seconds"]/status["frames"]*size_ratio
    checkpoint_seconds = timing["checkpoint_wall_seconds"]*size_ratio
    startup_seconds = max(pilot_wall-status["elapsed_seconds"], 0.)
    setup_seconds = timing["setup_wall_seconds"]
    estimated_total = 18*(steps_per_case*seconds_per_step+51*scalar_seconds
                          +21*frame_seconds+checkpoint_seconds+setup_seconds+startup_seconds)
    checks["throughput"] = {"seconds_per_full_rk4_step": seconds_per_step,
                            "seconds_per_diagnostic": seconds_per_diagnostic,
                            "actual_driver_pilot_seconds": pilot_wall,
                            "production_dt": production_dt,
                            "steps_per_case_including_event_alignment": steps_per_case,
                            "frame_output_seconds_per_frame": frame_seconds,
                            "checkpoint_seconds_per_case": checkpoint_seconds,
                            "projected_18_case_cpu_seconds": estimated_total,
                            "release_limit_cpu_seconds": 1800.}
    cost_pass = estimated_total <= 1800.
    report = {"status": "PASS" if cost_pass else "COST_LIMIT_FAIL", "compute_location": "Unity", "checks": checks,
              "numerical_checks_status": "PASS", "cost_gate_pass": cost_pass,
              "basis_sha256": hashlib.sha256(basis.read_bytes()).hexdigest(),
              "code_sha256": hashes(), "elapsed_seconds": time.monotonic()-started,
              "budget_semantics": "Separate kinetic-energy and buoyancy-variance balances; RK-stage work integrals have a measured finite-step residual.",
              "maximum_campaign_time": 10., "maximum_case_count": 18}
    (out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "projected_18_case_cpu_seconds": estimated_total,
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    assert cost_pass, checks["throughput"]
    (out/"PASS.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("basis")
    parser.add_argument("out")
    args = parser.parse_args()
    main(args.basis, args.out)
