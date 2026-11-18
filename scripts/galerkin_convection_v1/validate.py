#!/usr/bin/env python3
"""Independent release checks for the nonlinear convection Galerkin model.

Run on Unity only.  This program checks physical identities and analytic mode
solutions, rather than trusting the engine's diagnostic implementations.  A
PASS.json is written only after the algebra, time refinement, and bounded cost
checks pass.  The result is a numerical-implementation gate, not evidence of a
large-scale vortex or of adequacy for the capped Breeze atmosphere.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np

from model import GalerkinModel


def physical(coefficients):
    return np.fft.ifftn(coefficients, axes=(-3, -2, -1), norm="forward").real


def spectral(fields):
    return np.fft.fftn(fields, axes=(-3, -2, -1), norm="forward")


def norm(array):
    return float(np.linalg.norm(np.asarray(array).ravel()))


def relative(actual, expected, floor=1e-30):
    return norm(np.asarray(actual)-np.asarray(expected))/max(norm(expected), floor)


def inner(left, right):
    """Volume mean, summing vector components if supplied."""
    return float(np.mean(np.sum(left*right, axis=0))) if left.ndim == 4 else float(np.mean(left*right))


def potential_moment(model, theta):
    """Independent exact half-layer integral of (1-z)*theta.

    Full complex exponential integral; unlike the engine's sine-specific
    implementation this also defines the n=0 value and even/odd endpoint term.
    """
    theta = theta[3] if theta.ndim == 4 else theta
    numbers = model.iz.astype(float)
    weights = np.zeros(numbers.shape, dtype=complex)
    nonzero = numbers != 0
    q = np.pi*numbers[nonzero]
    weights[nonzero] = 1j/q + (1-np.cos(q))/(q*q)
    weights[~nonzero] = .5
    return float(np.real(np.sum(theta[0, 0]*weights)))


def half_average(model, coefficients):
    """Exact physical-domain average; handles odd z parity as well as even."""
    numbers = model.iz.astype(float)
    weights = np.ones(numbers.shape, dtype=complex)
    nonzero = numbers != 0
    q = np.pi*numbers[nonzero]
    weights[nonzero] = 1j*(1-np.cos(q))/q
    return np.real(np.sum(coefficients[..., 0, 0, :]*weights, axis=-1))


def curl(model, velocity):
    kx, ky, kz = model.k
    u, v, w = velocity
    return 1j*np.stack((ky*w-kz*v, kz*u-kx*w, kx*v-ky*u))


def transfer(source, target, state):
    """Embed identical normalized Fourier amplitudes; never reseed refinement."""
    out = np.zeros((4,)+target.shape, dtype=complex)
    keep = np.any(abs(state) > 0, axis=0)
    for i, j, n in zip(*np.nonzero(keep)):
        ix, iy, iz = source.ix[i], source.iy[j], source.iz[n]
        if not (abs(ix) < target.shape[0]/2 and abs(iy) < target.shape[1]/2
                and abs(iz) < target.shape[2]/2):
            raise ValueError("Target grid cannot represent a source coefficient")
        out[:, ix % target.shape[0], iy % target.shape[1], iz % target.shape[2]] = state[:, i, j, n]
    return target.project(out)


def random_state(model, seed=92751, velocity_rms=.01, temperature_rms=.001):
    """Broad independent test field including mean shear and barotropic modes."""
    rng = np.random.default_rng(seed)
    state = model.project(spectral(rng.normal(size=(4,)+model.shape)))
    field = physical(state)
    state[:3] *= velocity_rms/math.sqrt(inner(field[:3], field[:3]))
    state[3] *= temperature_rms/math.sqrt(inner(field[3], field[3]))
    return model.project(state)


def copy_model(model, **changes):
    values = {key: getattr(model, key) for key in ("K", "Nv", "kfund", "nu", "kappa", "G", "f")}
    values.update(changes)
    return GalerkinModel(**values)


class Checks:
    def __init__(self):
        self.results = []

    def add(self, name, value, tolerance, **detail):
        value = float(value)
        self.results.append({"name": name, "error": value if np.isfinite(value) else None,
                             "tolerance": float(tolerance),
                             "passed": bool(np.isfinite(value) and value <= tolerance),
                             **detail})

    @property
    def passed(self):
        return all(item["passed"] for item in self.results)


def algebra_checks(model, checks, prefix):
    state = random_state(model)
    fields = physical(state)
    derivatives = physical(1j*model.k[:, None]*state[None])
    terms = model.tendencies(state)
    tendencies = {key: physical(value) for key, value in terms.items()}
    kscale = max(norm(model.k*state[:3]), 1e-30)
    checks.add(prefix+"/divergence", norm(np.sum(model.k*state[:3], axis=0))/kscale, 2e-13)
    checks.add(prefix+"/projection_idempotent", relative(model.project(state), state), 2e-13)
    complex_fields = np.fft.ifftn(state, axes=(-3, -2, -1), norm="forward")
    checks.add(prefix+"/reality", norm(complex_fields.imag)/norm(fields), 2e-13)
    # z=0 and z=1 are exact collocation locations on the even-sized grid.
    walls = (0, model.shape[2]//2)
    checks.add(prefix+"/impermeable_fixed_temperature", norm(fields[2:, :, :, walls])/norm(fields), 2e-13)
    checks.add(prefix+"/stress_free", norm(derivatives[2, :2, :, :, walls])/norm(derivatives), 2e-13)
    checks.add(prefix+"/rhs_term_sum", relative(model.rhs(state), sum(terms.values())), 3e-13)

    advection = tendencies["advective"]
    ke_scale = max(norm(state[:3])*norm(terms["advective"][:3]), 1e-30)
    q_scale = max(norm(state[3])*norm(terms["advective"][3]), 1e-30)
    checks.add(prefix+"/nonlinear_kinetic_conservation", abs(inner(fields[:3], advection[:3]))/ke_scale, 3e-13)
    checks.add(prefix+"/nonlinear_thermal_variance_conservation", abs(inner(fields[3], advection[3]))/q_scale, 3e-13)
    rotation_scale = max(norm(state[:3])*norm(terms["coriolis"][:3]), 1e-30)
    checks.add(prefix+"/coriolis_no_work", abs(inner(fields[:3], tendencies["coriolis"][:3]))/rotation_scale, 3e-13)
    heat = inner(fields[2], fields[3])
    heat_scale = max(norm(state[2])*norm(state[3]), 1e-30)
    checks.add(prefix+"/buoyancy_conversion", abs(inner(fields[:3], tendencies["buoyancy"][:3])-heat)/heat_scale, 3e-13)
    expected_viscous = -model.nu*float(np.sum(abs(1j*model.k[:, None]*state[None, :3])**2))
    expected_thermal = -model.kappa*float(np.sum(abs(1j*model.k*state[3])**2))
    checks.add(prefix+"/viscous_dissipation", relative(inner(fields[:3], tendencies["diffusion"][:3]), expected_viscous), 3e-13)
    checks.add(prefix+"/thermal_dissipation", relative(inner(fields[3], tendencies["diffusion"][3]), expected_thermal), 3e-13)
    checks.add(prefix+"/background_thermal_power", abs(inner(fields[3], tendencies["background_gradient"][3])-model.G*heat)/max(model.G*heat_scale, 1e-30), 3e-13)

    raw_advection = -np.sum(fields[:3, None]*derivatives, axis=0)
    raw_hat = spectral(raw_advection)
    raw_masked = raw_hat[:3]*model.mask_u
    pressure = raw_masked-terms["advective"][:3]
    checks.add(prefix+"/pressure_no_work", abs(inner(fields[:3], physical(pressure)))/ke_scale, 4e-13)
    checks.add(prefix+"/pressure_curl_free", norm(curl(model, pressure))/max(norm(model.k*raw_masked), 1e-30), 4e-13)
    checks.add(prefix+"/nonlinear_potential_exchange", abs(potential_moment(model, terms["advective"])+heat)/heat_scale, 5e-13)
    checks.add(prefix+"/engine_potential_moment", abs(model.potential_moment(state)-potential_moment(model, state))/max(norm(state[3]), 1e-30), 3e-13)

    # Vertically averaged Reynolds-stress work is independently reconstructed.
    bt = np.zeros_like(state[:3])
    bt[:2, :, :, 0] = state[:2, :, :, 0]
    U = physical(bt)
    fluctuation = fields[:3]-U
    stress = np.mean(fluctuation[:2, None]*fluctuation[None, :2], axis=-1)
    grad_U = physical(1j*model.k[:2, None]*bt[None, :2])
    stress_work = float(np.mean(np.sum(stress*grad_U.transpose(1, 0, 2, 3, 4).mean(axis=-1), axis=(0, 1))))
    mean_work = inner(U, advection[:3])
    fluctuation_work = inner(fluctuation, advection[:3])
    checks.add(prefix+"/barotropic_stress_exchange", abs(mean_work-stress_work)/ke_scale, 4e-13)
    checks.add(prefix+"/barotropic_fluctuation_exchange", abs(mean_work+fluctuation_work)/ke_scale, 4e-13)

    # Whole-layer curl balance also tests the extra mean-shear harmonics.
    # An average on the doubled periodic extension would wrongly erase the
    # horizontal vorticity components; integrate only the physical half layer.
    omega_hat = curl(model, state[:3])
    omega = physical(omega_hat)
    stretching = np.sum(omega[:, None]*derivatives[:, :3], axis=0)
    curl_rhs_average = half_average(model, curl(model, model.rhs(state)[:3]))
    stretch_average = half_average(model, spectral(stretching))
    rotation_average = half_average(model, model.f*1j*model.k[2]*state[:3])
    viscosity_average = half_average(model, -model.nu*model.k2*omega_hat)
    curl_expected = stretch_average+rotation_average+viscosity_average
    curl_scale = max(norm(stretch_average)+norm(rotation_average)+norm(viscosity_average), 1e-30)
    checks.add(prefix+"/whole_layer_vorticity_balance", norm(curl_rhs_average-curl_expected)/curl_scale, 3e-12,
               actual=curl_rhs_average.tolist(), expected=curl_expected.tolist())

    diagnosis = model.diagnostics(state)
    independent = {
        "kinetic": .5*inner(fields[:3], fields[:3]),
        "thermal_variance": .5*inner(fields[3], fields[3]),
        "potential_moment": potential_moment(model, state),
        "buoyancy_power": heat,
        "barotropic_power": mean_work,
        "barotropic_stress_power": stress_work,
    }
    for key, expected in independent.items():
        checks.add(prefix+"/diagnostic_"+key, abs(diagnosis[key]-expected)/max(abs(expected), 1e-16), 3e-11)
    zeta_bt = physical(curl(model, bt))[2]
    zeta_rate = physical(curl(model, model.rhs(state)[:3]))[2].mean(axis=-1)
    third_rate = float(3*np.mean(zeta_bt[:, :, 0]**2*zeta_rate))
    third_scale = max(float(3*np.mean(abs(zeta_bt[:, :, 0]**2*zeta_rate))), 1e-30)
    checks.add(prefix+"/third_moment_actual_projected_tendency", abs(diagnosis["barotropic_vorticity_third_moment_tendency"]-third_rate)/third_scale, 3e-12)
    third_sources = sum(diagnosis[key] for key in ("barotropic_third_moment_self_advective",
                        "barotropic_third_moment_convective", "barotropic_third_moment_viscous",
                        "barotropic_third_moment_coriolis", "barotropic_third_moment_buoyancy"))
    checks.add(prefix+"/third_moment_source_sum", abs(third_sources-third_rate)/third_scale, 3e-12,
               note="Projected barotropic self-advection is retained; cubic vorticity is not assumed invariant.")
    return state


def analytic_modes(model, checks, prefix):
    x = model.kfund*model.x[:, None, None]
    q = np.pi*model.Nv
    z = model.z[None, None, :]
    W, T = .037, .021
    fields = np.zeros((4,)+model.shape)
    fields[0] = -q/model.kfund*W*np.sin(x)*np.cos(q*z)
    fields[2] = W*np.cos(x)*np.sin(q*z)
    fields[3] = T*np.cos(x)*np.sin(q*z)
    state = model.project(spectral(fields))
    actual = physical(model.tendencies(state)["advective"])[3]
    expected = np.broadcast_to(-.5*q*W*T*np.sin(2*q*z), model.shape)
    checks.add(prefix+"/roll_generates_2Nv_mean_temperature", relative(actual, expected), 2e-12)
    moment_expected = -.25*W*T
    checks.add(prefix+"/roll_potential_conversion", relative(potential_moment(model, model.tendencies(state)["advective"]), moment_expected), 2e-12)
    # The generated mean profile feeds back on the retained convective roll.
    M = -.016
    fields[3] = M*np.sin(2*q*z)
    feedback = physical(model.tendencies(model.project(spectral(fields)))["advective"])[3]
    expected = np.broadcast_to(q*W*M*np.cos(x)*np.sin(q*z), model.shape)
    checks.add(prefix+"/mean_temperature_feeds_back_on_roll", relative(feedback, expected), 3e-12)
    checks.add(prefix+"/saturation_feedback_sign", max(float(np.mean(feedback*np.cos(x)*np.sin(q*z))), 0.), 0.,
               feedback_coefficient=q*W*M)

    # Three independent real basis states span a linear convective wave.
    q = np.pi
    k = model.kfund
    d = k*k+q*q
    basis_fields = []
    for component in range(3):
        value = np.zeros((4,)+model.shape)
        if component == 0:
            value[0] = -q/k*np.sin(x)*np.cos(q*z)
            value[2] = np.cos(x)*np.sin(q*z)
        elif component == 1:
            value[1] = np.sin(x)*np.cos(q*z)
        else:
            value[3] = np.cos(x)*np.sin(q*z)
        basis_fields.append(value)
    basis = [model.project(spectral(value)) for value in basis_fields]
    matrix = np.empty((3, 3))
    leakage = []
    for column, state in enumerate(basis):
        rhs = model.rhs(state, nonlinear=False)
        for row, vector in enumerate(basis):
            matrix[row, column] = np.vdot(vector, rhs).real/np.vdot(vector, vector).real
        reconstruction = sum(matrix[row, column]*basis[row] for row in range(3))
        leakage.append(relative(rhs, reconstruction))
    expected_matrix = np.array([[-model.nu*d, -model.f*q*k/d, k*k/d],
                                [model.f*q/k, -model.nu*d, 0.],
                                [model.G, 0., -model.kappa*d]])
    checks.add(prefix+"/linear_wave_operator", relative(matrix, expected_matrix), 3e-13,
               observed=matrix.tolist(), analytic=expected_matrix.tolist())
    checks.add(prefix+"/linear_wave_subspace", max(leakage), 3e-13)
    if math.isclose(model.nu, model.kappa, rel_tol=1e-14, abs_tol=0):
        root = np.lib.scimath.sqrt((model.G*k*k-model.f**2*q*q)/d)
        exact = np.array([-model.nu*d, -model.nu*d+root, -model.nu*d-root])
        observed = np.linalg.eigvals(matrix)
        error = min(norm(observed-exact[list(order)]) for order in itertools.permutations(range(3)))
        checks.add(prefix+"/Pr1_analytic_linear_dispersion", error/max(norm(exact), 1e-30), 4e-12,
                   eigenvalues=[[float(v.real), float(v.imag)] for v in observed],
                   analytic=[[float(v.real), float(v.imag)] for v in exact])

    # A three-wave barotropic field has a known nonzero vorticity skewness.
    # No assumption of cubic-Casimir conservation under finite truncation.
    y = model.kfund*model.y[None, :, None]
    c = .5
    fields = np.zeros((4,)+model.shape)
    fields[0] = k*(np.sin(y)+c*np.sin(x+y))
    fields[1] = -k*(np.sin(x)+c*np.sin(x+y))
    state = model.project(spectral(fields))
    zeta = physical(curl(model, state[:3]))[2]
    skew = float(np.mean(zeta**3)/np.mean(zeta*zeta)**1.5)
    expected_skew = -3*c/(1+2*c*c)**1.5
    checks.add(prefix+"/signed_vorticity_is_representable", abs(skew-expected_skew), 3e-13,
               observed_skewness=skew, analytic_skewness=expected_skew)
    reflected = physical(curl(model, model.reflect(state)[:3]))[2]
    reflected_skew = float(np.mean(reflected**3)/np.mean(reflected**2)**1.5)
    checks.add(prefix+"/reflection_flips_vorticity_skewness", abs(reflected_skew+skew), 3e-13)


def symmetry_and_quadrature(model, state, checks, prefix):
    opposite = copy_model(model, f=-model.f)
    reflected = model.reflect(state)
    checks.add(prefix+"/reflection_involution", relative(model.reflect(reflected), state), 2e-13)
    checks.add(prefix+"/f_reversal_equivariance", relative(opposite.rhs(reflected), model.reflect(model.rhs(state))), 5e-13)
    checks.add(prefix+"/reflection_energy", abs(norm(reflected)-norm(state))/norm(state), 2e-13)
    finer = copy_model(model, grid_shape=tuple(2*n for n in model.shape))
    finer_state = transfer(model, finer, state)
    expected = transfer(model, finer, model.rhs(state))
    checks.add(prefix+"/quadrature_doubling_rhs", relative(finer.rhs(finer_state), expected), 2e-12,
               original_grid=list(model.shape), doubled_grid=list(finer.shape))
    checks.add(prefix+"/quadrature_doubling_potential", abs(potential_moment(finer, finer_state)-potential_moment(model, state))/max(norm(state[3]), 1e-30), 3e-13)


def rates(model, state, rhs):
    """Independent spectral scalar derivatives, with half-layer M integral."""
    return np.array([np.vdot(state[:3], rhs[:3]).real,
                     np.vdot(state[3], rhs[3]).real,
                     potential_moment(model, rhs)])


def values(model, state):
    return np.array([.5*np.vdot(state[:3], state[:3]).real,
                     .5*np.vdot(state[3], state[3]).real,
                     potential_moment(model, state)])


def evolve(model, state, dt, steps):
    state = state.copy()
    accumulated = np.zeros(3)
    for _ in range(steps):
        state, stages = model.rk4_step(state, dt, return_stages=True)
        accumulated += dt*sum(weight*rates(model, value, rhs) for value, rhs, weight in stages)
        if not np.all(np.isfinite(state)):
            raise FloatingPointError("Non-finite state in bounded time-refinement validation")
    return state, accumulated


def temporal_checks(model, state, dt, checks, prefix):
    # Fixed physical endpoint and 1:2:4 steps.  This exercises nonlinearity;
    # initial production perturbations alone would make this a linear test.
    steps = 24
    solutions, budgets = [], []
    start = values(model, state)
    residuals = []
    for divisor in (1, 2, 4):
        solution, budget = evolve(model, state, dt/divisor, steps*divisor)
        solutions.append(solution)
        budgets.append(budget)
        residuals.append(values(model, solution)-start-budget)
    e12 = relative(solutions[0], solutions[1])
    e24 = relative(solutions[1], solutions[2])
    observed_ratio = e12/max(e24, 1e-30)
    checks.add(prefix+"/time_refinement_accuracy", e12, 1e-5,
               dt=dt, endpoint=steps*dt, full_half=e12, half_quarter=e24,
               observed_error_ratio=observed_ratio)
    # The floor avoids rejecting a correct solver once differences round off.
    checks.add(prefix+"/time_refinement_fourth_order", e24/max(e12, 2e-13), .15,
               expected_ratio=1/16, observed_ratio=observed_ratio)
    norms = np.maximum(np.abs(start), 1e-12)
    checks.add(prefix+"/rk_stage_budget_closure", norm(residuals[0]/norms), 1e-5,
               kinetic_thermal_potential_residual=residuals[0].tolist())
    checks.add(prefix+"/rk_stage_budget_refinement", norm(residuals[1]/norms)/max(norm(residuals[0]/norms), 2e-13), .2)
    checks.add(prefix+"/linear_potential_time_integral", abs(residuals[0][2])/max(abs(start[2]), 1e-12), 3e-12)
    opposite = copy_model(model, f=-model.f)
    mirrored, _ = evolve(opposite, model.reflect(state), dt, steps)
    checks.add(prefix+"/f_reversed_trajectory", relative(mirrored, model.reflect(solutions[0])), 3e-12)
    # A full trajectory must differ from the deliberately linearized equation.
    linear = state.copy()
    for _ in range(steps):
        linear = model.rk4_step(linear, dt, nonlinear=False)
    nonlinear_difference = relative(solutions[0], linear)
    checks.add(prefix+"/nonlinear_feedback_is_active", 0. if nonlinear_difference > 1e-7 else 1., 0.,
               relative_full_vs_linear_difference=nonlinear_difference)
    next_state, stages = model.rk4_step(state, dt, return_stages=True)
    accounted_state, work = model.rk4_step_with_work(state, dt)
    checks.add(prefix+"/work_api_same_trajectory", relative(accounted_state, next_state), 2e-13)
    expected = dt*sum(weight*rates(model, value, rhs) for value, rhs, weight in stages)
    for index, key in enumerate(("kinetic", "thermal_variance", "potential_moment")):
        checks.add(prefix+"/work_api_"+key, abs(work[key]-expected[index])/max(abs(expected[index]), 1e-20), 3e-12)
    bt_expected = dt*sum(weight*np.vdot(value[:2, :, :, 0], rhs[:2, :, :, 0]).real for value, rhs, weight in stages)
    checks.add(prefix+"/work_api_barotropic", abs(work["barotropic_kinetic"]-bt_expected)/max(abs(bt_expected), 1e-20), 3e-12)
    source_sums = {
        "kinetic": work["buoyancy"]+work["viscous"]+work["kinetic_advective_and_coriolis"],
        "thermal_variance": work["thermal_gradient"]+work["thermal_diffusive"]+work["thermal_advective"],
        "potential_moment": work["potential_advective"]+work["potential_diffusive"],
        "barotropic_kinetic": work["barotropic_advective"]+work["barotropic_viscous"],
    }
    for key, expected in source_sums.items():
        checks.add(prefix+"/work_sources_"+key, abs(work[key]-expected)/max(abs(work[key]), 1e-20), 3e-12)


def timed(call, repetitions):
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    for _ in range(repetitions):
        call()
    return {"wall_seconds": (time.perf_counter()-wall_start)/repetitions,
            "cpu_seconds": (time.process_time()-cpu_start)/repetitions}


def estimate_cost(model, args, count):
    state = model.initial_state(seed=args.seed, velocity_rms=args.velocity_rms,
                                temperature_rms=args.temperature_rms)
    # Include the work-accounting API if available in the final engine.
    step_method = getattr(model, "rk4_step_with_work", None)
    if step_method is None:
        def step():
            result, stages = model.rk4_step(state, args.dt, return_stages=True)
            sum(weight*rates(model, value, rhs) for value, rhs, weight in stages)
            return result
    else:
        def step():
            return step_method(state, args.dt)
    step()
    timing = timed(step, 32)
    diagnostic = timed(lambda: model.diagnostics(state), 8)
    plane = timed(lambda: model.planes(state, nxy=args.snapshot_nxy), 4)
    nsteps = math.ceil(args.duration/args.dt)
    ndiagnostics = math.ceil(args.duration/args.diagnostic_interval)+1
    nplanes = math.ceil(args.duration/args.snapshot_interval)+1
    cpu = count*(nsteps*timing["cpu_seconds"] + ndiagnostics*diagnostic["cpu_seconds"] + nplanes*plane["cpu_seconds"])
    wall = count*(nsteps*timing["wall_seconds"] + ndiagnostics*diagnostic["wall_seconds"] + nplanes*plane["wall_seconds"])
    return {"K": model.K, "Nv": model.Nv, "runs": count,
            "steps_per_run": nsteps, "diagnostics_per_run": ndiagnostics,
            "plane_exports_per_run": nplanes, "step": timing,
            "diagnostics": diagnostic, "planes": plane,
            "estimated_cpu_seconds": cpu, "estimated_serial_wall_seconds": wall}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--K", "--k", dest="K", type=int, default=4)
    parser.add_argument("--Nv", "--nv", dest="Nv", type=int, default=2)
    parser.add_argument("--refined-K", "--refined-k", dest="refined_K", type=int, default=6)
    parser.add_argument("--refined-Nv", "--refined-nv", dest="refined_Nv", type=int, default=3)
    parser.add_argument("--kfund", type=float, default=20.)
    parser.add_argument("--nu", type=float, default=1e-5)
    parser.add_argument("--kappa", type=float, default=1e-5)
    parser.add_argument("--G", "--g", dest="G", type=float, required=True)
    parser.add_argument("--f", type=float, default=1.)
    parser.add_argument("--dt", type=float, default=.02)
    parser.add_argument("--duration", type=float, default=200.)
    parser.add_argument("--base-runs", type=int, default=4)
    parser.add_argument("--refined-runs", type=int, default=1)
    parser.add_argument("--cpu-budget-seconds", type=float, default=1800.)
    parser.add_argument("--cost-safety-factor", type=float, default=2.)
    parser.add_argument("--diagnostic-interval", type=float, default=.5)
    parser.add_argument("--snapshot-interval", type=float, default=1.)
    parser.add_argument("--snapshot-nxy", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--velocity-rms", type=float, default=1e-4)
    parser.add_argument("--temperature-rms", type=float, default=1e-5)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    pass_path = args.output/"PASS.json"
    # An earlier successful marker cannot release a changed or failed check.
    pass_path.unlink(missing_ok=True)
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    checks = Checks()
    report = {"purpose": "Independent implementation and bounded-cost release gate",
              "scope": "Nonlinear stress-free fixed-temperature reference layer; not Breeze validation",
              "configuration": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
              "environment": {"hostname": platform.node(), "python": sys.version,
                              "numpy": np.__version__, "platform": platform.platform()},
              "checks": checks.results,
              "source_sha256": {name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                                for name in ("model.py", "validate.py")}}
    try:
        checks.add("configuration/positive_step_and_budget", 0. if args.dt > 0 and args.duration > 0 and args.cpu_budget_seconds > 0 else 1., 0.)
        models = [GalerkinModel(K=args.K, Nv=args.Nv, kfund=args.kfund,
                               nu=args.nu, kappa=args.kappa, G=args.G, f=args.f),
                  GalerkinModel(K=args.refined_K, Nv=args.refined_Nv, kfund=args.kfund,
                               nu=args.nu, kappa=args.kappa, G=args.G, f=args.f)]
        report["models"] = [model.metadata() for model in models]
        for prefix, model in zip(("base", "refined"), models):
            state = algebra_checks(model, checks, prefix)
            analytic_modes(model, checks, prefix)
            symmetry_and_quadrature(model, state, checks, prefix)
            temporal_checks(model, state, args.dt, checks, prefix)
        # An exact refinement of the seeded production IC must preserve its
        # resolved physical field, including the zero initial barotropic mode.
        base_initial = models[0].initial_state(seed=args.seed, velocity_rms=args.velocity_rms,
                                               temperature_rms=args.temperature_rms)
        embedded = transfer(models[0], models[1], base_initial)
        checks.add("initial_condition/refinement_preserves_norm", abs(norm(embedded)-norm(base_initial))/norm(base_initial), 3e-13)
        checks.add("initial_condition/zero_barotropic_velocity", norm(base_initial[:3, :, :, 0]), 1e-20)
        checks.add("initial_condition/zero_mean_temperature_and_shear", norm(base_initial[:, 0, 0, :]), 1e-20)
        cost = [estimate_cost(model, args, count) for model, count in zip(models, (args.base_runs, args.refined_runs))]
        validation_cpu = time.process_time()-cpu_start
        # Rendering, compression, and queue delay are not model-step costs.
        # Reserve a fixed two minutes of CPU for small publication outputs;
        # the campaign wall/CPU monitor must enforce the remaining hard cap.
        estimate = args.cost_safety_factor*sum(item["estimated_cpu_seconds"] for item in cost)+validation_cpu+120.
        report["cost"] = {"measurements": cost, "safety_factor": args.cost_safety_factor,
                          "validation_cpu_seconds": validation_cpu,
                          "publication_cpu_reserve_seconds": 120.,
                          "conservative_total_cpu_seconds": estimate,
                          "cpu_budget_seconds": args.cpu_budget_seconds,
                          "excludes": "scheduler wait; scientific publication must remain on Unity"}
        checks.add("cost/bounded_campaign", estimate, args.cpu_budget_seconds)
    except Exception as error:
        report["exception"] = {"type": type(error).__name__, "message": str(error),
                               "traceback": traceback.format_exc()}
        checks.add("execution/exception", 1., 0.)
    report["passed"] = checks.passed
    report["validation_wall_seconds"] = time.perf_counter()-wall_start
    report["validation_cpu_seconds"] = time.process_time()-cpu_start
    (args.output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    if report["passed"]:
        pass_path.write_text(json.dumps({"passed": True,
                                        "source_sha256": report["source_sha256"],
                                        "configuration": report["configuration"],
                                        "conservative_total_cpu_seconds": report["cost"]["conservative_total_cpu_seconds"],
                                        "validation": "validation.json"}, indent=2)+"\n")
    failures = [item["name"] for item in checks.results if not item["passed"]]
    print(json.dumps({"passed": report["passed"], "checks": len(checks.results),
                      "failures": failures, "report": str(args.output/"validation.json")}), flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
