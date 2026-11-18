"""Complete-mode adequacy and conservative-operator release gate; Unity only."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment

from model import (F_ABS, HEIGHT_M, full_system, read_reference,
                   reference_arrays, source_hashes)
from optimized import AugmentedSystem, complete_eigenmode_bases


def frozen_hashes():
    return {k: v for k, v in source_hashes().items()
            if k in ("model.py", "optimized.py", "optimized_screen.py")}


def relative(a, b, floor=1e-14):
    return float(la.norm(a-b)/max(la.norm(b), floor))


def algebra_checks():
    """Checks deliberately include nonorthogonal kinetic mass matrices."""
    rng = np.random.default_rng(2026111804)
    nv = 8
    gc = la.qr(rng.normal(size=(nv, 2)), mode="economic")[0]
    gc = gc@np.array([[1., .23], [0., .83]])
    q = np.arange(1, nv+1)*np.pi
    hc = la.qr(np.column_stack((q[:, None]*gc, rng.normal(size=(nv, 2)))),
               mode="economic")[0]
    hc = hc@(np.eye(4)+.12*np.triu(np.ones((4, 4)), 1))
    stable = np.full(64, .7*F_ABS**2)
    s = AugmentedSystem(gc, hc, 18, stable, thermal_diffusion=False,
                        thermal_relaxation=False, nu=0.)
    checks = {"horizontal_mass_condition": float(np.linalg.cond(s.H)),
              "vertical_mass_condition": float(np.linalg.cond(s.W)),
              "derivative_span_relative": s.derivative_span_error}
    residuals = {}
    for label, kx, ky in (("nonzero", 2.4, -.7), ("zero", 0., 0.)):
        p = s.projector(kx, ky)
        d = s.divergence(kx, ky)
        residuals[label+"_divergence"] = float(la.norm(d@p)/la.norm(d))
        residuals[label+"_idempotence"] = float(la.norm(p@p-p)/la.norm(p))
        residuals[label+"_energy_adjoint"] = float(la.norm(s.M@p-p.conj().T@s.M)/la.norm(s.M))
        y = rng.normal(size=s.nvelocity+s.nb+1)+1j*rng.normal(size=s.nvelocity+s.nb+1)
        y[:s.nvelocity] = p@y[:s.nvelocity]
        dy = s.full_operator(kx, ky)@y
        kp = np.real(np.vdot(y[:s.nvelocity], s.M@dy[:s.nvelocity]))
        bp = np.real(np.vdot(y[s.nvelocity:], dy[s.nvelocity:]))/.7
        energy_scale = (la.norm(y[:s.nvelocity])*la.norm(s.M@dy[:s.nvelocity])
                        + la.norm(y[s.nvelocity:])*la.norm(dy[s.nvelocity:])/.7)
        residuals[label+"_uniform_N2_energy"] = float(abs(kp+bp)/max(energy_scale, 1e-14))
        force = rng.normal(size=s.nvelocity)+1j*rng.normal(size=s.nvelocity)
        pw = np.real(np.vdot(y[:s.nvelocity], s.M@(p@force-force)))
        residuals[label+"_pressure_work"] = float(abs(pw)/max(la.norm(y)*la.norm(s.M@force), 1e-14))
        if label == "zero":
            residuals["zero_horizontal_preserved"] = relative(p[:2*s.nh, :2*s.nh], np.eye(2*s.nh))
            residuals["zero_vertical_removed"] = float(la.norm(p[2*s.nh:]))
    # Uniform N2, one sine velocity mode, finite cosine buoyancy rank.
    dispersion = []
    for n2 in (.7, -2.):
        ss = AugmentedSystem(np.ones((1, 1)), np.ones((1, 1)), 16,
                             np.full(64, n2*F_ABS**2), thermal_diffusion=False,
                             thermal_relaxation=False, nu=0.)
        kh = 3.
        omega2 = (np.pi**2+n2*float(np.sum(ss.C**2))*kh**2)/(np.pi**2+kh**2)
        expected = np.array([0., np.sqrt(complex(-omega2)), -np.sqrt(complex(-omega2))])
        ev = ss.spectrum(kh)
        error = max(float(np.min(abs(ev-v))) for v in expected)/max(1., np.max(abs(expected)))
        residuals[f"finite_thermal_dispersion_N2_{n2}"] = error
        dispersion.append({"N2_over_f2": n2, "overlap_squared": float(np.sum(ss.C**2)),
                           "omega_squared": float(omega2), "relative_error": error})
    # Independently reconstruct physical fields to check the dissipative work.
    sd = AugmentedSystem(gc, hc, 18, np.linspace(-.2, .7, 64)*F_ABS**2)
    kx, ky = 2.4, -.7
    y = rng.normal(size=sd.nvelocity+sd.nb+1)+1j*rng.normal(size=sd.nvelocity+sd.nb+1)
    y[:sd.nvelocity] = sd.projector(kx, ky)@y[:sd.nvelocity]
    nh, vel = sd.nh, sd.nvelocity
    ux, uy, ww, bb = y[:nh], y[nh:2*nh], y[2*nh:vel], y[vel:]
    z, weight = sd.z, sd.weight
    cos = np.sqrt(2.)*np.cos(z[:, None]*q)
    sin = np.sqrt(2.)*np.sin(z[:, None]*q)
    dh = -sin@(q[:, None]*hc)
    dg = cos@(q[:, None]*gc)
    rz = np.arange(sd.nb+1)*np.pi
    dbeta = -np.sqrt(2.)*np.sin(z[:, None]*rz)*rz
    fields = (sd.h@ux, sd.h@uy, sd.g@ww)
    derivatives = (dh@ux, dh@uy, dg@ww)
    bz, dbz = sd.beta@bb, dbeta@bb
    k2 = kx*kx+ky*ky
    viscous = sd.nu*sum(float(np.sum(weight*(k2*abs(v)**2+abs(dv)**2)))
                        for v, dv in zip(fields, derivatives))
    buoyancy = float(np.real(np.sum(weight*np.conj(fields[2])*bz)))
    n2, kappa, relaxation = reference_arrays(sd.n2_si, z)
    thermal_gradient = -float(np.real(np.sum(weight*np.conj(bz)*n2*fields[2])))
    thermal_loss = float(np.sum(weight*(kappa*(k2*abs(bz)**2+abs(dbz)**2)+relaxation*abs(bz)**2)))
    dy = sd.full_operator(kx, ky)@y
    kinetic = float(np.real(np.vdot(y[:vel], sd.M@dy[:vel])))
    variance = float(np.real(np.vdot(bb, dy[vel:])))
    residuals["independent_physical_kinetic_work"] = abs(kinetic-buoyancy+viscous)/max(abs(kinetic)+abs(buoyancy)+viscous, 1e-14)
    residuals["independent_physical_thermal_work"] = abs(variance-thermal_gradient+thermal_loss)/max(abs(variance)+abs(thermal_gradient)+thermal_loss, 1e-14)
    # Local mean strain: the two modal work evaluations have opposite signs.
    grad = np.array([[.3, -.8], [.2, -.3]])
    u = np.stack((ux, uy))
    stress = np.real(np.einsum("ia,ab,jb->ij", u.conj(), sd.H, u))
    pbt = float(np.sum(stress*grad))
    strain_force = -grad@u
    strain_work = float(np.real(sum(np.vdot(u[i], sd.H@strain_force[i]) for i in range(2))))
    residuals["equal_opposite_mean_strain_work"] = abs(strain_work+pbt)/max(abs(strain_work)+abs(pbt), 1e-14)
    checks.update(residuals=residuals, dispersion=dispersion,
                  kinetic_work_components={"kinetic": kinetic, "buoyancy": buoyancy, "viscous_loss": viscous},
                  thermal_work_components={"variance": variance, "gradient": thermal_gradient, "loss": thermal_loss})
    assert max(residuals.values()) < 1e-9, checks
    assert s.derivative_span_error < 1e-10
    return checks


def eigenstate_check(reference, reduced, kh, selected):
    values, states, expansion = reference.spectrum(kh, vectors=True)
    rows = []
    for entry in selected:
        wanted = complex(entry["eigenvalue_real"], entry["eigenvalue_imag"])
        i = int(np.argmin(abs(values-wanted)))
        value = values[i]
        full = expansion@states[:, i]
        nr = reference.m
        ux = reference.h_coeff@full[:nr]
        uy = reference.h_coeff@full[nr:2*nr]
        ww = reference.g_coeff@full[2*nr:3*nr]
        bb = full[3*nr:]
        a = reduced.Hinv@(reduced.h_coeff.T@ux)
        b = reduced.Hinv@(reduced.h_coeff.T@uy)
        w = reduced.Winv@(reduced.g_coeff.T@ww)
        contained = np.r_[a, b, w, bb]
        original_velocity = np.r_[ux, uy, ww]
        reconstructed_velocity = np.r_[reduced.h_coeff@a, reduced.h_coeff@b, reduced.g_coeff@w]
        containment = relative(reconstructed_velocity, original_velocity)
        residual = la.norm(reduced.full_operator(kh)@contained-value*contained)/max(abs(value)*la.norm(contained), 1e-14)
        ev = reduced.spectrum(kh)
        value_error = float(np.min(abs(ev-value))/max(abs(value), .01))
        rows.append({"reference_eigenvalue_real": float(value.real), "reference_eigenvalue_imag": float(value.imag),
                     "velocity_containment_relative": containment,
                     "retained_eigenstate_residual_relative": float(residual),
                     "nearest_eigenvalue_relative_difference": value_error})
    return rows


def run(source, reference_report, out):
    started = time.monotonic()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    hashes = frozen_hashes()
    old = json.loads(Path(reference_report).read_text())
    n2, provenance = read_reference(source)
    assert old["status"] == "REFERENCE_SCREEN_COMPLETE"
    assert old["provenance"]["source_sha256"] == provenance["source_sha256"]
    assert old["code_sha256"]["model.py"] == hashes["model.py"]
    kstar = float(old["selected_kh_H"])
    band = np.array([.8, 1., 1.2])*kstar
    checks = algebra_checks()
    fine, coarse = full_system(64, 96, n2), full_system(32, 64, n2)
    reference_growth = np.array([np.max(fine.spectrum(k).real) for k in band], dtype=float)
    coarse_growth = np.array([np.max(coarse.spectrum(k).real) for k in band], dtype=float)
    assert abs(reference_growth[1]-old["selected_growth_over_f"]) < 1e-10
    scale = np.maximum(abs(reference_growth), .01)
    reference_error = float(np.max(abs(reference_growth-coarse_growth)/scale))
    bases = complete_eigenmode_bases(fine, kstar, counts=(1, 2))
    spectra, controls, eigenstates, errors, thermal, reflection, dimensions = {}, {}, {}, {}, {}, {}, {}
    all_controls = {}
    for modes in (1, 2):
        basis = bases[modes]
        g, h = basis["g_coeff"], basis["h_coeff"]
        dimensions[str(modes)] = {"selected_convective_modes": modes, "w_functions": int(g.shape[1]),
                                  "horizontal_functions": int(h.shape[1]), "selected_modes": basis["selected_modes"]}
        for nb in (32, 64, 96):
            s = AugmentedSystem(g, h, nb, n2)
            key = f"modes_{modes}_thermal_{nb+1}"
            spectra[key] = [float(np.max(s.spectrum(k).real)) for k in band]
            ctrl_nb = AugmentedSystem(g, h, nb, n2, control=True)
            all_controls[key] = [float(np.max(ctrl_nb.spectrum(k).real)) for k in band]
        s = AugmentedSystem(g, h, 96, n2)
        ctrl = AugmentedSystem(g, h, 96, n2, control=True)
        controls[str(modes)] = [float(np.max(ctrl.spectrum(k).real)) for k in band]
        eigenstates[str(modes)] = eigenstate_check(fine, s, kstar, basis["selected_modes"])
        errors[str(modes)] = float(np.max(abs(np.asarray(spectra[f"modes_{modes}_thermal_97"])-reference_growth)/scale))
        thermal[str(modes)] = {"33_to_65": float(np.max(abs(np.asarray(spectra[f"modes_{modes}_thermal_33"])-spectra[f"modes_{modes}_thermal_65"])/scale)),
                               "65_to_97": float(np.max(abs(np.asarray(spectra[f"modes_{modes}_thermal_65"])-spectra[f"modes_{modes}_thermal_97"])/scale))}
        rows = []
        for k in band:
            plus, minus = s.spectrum(k, f=1.), s.spectrum(k, f=-1.)
            cost = abs(plus[:, None]-minus[None, :])
            ii, jj = linear_sum_assignment(cost)
            spectral_error = float(np.max(cost[ii, jj])/max(1., np.max(abs(plus))))
            transform = np.r_[-np.ones(s.nh), np.ones(s.nh+s.nw+s.nb+1)]
            ap, am = s.full_operator(k, f=1.), s.full_operator(-k, f=-1.)
            operator_error = float(la.norm(am*transform[None, :]-transform[:, None]*ap)/max(la.norm(ap), 1.))
            rows.append({"kh_H": float(k), "matched_spectrum_relative": spectral_error,
                         "reflected_operator_relative": operator_error})
        reflection[str(modes)] = rows
        checks[f"basis_{modes}_derivative_span"] = s.derivative_span_error
    one_two = float(np.max(abs(np.asarray(spectra["modes_1_thermal_97"])-spectra["modes_2_thermal_97"])/scale))
    thermal_candidates = {}
    for nb in (32, 64, 96):
        rows = []
        for modes in (1, 2):
            key, refkey = f"modes_{modes}_thermal_{nb+1}", f"modes_{modes}_thermal_97"
            error = float(np.max(abs(np.asarray(spectra[key])-spectra[refkey])/scale))
            ctrlscale = np.maximum(abs(np.asarray(all_controls[refkey])), .01)
            ctrlerror = float(np.max(abs(np.asarray(all_controls[key])-all_controls[refkey])/ctrlscale))
            rows.append({"modes": modes, "heated_relative_error": error,
                         "control_relative_error": ctrlerror,
                         "control_max_growth": max(all_controls[key])})
        thermal_candidates[str(nb)] = rows
    selected_nb = next((nb for nb in (32, 64, 96)
                       if all(max(r["heated_relative_error"], r["control_relative_error"]) < .05
                              and r["control_max_growth"] < .001 for r in thermal_candidates[str(nb)])), 96)
    eig_ok = all(max(row[k] for k in ("velocity_containment_relative", "retained_eigenstate_residual_relative", "nearest_eigenvalue_relative_difference")) < 1e-7
                 for rows in eigenstates.values() for row in rows)
    mirror_ok = all(max(row["matched_spectrum_relative"], row["reflected_operator_relative"]) < 1e-9
                    for rows in reflection.values() for row in rows)
    gates = {"reference_has_positive_growth": bool(reference_growth[1] > .001),
             "full_reference_band_convergence_below_10pct": bool(reference_error < .1),
             "one_mode_band_growth_error_below_20pct": bool(errors["1"] < .2),
             "two_mode_band_growth_error_below_20pct": bool(errors["2"] < .2),
             "one_two_band_growth_difference_below_20pct": bool(one_two < .2),
             "thermal_65_97_difference_below_5pct": bool(max(t["65_to_97"] for t in thermal.values()) < .05),
             "full_eigenstate_containment_and_growth": bool(eig_ok),
             "reflection_operator_and_spectrum": bool(mirror_ok),
             "gradient_removed_control_growth_below_001": bool(max(max(v) for v in controls.values()) < .001),
             "frozen_source_unchanged": bool(hashes == frozen_hashes())}
    report = {"status": "PASS" if all(gates.values()) else "STOP_FOR_MODEL_ADEQUACY",
              "compute_location": "Unity", "provenance": provenance, "code_sha256": hashes,
              "reference_report_sha256": hashlib.sha256(Path(reference_report).read_bytes()).hexdigest(),
              "checks": checks, "selected_kh_H": kstar,
              "selected_wavelength_m": float(2*np.pi*HEIGHT_M/kstar), "band_kh_H": band.tolist(),
              "reference_growth_over_f": reference_growth.tolist(), "coarse_reference_growth_over_f": coarse_growth.tolist(),
              "full_reference_relative_difference": reference_error,
              "basis_dimensions": dimensions, "growth_over_f": spectra,
              "relative_growth_errors": errors, "one_two_relative_difference": one_two,
              "thermal_refinement_relative": thermal, "eigenstate_checks": eigenstates,
              "thermal_basis_candidates": thermal_candidates, "selected_thermal_nb": selected_nb,
              "reflection_checks": reflection, "gradient_removed_control_growth_over_f": controls,
              "gates": gates, "elapsed_seconds": time.monotonic()-started,
              "limitations": ["Frozen measured profile; finite-window linear response only.",
                              "Boussinesq approximation omits the original density contrast.",
                              "One/two selected convective modes retain their required horizontal poloidal/toroidal shapes.",
                              "Gradient-removed control is not switching off heating with an unchanged profile.",
                              "No nonlinear or saturated convection claim follows from this gate."]}
    if all(gates.values()):
        np.savez_compressed(out/"basis.npz", g1=bases[1]["g_coeff"], h1=bases[1]["h_coeff"],
                            g2=bases[2]["g_coeff"], h2=bases[2]["h_coeff"], n2_si=n2,
                            kh_H=kstar, thermal_nb=selected_nb, reference_nv=64,
                            f_abs_s_inv=F_ABS, height_m=HEIGHT_M)
        report["basis_sha256"] = hashlib.sha256((out/"basis.npz").read_bytes()).hexdigest()
        report["individual_basis_sha256"] = {}
        for modes in (1, 2):
            path = out/f"basis_m{modes}.npz"
            np.savez_compressed(path, g_coeff=bases[modes]["g_coeff"], h_coeff=bases[modes]["h_coeff"],
                                n2_si=n2, thermal_nb=selected_nb, khstar=kstar, kh_H=kstar,
                                source_sha256=provenance["source_sha256"])
            report["individual_basis_sha256"][str(modes)] = hashlib.sha256(path.read_bytes()).hexdigest()
    (out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    if all(gates.values()):
        (out/"PASS.json").write_text(json.dumps({k: report[k] for k in
            ("status", "code_sha256", "provenance", "reference_report_sha256", "basis_sha256",
             "individual_basis_sha256", "selected_thermal_nb", "gates")}, indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ("status", "relative_growth_errors", "thermal_refinement_relative", "gates", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("reference_report")
    parser.add_argument("out")
    args = parser.parse_args()
    run(args.source, args.reference_report, args.out)
