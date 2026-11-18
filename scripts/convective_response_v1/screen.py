"""Bounded vertical stability and truncation screen. Run on Unity only."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from model import (F_ABS, HEIGHT_M, VerticalSystem, algebra_checks, full_system,
                   read_reference, real_velocity_basis, source_hashes)


def run(source, out, reference_only=False):
    started = time.monotonic()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    n2, provenance = read_reference(source)
    checks = algebra_checks()
    kh = np.geomspace(.5, 64., 33)
    references = {}
    objects = {}
    for nv, nb in ((32, 64), (64, 96)):
        key = f"velocity_{nv}_thermal_{nb+1}"
        system = full_system(nv, nb, n2)
        objects[key] = system
        rows = []
        for k in kh:
            values = system.spectrum(k)
            leading = values[np.argmax(values.real)]
            rows.append({"kh_H": float(k), "growth_over_f": float(leading.real),
                         "frequency_over_f": float(leading.imag)})
        references[key] = rows
        print(key, "max growth", max(r["growth_over_f"] for r in rows), flush=True)
    fine_key = "velocity_64_thermal_97"
    coarse_key = "velocity_32_thermal_65"
    fine = references[fine_key]
    ix = int(np.argmax([r["growth_over_f"] for r in fine]))
    kstar = float(kh[ix])
    growth = fine[ix]["growth_over_f"]
    coarse_growth = references[coarse_key][ix]["growth_over_f"]
    reference_difference = abs(coarse_growth-growth)/max(abs(growth), .01)
    if reference_only:
        report = {"status": "REFERENCE_SCREEN_COMPLETE", "compute_location": "Unity",
                  "provenance": provenance, "code_sha256": source_hashes(), "checks": checks,
                  "reference_scan": references, "selected_kh_H": kstar,
                  "selected_growth_over_f": growth,
                  "reference_growth_relative_difference": reference_difference,
                  "elapsed_seconds": time.monotonic()-started,
                  "release": "No response simulation is authorized by this reference-only screen."}
        (out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps({k: report[k] for k in ("status", "selected_kh_H", "selected_growth_over_f", "reference_growth_relative_difference", "elapsed_seconds")}), flush=True)
        return
    g, selection = real_velocity_basis(objects[fine_key], kstar, 2)
    np.savez_compressed(out/"selected_velocity_basis.npz", g_coeff=g, n2_si=n2,
                        kh_H=kstar, f_abs_s_inv=F_ABS, height_m=HEIGHT_M)
    band = kstar*np.array([.8, 1., 1.2])
    reduced = {}
    for nm in (1, 2):
        for nb in (32, 64, 96):
            system = VerticalSystem(g[:, :nm], nb, n2)
            rows = []
            for k in band:
                vals = system.spectrum(float(k))
                rows.append(float(np.max(vals.real)))
            reduced[f"velocity_{nm}_thermal_{nb+1}"] = rows
    fine_band = [float(np.max(objects[fine_key].spectrum(float(k)).real)) for k in band]
    two = reduced["velocity_2_thermal_97"]
    one = reduced["velocity_1_thermal_97"]
    scale = np.maximum(np.abs(fine_band), .01)
    two_error = float(np.max(np.abs(np.asarray(two)-fine_band)/scale))
    one_two_difference = float(np.max(np.abs(np.asarray(two)-one)/scale))
    thermal_difference = float(np.max(np.abs(np.asarray(reduced["velocity_2_thermal_65"])-two)/scale))
    control = VerticalSystem(g, 96, n2, control=True)
    control_growth = [float(np.max(control.spectrum(float(k)).real)) for k in band]
    gates = {"reference_has_resolved_growth": bool(growth > .001),
             "vertical_reference_growth_difference_below_10pct": bool(reference_difference < .1),
             "two_mode_growth_error_below_10pct": bool(two_error < .1),
             "one_two_growth_difference_below_10pct": bool(one_two_difference < .1),
             "thermal_65_97_difference_below_5pct": bool(thermal_difference < .05),
             "leading_mode_real": bool(not selection[0]["requires_real_imag_pair"]),
             "neutral_control_no_resolved_growth": bool(max(control_growth) < .001),
             "selected_scale_inside_scan": bool(0 < ix < len(kh)-1)}
    report = {"status": "PASS" if all(gates.values()) else "STOP_FOR_MODEL_ADEQUACY",
              "compute_location": "Unity", "model": "Frozen-profile finite-f Boussinesq linear response",
              "provenance": provenance, "code_sha256": source_hashes(), "checks": checks,
              "reference_scan": references, "selected_kh_H": kstar,
              "selected_horizontal_wavelength_m": float(2*np.pi*HEIGHT_M/kstar),
              "selected_growth_over_f": growth, "reference_growth_relative_difference": reference_difference,
              "selected_vertical_functions": selection, "validation_band_kh_H": band.tolist(),
              "fine_reference_growth_on_band": fine_band, "reduced_growth_on_band": reduced,
              "two_mode_relative_growth_error": two_error, "one_two_relative_growth_difference": one_two_difference,
              "thermal_relative_growth_difference": thermal_difference, "gradient_removed_control_growth": control_growth,
              "gates": gates, "elapsed_seconds": time.monotonic()-started,
              "limitations": ["Boussinesq reduction omits the original roughly tenfold density contrast.",
                              "Frozen measured stratification; no saturation or self-maintained convection is modeled.",
                              "Gradient-removed control is not instantaneous bottom-heating switch-off.",
                              "Vertical velocity sponge and global thermal controller are omitted.",
                              "No nonlinear vortex-response experiment is released by a failed adequacy gate."]}
    (out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    if report["status"] == "PASS":
        (out/"PASS.json").write_text(json.dumps({"status": "PASS", "code_sha256": report["code_sha256"]})+"\n")
    print(json.dumps({k: report[k] for k in ("status", "selected_kh_H", "selected_growth_over_f", "gates", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("out")
    ap.add_argument("--reference-only", action="store_true")
    args = ap.parse_args()
    run(args.source, args.out, reference_only=args.reference_only)
