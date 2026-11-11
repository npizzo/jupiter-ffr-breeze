#!/usr/bin/env python3
"""Fixed-region vorticity diagnostics; execute only in an allocated Unity job."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
PUBLICATION = REPO / "scripts/gamma2026_publication_v2/publish.py"
EXPECTED_PUBLICATION_SHA = "7570daa965e55daace0e1c01e03a541989910d6fb6eaba6acc4cb96dcc46e374"
EXPECTED_CORE_SHA = "4f42089927e37b6ecb75de4d4b0bf337930550b2a7f0b83cfa9f89fb5bbbc6cc"
MOVIES = "https://npizzo.github.io/jupiter-ffr-breeze/#gamma75"
DISJOINT = ("core_lt1", "annulus_1_3", "annulus_3_4", "annulus_4_5", "exterior_ge5")
CUMULATIVE = {"interior_lt3": DISJOINT[:2], "interior_lt4": DISJOINT[:3], "trap_lt5": DISJOINT[:4]}
REGIONS = DISJOINT + tuple(CUMULATIVE)
DEFINITIONS = {"core_lt1": "0 <= r/Lgamma < 1", "annulus_1_3": "1 <= r/Lgamma < 3",
               "annulus_3_4": "3 <= r/Lgamma < 4", "annulus_4_5": "4 <= r/Lgamma < 5",
               "exterior_ge5": "r/Lgamma >= 5, within the periodic computational square",
               "interior_lt3": "r/Lgamma < 3", "interior_lt4": "r/Lgamma < 4", "trap_lt5": "r/Lgamma < 5"}
ADDITIVE = ("cell_count", "area_m2", "signed_circulation_m2_s", "positive_circulation_m2_s",
            "negative_circulation_magnitude_m2_s", "relative_enstrophy_integral_m2_s2")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def masks(radius, scale=1.0):
    r = radius
    return dict(zip(REGIONS, (r < scale, (r >= scale) & (r < 3 * scale), (r >= 3 * scale) & (r < 4 * scale),
        (r >= 4 * scale) & (r < 5 * scale), r >= 5 * scale, r < 3 * scale, r < 4 * scale, r < 5 * scale)))


def region_stats(zeta, x, y, radius, mask, cell_area, global_values):
    """All radii/centroids are fixed-origin Eulerian moments, never vortex cores."""
    q, r, xx, yy = zeta[mask], radius[mask], x[mask], y[mask]
    count = int(q.size)
    integral = 0.5 * float(np.sum(q * q)) * cell_area
    result = {"cell_count": count, "area_m2": count * cell_area,
        "signed_circulation_m2_s": float(np.sum(q)) * cell_area,
        "mean_vorticity_s_inv": float(np.mean(q)) if count else None,
        "relative_enstrophy_integral_m2_s2": integral,
        "relative_enstrophy_fraction_of_domain": integral / global_values["relative_z"] if global_values["relative_z"] else None}
    for name, sign in (("positive", 1), ("negative", -1)):
        weight = np.maximum(sign * q, 0)
        total = float(np.sum(weight))
        circulation = total * cell_area
        key = "positive_circulation_m2_s" if sign == 1 else "negative_circulation_magnitude_m2_s"
        result[key] = circulation
        result[f"{name}_circulation_fraction_of_domain_sign"] = circulation / global_values[name] if global_values[name] else None
        result[f"{name}_mean_radius_Lgamma"] = float(np.sum(weight * r) / total) if total else None
        result[f"{name}_rms_radius_Lgamma"] = float(np.sqrt(np.sum(weight * r * r) / total)) if total else None
        result[f"{name}_centroid_x_Lgamma"] = float(np.sum(weight * xx) / total) if total else None
        result[f"{name}_centroid_y_Lgamma"] = float(np.sum(weight * yy) / total) if total else None
    a, b = result["negative_mean_radius_Lgamma"], result["positive_mean_radius_Lgamma"]
    result["radial_separation_Lgamma"] = a - b if a is not None and b is not None else None
    return result


def measure(zeta, xx, yy, radius, cell_area, *, mask_radius=None, mask_scale=1.0):
    if zeta.shape != xx.shape or zeta.shape != yy.shape or zeta.shape != radius.shape:
        raise ValueError("Inconsistent regional geometry")
    if not all(np.isfinite(a).all() for a in (zeta, xx, yy, radius)) or not math.isfinite(cell_area) or cell_area <= 0:
        raise ValueError("Nonfinite field or invalid cell area")
    if np.any(radius < 0):
        raise ValueError("Negative radius")
    physical_radius = radius if mask_radius is None else mask_radius
    if physical_radius.shape != radius.shape or not np.isfinite(physical_radius).all() or np.any(physical_radius < 0) or not math.isfinite(mask_scale) or mask_scale <= 0:
        raise ValueError("Invalid physical mask geometry")
    if not np.allclose(physical_radius / mask_scale, radius, rtol=32 * np.finfo(float).eps, atol=0):
        raise ValueError("Mask and reported radii use inconsistent scales")
    mm = masks(physical_radius, mask_scale)
    counts = sum(mm[name].astype(np.int8) for name in DISJOINT)
    if not np.all(counts == 1):
        raise ValueError("Fixed regions do not partition the domain")
    for name, parts in CUMULATIVE.items():
        if not np.array_equal(mm[name], np.logical_or.reduce([mm[p] for p in parts])):
            raise ValueError("Cumulative mask identity failed")
    totals = {"relative_z": 0.5 * float(np.sum(zeta * zeta)) * cell_area,
              "positive": float(np.sum(np.maximum(zeta, 0))) * cell_area,
              "negative": float(np.sum(np.maximum(-zeta, 0))) * cell_area}
    values = {name: region_stats(zeta, xx, yy, radius, mask, cell_area, totals) for name, mask in mm.items()}
    whole = region_stats(zeta, xx, yy, radius, np.ones(zeta.shape, dtype=bool), cell_area, totals)
    # Pairwise NumPy reductions plus five-region recombination. For signed
    # cancellation, use absolute circulation throughput, not near-zero net flow.
    tolerance = np.finfo(float).eps * (64 + 8 * math.ceil(math.log2(max(2, zeta.size))))
    worst_scaled = 0.0
    for target, parts in [(whole, DISJOINT)] + [(values[k], v) for k, v in CUMULATIVE.items()]:
        for key in ADDITIVE:
            observed = sum(values[p][key] for p in parts)
            expected = target[key]
            scale = abs(expected) + sum(abs(values[p][key]) for p in parts)
            if key == "signed_circulation_m2_s":
                scale = totals["positive"] + totals["negative"]
            delta = abs(observed - expected)
            if delta > tolerance * scale + np.finfo(float).tiny:
                raise ValueError(f"Region partition identity failed: {key}")
            worst_scaled = max(worst_scaled, delta / scale if scale else 0.0)
        for sign, circulation in (("positive", "positive_circulation_m2_s"),
                                  ("negative", "negative_circulation_magnitude_m2_s")):
            for moment in ("mean_radius_Lgamma", "centroid_x_Lgamma", "centroid_y_Lgamma", "rms_radius_Lgamma"):
                key = f"{sign}_{moment}"
                power = 2 if moment == "rms_radius_Lgamma" else 1
                terms = [values[p][circulation] * values[p][key] ** power for p in parts if values[p][key] is not None]
                expected = target[circulation] * target[key] ** power if target[key] is not None else 0.0
                scale = abs(expected) + sum(abs(v) for v in terms)
                # Centroid numerators can cancel within every subregion too.
                # Bound reduction error with absolute coordinate throughput.
                coordinate = xx if moment == "centroid_x_Lgamma" else yy if moment == "centroid_y_Lgamma" else radius
                scale = max(scale, totals[sign] * float(np.max(np.abs(coordinate))) ** power)
                delta = abs(sum(terms) - expected)
                if delta > tolerance * scale + np.finfo(float).tiny:
                    raise ValueError(f"Weighted regional moment identity failed: {key}")
                worst_scaled = max(worst_scaled, delta / scale if scale else 0.0)
    for row in values.values():
        if any(v is not None and not math.isfinite(v) for v in row.values()):
            raise ValueError("Nonfinite regional statistic")
        for key, value in row.items():
            if "fraction_of" in key and value is not None and not -tolerance <= value <= 1 + tolerance:
                raise ValueError("Regional fraction outside [0,1]")
        scale = row["positive_circulation_m2_s"] + row["negative_circulation_magnitude_m2_s"]
        if abs(row["signed_circulation_m2_s"] - row["positive_circulation_m2_s"] + row["negative_circulation_magnitude_m2_s"]) > tolerance * scale + np.finfo(float).tiny:
            raise ValueError("Signed circulation identity failed")
    return values, {"identity_relative_tolerance": tolerance, "max_scaled_additive_residual": worst_scaled}


def analytical_checks():
    """Exact discrete constructed fields, not a continuum quadrature assertion."""
    edge = np.array([0., 1., 3., 4., 5., 6.])
    mm = masks(edge)
    for name, expected in zip(DISJOINT, ([0], [1], [2], [3], [4, 5])):
        if list(np.flatnonzero(mm[name])) != expected:
            raise AssertionError("Boundary inclusivity test failed")
    x = np.array([.5, -.5, 2., -2., 3.5, -3.5, 4.5, -4.5, 5.5, -5.5])
    y = np.zeros(10)
    q = np.array([2., -1., 3., -2., 4., -3., 5., -4., 6., -5.])
    values, report = measure(q, x, y, np.abs(x), 2.)
    physical, _ = measure(q, x, y, np.abs(x), 2., mask_radius=2.5 * np.abs(x), mask_scale=2.5)
    if physical != values:
        raise AssertionError("Physical-threshold/normalized-moment equivalence failed")
    for i, (name, radius, zint) in enumerate(zip(DISJOINT, (.5, 2., 3.5, 4.5, 5.5), (5., 13., 25., 41., 61.))):
        row = values[name]
        expected = {"cell_count": 2, "area_m2": 4., "positive_circulation_m2_s": 2. * (i + 2),
            "negative_circulation_magnitude_m2_s": 2. * (i + 1), "signed_circulation_m2_s": 2.,
            "relative_enstrophy_integral_m2_s2": zint, "relative_enstrophy_fraction_of_domain": zint / 145.,
            "positive_mean_radius_Lgamma": radius, "negative_mean_radius_Lgamma": radius,
            "positive_centroid_x_Lgamma": radius, "negative_centroid_x_Lgamma": -radius,
            "radial_separation_Lgamma": 0.}
        for key, value in expected.items():
            if not math.isclose(row[key], value, rel_tol=32 * np.finfo(float).eps, abs_tol=1e-14):
                raise AssertionError(f"Constructed sign-field test failed: {name}/{key}")
    uniform, _ = measure(np.full(10, 2.), x, y, np.abs(x), 2.)
    for name, expected in zip(CUMULATIVE, (1.25, 2., 2.625)):
        if uniform[name]["positive_mean_radius_Lgamma"] != expected:
            raise AssertionError("Constant-annulus radial moment failed")
    for row in uniform.values():
        if row["mean_vorticity_s_inv"] != 2. or row["negative_mean_radius_Lgamma"] is not None or row["radial_separation_Lgamma"] is not None:
            raise AssertionError("Absent-sign convention failed")
    zero, _ = measure(np.zeros(10), x, y, np.abs(x), 2.)
    for row in zero.values():
        if row["relative_enstrophy_fraction_of_domain"] is not None or row["positive_mean_radius_Lgamma"] is not None:
            raise AssertionError("Zero-field undefined-moment convention failed")
    sparse, _ = measure(np.ones(1), np.array([.5]), np.zeros(1), np.array([.5]), 2.)
    if sparse["annulus_1_3"]["mean_vorticity_s_inv"] is not None:
        raise AssertionError("Empty-mask convention failed")
    # Independent two-dimensional fixture: the annulus 3<=r<4 has20 lattice
    # points, with radii3(4points),sqrt10(8points),sqrt13(8points).
    gx, gy = np.meshgrid(np.arange(-6., 7.), np.arange(-6., 7.))
    gr = np.hypot(gx, gy)
    annulus = (gr >= 3) & (gr < 4)
    fixture, _ = measure(np.where(annulus, 3., 0.), gx, gy, gr, 1.)
    expected_radius = (12 + 8 * math.sqrt(10) + 8 * math.sqrt(13)) / 20
    row = fixture["annulus_3_4"]
    for key, expected in (("positive_circulation_m2_s", 60.), ("relative_enstrophy_integral_m2_s2", 90.),
                          ("positive_centroid_x_Lgamma", 0.), ("positive_centroid_y_Lgamma", 0.),
                          ("positive_mean_radius_Lgamma", expected_radius)):
        if not math.isclose(row[key], expected, rel_tol=32 * np.finfo(float).eps, abs_tol=1e-14):
            raise AssertionError(f"Two-dimensional constant-annulus test failed: {key}")
    rotated, _ = measure(-q, -y, x, np.abs(x), 2.)
    for name in REGIONS:
        if rotated[name]["positive_circulation_m2_s"] != values[name]["negative_circulation_magnitude_m2_s"]:
            raise AssertionError("Sign reversal circulation failed")
        if not math.isclose(rotated[name]["positive_centroid_y_Lgamma"], values[name]["negative_centroid_x_Lgamma"], rel_tol=1e-14, abs_tol=1e-14):
            raise AssertionError("Rotated sign-weighted centroid failed")
    return {"status": "PASS", "tests": ["exact half-open radial boundaries", "constructed two-sign constant annuli",
        "physical-threshold/normalized-moment equivalence", "constant positive field moments", "two-dimensional constant annulus", "rotation and sign reversal",
        "absent sign", "zero field", "empty masks", "partition/cumulative integral and weighted-moment identities"], **report}


def observer_case(publication, root, name, audit, output):
    config_path = REPO / "configs/gamma2026_v2" / f"{name}.toml"
    config_sha = sha(config_path)
    config = tomllib.loads(config_path.read_text())
    metadata_path = root / name / "metadata.toml"
    metadata = tomllib.loads(metadata_path.read_text())
    for key, value in config.items():
        if metadata.get(key) != value:
            raise ValueError(f"Frozen configuration mismatch: {name}/{key}")
    if metadata["source_hash"] != EXPECTED_CORE_SHA or metadata["trap_Lgamma"] != 5.0:
        raise ValueError("Wrong frozen model or trap definition")
    if metadata["model"] != "unforced barotropic QG" or metadata["output_dimension"] != 2 or metadata["geophysicalflows_commit"] != "51ca98a0d471ce4ef8bbba8b73daf2bad54d53fb":
        raise ValueError("Wrong physical model or pinned solver revision")
    case = publication.case_analysis(root, name, audit)
    initial = case["initial"]
    xm, ym = np.meshgrid(initial["x_m"], initial["y_m"])
    physical_radius = np.hypot(xm, ym)
    xx, yy, radius = xm / case["Lgamma"], ym / case["Lgamma"], physical_radius / case["Lgamma"]
    area = float(initial["x_m"][1] - initial["x_m"][0]) * float(initial["y_m"][1] - initial["y_m"][0])
    rows, identity_reports = [], []
    for i, path in enumerate(case["files"]):
        field = publication.read_frame(path, metadata["n"])
        if sha(path) != case["summary"]["frames"][i]["sha256"]:
            raise ValueError("Snapshot changed after native diagnostic gate")
        if not np.array_equal(field["x_m"], initial["x_m"]) or not np.array_equal(field["y_m"], initial["y_m"]):
            raise ValueError("Coordinates changed")
        values, identities = measure(field["vorticity_s_inv"], xx, yy, radius, area,
                                     mask_radius=physical_radius, mask_scale=case["Lgamma"])
        identity_reports.append(identities)
        time = field["time_seconds"] / publication.DAY
        if abs(time - i) > 1e-7:
            raise ValueError("Require all76 daily fields")
        # Independently verify this observer's trap moments against the frozen
        # publication, not merely its own partition identities.
        for key in ("positive_mean_radius_Lgamma", "negative_mean_radius_Lgamma", "radial_separation_Lgamma"):
            if not math.isclose(values["trap_lt5"][key], case["rows"][i][key], rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Observer disagrees with existing trap statistic")
        for region, value in values.items():
            rows.append({"case": name, "grid": metadata["n"], "kind": metadata["kind"], "seed": metadata["seed"],
                "time_days": time, "region": region, **value,
                "pv_snapshot_crosscheck_passed": case["rows"][i]["pv_enstrophy_crosscheck_passed"]})
    for path, key in ((metadata_path, "metadata_sha256"), (root / name / "diagnostics.csv", "simulation_diagnostics_sha256")):
        if sha(path) != case["summary"][key]:
            raise ValueError("Input provenance changed during observer")
    if sha(config_path) != config_sha:
        raise ValueError("Frozen configuration changed during observer")
    csv_path = output / f"{name}_regions.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    case["region_rows"] = rows
    case["region_summary"] = {"case": name, "grid": metadata["n"], "kind": metadata["kind"], "seed": metadata["seed"],
        "Lgamma_m": case["Lgamma"], "frames": len(case["files"]), "csv": csv_path.name,
        "configuration_file_sha256": config_sha,
        "initial": {row["region"]: row for row in rows[:len(REGIONS)]}, "final": {row["region"]: row for row in rows[-len(REGIONS):]},
        "snapshot_crosschecks": case["summary"]["snapshot_crosschecks"],
        "max_scaled_additive_residual": max(r["max_scaled_additive_residual"] for r in identity_reports),
        "native_audit_summary": f"native_audit/{name}_summary.json"}
    return case


def plot(cases, output):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    panels = [(region, "radial_separation_Lgamma", f"r < {limit} Lγ: negative − positive mean radius / Lγ")
              for region, limit in zip(CUMULATIVE, (3, 4, 5))]
    panels += [("core_lt1", "positive_circulation_fraction_of_domain_sign", "r < Lγ: fraction of domain positive circulation"),
               ("annulus_4_5", "positive_circulation_fraction_of_domain_sign", "4–5 Lγ: fraction of domain positive circulation"),
               ("exterior_ge5", "relative_enstrophy_fraction_of_domain", "r ≥ 5 Lγ: fraction of domain relative enstrophy")]
    palette = plt.get_cmap("tab10")
    color_pairs = sorted({(case["metadata"]["n"], case["metadata"]["seed"]) for case in cases})
    for case in cases:
        m = case["metadata"]
        label = f"{m['kind']} {m['n']}² s{m['seed']}"
        for ax, (region, key, title) in zip(axes.flat, panels):
            rows = [r for r in case["region_rows"] if r["region"] == region]
            ax.plot([r["time_days"] for r in rows], [np.nan if r[key] is None else r[key] for r in rows],
                    label=label, color=palette(color_pairs.index((m["n"], m["seed"]))),
                    linestyle="-" if m["kind"] == "polar" else "--", linewidth=1.2)
            ax.set_title(title, fontsize=10); ax.set_xlabel("Earth day"); ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("Fixed physical masks: Eulerian vorticity moments, not tracked vortices or cascade fluxes", fontsize=13)
    fig.savefig(output / "gamma_regions.png", dpi=140); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", nargs="+", default=["1024:11", "1024:22", "1024:33"])
    parser.add_argument("--self-test-only", action="store_true", help="Run allocated analytical checks without reading model fields")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("All tests and regional analysis require a Unity allocation")
    if args.output.exists():
        raise ValueError("Use a fresh output directory")
    if sha(PUBLICATION) != EXPECTED_PUBLICATION_SHA or sha(REPO / "scripts/gamma2026_v3/core.jl") != EXPECTED_CORE_SHA:
        raise ValueError("Frozen publication/model source changed")
    source_paths = [Path(__file__).resolve(), PUBLICATION, REPO / "scripts/gamma2026_v3/core.jl", REPO / "slurm/gamma2026-regions.sbatch"]
    source_hashes = {str(p.relative_to(REPO)): sha(p) for p in source_paths}
    spec = importlib.util.spec_from_file_location("gamma_publication_v2", PUBLICATION)
    publication = importlib.util.module_from_spec(spec); spec.loader.exec_module(publication)
    if args.self_test_only:
        args.output.mkdir(parents=True)
        checks = analytical_checks()
        if any(sha(p) != source_hashes[str(p.relative_to(REPO))] for p in source_paths):
            raise ValueError("Observer or frozen source changed during analytical checks")
        write_json(args.output / "analytical_checks.json", checks)
        write_json(args.output / "self_test.json", {"status": "PASS_ANALYTICAL_CHECKS_ONLY",
            "compute_location": "Unity allocated CPU job", "slurm_job_id": os.environ["SLURM_JOB_ID"],
            "source_hashes": source_hashes, "model_fields_read": False})
        print("REGIONS_ANALYTICAL_PASS", args.output, flush=True)
        return
    requested = [tuple(map(int, pair.split(":"))) for pair in args.pairs]
    if not requested or len(set(requested)) != len(requested) or any(p not in publication.EXPECTED_PAIRS for p in requested):
        raise ValueError("Unknown or duplicate pair")
    if not all(publication.case_ready(args.root, f"{kind}_n{n}_seed{seed}") for n, seed in requested for kind in publication.KINDS):
        raise ValueError("Every requested pair must be complete; no partial analysis")
    args.output.mkdir(parents=True)
    audit = args.output / "native_audit"; audit.mkdir()
    checks = analytical_checks(); write_json(args.output / "analytical_checks.json", checks)
    cases, effects = [], []
    for n, seed in requested:
        pair = [observer_case(publication, args.root, f"{kind}_n{n}_seed{seed}", audit, args.output) for kind in publication.KINDS]
        publication.pair_check(pair); cases.extend(pair)
        for region in CUMULATIVE:
            changes = []
            for case in pair:
                s = case["region_summary"]
                start, end = (s[time][region]["radial_separation_Lgamma"] for time in ("initial", "final"))
                changes.append(end - start if start is not None and end is not None else None)
            effects.append({"grid": n, "seed": seed, "region": region, "polar_change_Lgamma": changes[0],
                "flat_change_Lgamma": changes[1], "polar_minus_flat_change_Lgamma": changes[0] - changes[1] if None not in changes else None})
    plot(cases, args.output)
    if any(sha(p) != source_hashes[str(p.relative_to(REPO))] for p in source_paths):
        raise ValueError("Observer, gate or frozen source changed during analysis")
    failed = {c["name"]: c["summary"]["snapshot_crosschecks"] for c in cases if not c["summary"]["snapshot_crosschecks"]["pv_enstrophy_all_passed"]}
    report = {"status": "COMPLETE_WITH_PV_DIAGNOSTIC_DISCREPANCY" if failed else "COMPLETE", "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "compute_location": "Unity allocated CPU job", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": sha(__file__), "publication_v2_sha256": sha(PUBLICATION), "core_sha256": EXPECTED_CORE_SHA,
        "source_hashes": source_hashes, "runtime": {"numpy": np.__version__, "matplotlib": matplotlib.__version__, "h5py": publication.h5py.__version__},
        "slurm_sha256": sha(REPO / "slurm/gamma2026-regions.sbatch"), "regions": DEFINITIONS,
        "mask_arithmetic": "Compare hypot(x_m,y_m) directly with each physical threshold a*Lgamma, matching the frozen r<5Lgamma gate; normalize only reported radii/centroids.",
        "fraction_denominators": "All-domain relative-enstrophy integral, or all-domain circulation magnitude of the same sign. No mass, regional energy or PV-invariant budget is inferred.",
        "undefined_values": "Null in JSON / empty in CSV when a mask is empty, sign is absent, or denominator is zero; never substituted with zero radius.",
        "requested_pairs": args.pairs, "inventory": publication.inventory(args.root), "cases": [c["region_summary"] for c in cases],
        "paired_changes": effects, "analytical_checks": checks, "vorticity_energy_relative_enstrophy_checks_passed": True,
        "pv_snapshot_crosscheck_passed": not failed, "pv_audit_failures": failed,
        "movies": MOVIES,
        "limitations": ["Fixed Eulerian origin/masks; moments describe signed vorticity, not identified or tracked vortices",
            "The trap is not a material or impermeable boundary; mask contents can change through advection and background-PV interactions",
            "Relative enstrophy is not an invariant for spatially varying background PV; regional changes are not a closed transport/dissipation budget",
            "No convection, mass, physical depths, regional-energy partition or nonlinear cascade flux is computed",
            "PV discrepancies retain the original tolerance and are not presumed to be Float32 roundoff",
            "Prespecified interior masks reduce interpretive ambiguity; they do not remove all trap or periodic-boundary effects"]}
    (args.output / "README.md").write_text("# Fixed-region gamma-plane diagnostics\n\n"
        "Completed two-dimensional, unforced barotropic QG cases; all76 daily fields are included. "
        "These fixed regions describe actual signed vorticity and relative enstrophy, not tracked vortices, convection or inverse energy transfer.\n\n"
        "Masks: r/Lγ <1, 1–3, 3–4, 4–5 and ≥5; cumulative interiors r/Lγ <3, <4 and <5. "
        "Each annulus includes its lower boundary and excludes its upper boundary. The exterior is restricted to the periodic computational square. "
        "Fractions use whole-domain denominators of the same sign/quantity. Undefined moments are blank/null.\n\n"
        "Native energy and relative-enstrophy checks retain the frozen publication thresholds. "
        + ("**Unresolved total-PV diagnostic failures remain**; see regions.json and native_audit. No tight PV-conservation claim is made.\n\n" if failed else "Total-PV saved-field comparisons pass the specified threshold in this selection; this is not conservation of total PV enstrophy.\n\n")
        + "The figure shows every selected case; CSVs contain all regions/times, JSON records initial/final values, matched changes and provenance. "
        "Positive paired change means polar minus flat change, not necessarily positive absolute separation. "
        "Changing contents of fixed masks does not measure individual vortex migration.\n\n"
        f"[Watch the already verified gamma comparison movies]({MOVIES}). "
        "This observer does not generate new movies.\n")
    assets = list(args.output.rglob("*"))
    report["assets"] = {str(p.relative_to(args.output)): {"bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(assets) if p.is_file()}
    write_json(args.output / "regions.json", report)
    print("REGIONS_READY", args.output, flush=True)


if __name__ == "__main__":
    main()
