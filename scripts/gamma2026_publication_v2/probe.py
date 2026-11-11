#!/usr/bin/env python3
"""Scan completed native snapshots without rendering or altering a simulation."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", default="flat_n1024_seed11")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("The probe requires an allocated Unity job")
    if args.output.exists():
        raise ValueError("Use a fresh probe directory")
    path = Path(__file__).with_name("publish.py")
    spec = importlib.util.spec_from_file_location("gamma_publication_v2", path)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    args.output.mkdir(parents=True)
    result = driver.case_analysis(args.root, args.case, args.output)
    checks = result["summary"]["snapshot_crosschecks"]
    report = {
        "status": "PASS" if checks["all_passed"] else "PV_DIAGNOSTIC_DISCREPANCY",
        "scope": "Native-grid reader and per-field diagnostic consistency; no integration or rendering",
        "case": result["name"],
        "compute_location": "Unity allocated CPU job",
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "saved_frames_checked": len(result["files"]),
        "snapshot_crosschecks": checks,
        "renderer_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "probe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "movie_rendering_tested": False,
        "note": "Completed probe execution does not imply all scientific checks passed. Original PV tolerance is unchanged and failed days are explicit. The generated numeric diagnostic CSV can be passed to the native checkpoint audit."
    }
    (args.output / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
