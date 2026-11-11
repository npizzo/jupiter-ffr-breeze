#!/usr/bin/env python3
"""Exercise the frozen publication reader and native budget checks on Unity."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", default="polar_n1024_seed11")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("The schema/budget probe requires an allocated Unity job")
    if args.output.exists():
        raise ValueError("Probe output already exists; choose a fresh directory")
    driver_path = Path(__file__).with_name("publish.py")
    driver_sha = hashlib.sha256(driver_path.read_bytes()).hexdigest()
    expected_sha = "915aae9946ed5bc0c7aa0aa3c53724b9f5a1c03bf1e8b0f31b5362ab224eb290"
    if driver_sha != expected_sha:
        raise ValueError("The probe must exercise the reviewed frozen publication driver")
    spec = importlib.util.spec_from_file_location("gamma_publication_driver", driver_path)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    if not driver.case_ready(args.root, args.case):
        raise ValueError("Probe case is not complete")
    args.output.mkdir(parents=True)
    print("Checking native JLD2 schema and every saved-frame budget:", args.case, flush=True)
    result = driver.case_analysis(args.root, args.case, args.output)
    summary = result["summary"]
    report = {
        "status": "PASS",
        "scope": "Publication reader, completed-case guards, native-grid Fourier energy, relative/PV enstrophy and diagnostic output schema",
        "compute_location": "Unity allocated CPU job",
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "case": result["name"],
        "saved_frames_checked": len(result["files"]),
        "completed_days": summary["completed_days"],
        "snapshot_crosschecks": summary["snapshot_crosschecks"],
        "frozen_driver_sha256": driver_sha,
        "probe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "movie_rendering_tested": False,
        "note": "Validation only. No simulation source changed and no pair movie or scientific interpretation is published by this probe."
    }
    (args.output / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
