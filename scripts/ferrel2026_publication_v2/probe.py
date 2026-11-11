"""Reader/renderer QA on one completed six-hour case, not campaign validation.

The parent task observed the individual case's post-run checks pass while the
multi-case pilot remains active. This observer does not create any science gate.
"""
import os
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("Run the schema probe in an allocated Unity job")

driver = Path(__file__).with_name("render_report.py")
expected_sha = "ad5fc63bf40dcccdcc7ce7ac88608f3e4545f52bdb043d35a25add181b407697"
actual_sha = hashlib.sha256(driver.read_bytes()).hexdigest()
if actual_sha != expected_sha:
    raise SystemExit("Frozen publication driver changed; refusing schema probe")

# Import only after the allocation and immutable-source guards.
from render_report import collect, movie, figures, digest
import numpy as np

source = Path("data/ferrel2026-preflight-gpu-v1/pilot512_b0")
out = Path("rendered/ferrel2026-schema-probe-v2")
if out.exists() and any(out.iterdir()):
    raise SystemExit("Use a new probe revision; refusing to overwrite probe assets")
out.mkdir(parents=True, exist_ok=True)
case = collect("b0", source, out)
m = case["m"]
if [m[k] for k in ("nx", "ny", "nz")] != [512, 512, 96] or m["stop_seconds"] != 21600:
    raise ValueError("Probe expected the completed 512²×96 six-hour unheated pilot")
limits = np.maximum(case["limits"], 1e-12)
media = movie([case], out / "b0_schema_vorticity_tracer.mp4", limits, case["tracer_range"], fps=6)
figures([case], out)
report = {
    "status": "SCHEMA_AND_ENCODING_QA_PASS",
    "scope": "Completed individual-case reader and renderer QA only; no full pilot PASS or physical interpretation",
    "source": str(source), "compute_location": "Unity allocated CPU job",
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "frozen_renderer_sha256": actual_sha,
    "probe_sha256": digest(__file__),
    "case": case["summary"], "media": [media],
    "parent_observation": "Individual pilot512_b0 post-run checks passed within multi-case job1409341; entire campaign still requires its own gate",
    "raw_fields_downloaded": False,
}
(out / "probe-report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
assets = [{"file": p.name, "bytes": p.stat().st_size, "sha256": digest(p)}
          for p in sorted(out.iterdir()) if p.is_file()]
(out / "probe-manifest.json").write_text(json.dumps({
    "status": "READY_FOR_VISUAL_QA", "scope": "reader/renderer QA only", "assets": assets}, indent=2) + "\n")
print("SCHEMA_AND_ENCODING_QA_PASS", out, flush=True)
