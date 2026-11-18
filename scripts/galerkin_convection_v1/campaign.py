"""Run exactly the declared five cases after a successful validation gate."""
import hashlib
import json
from pathlib import Path
import sys
import time

from run import run


root = Path(sys.argv[1])
assert (root/"validation/PASS.json").is_file(), "No passing release check"
manifest = json.loads((root/"source-manifest.json").read_text())
source = Path(__file__).parent
for name, checksum in manifest.items():
    assert hashlib.sha256((source/name).read_bytes()).hexdigest() == checksum, name
config = json.loads((root/"config.json").read_text())
start = time.monotonic()
for case in config["cases"]:
    elapsed = time.monotonic()-start
    if elapsed >= config["cpu_budget_seconds"]:
        raise RuntimeError("Campaign cost limit reached; no further cases released")
    run(root/"config.json", case["name"], root/"runs")
summary = {"status": "COMPLETE", "cases": [c["name"] for c in config["cases"]],
           "elapsed_seconds": time.monotonic()-start}
(root/"CAMPAIGN_COMPLETE.json").write_text(json.dumps(summary, indent=2)+"\n")
print(json.dumps(summary), flush=True)
