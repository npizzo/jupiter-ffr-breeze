"""Prepare protocol and execute independent release validation on Unity."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from run import prepare


root = Path(sys.argv[1])
prepare(root)
config = json.loads((root/"config.json").read_text())
source = Path(__file__).parent
command = [sys.executable, str(source/"validate.py"), "--output", str(root/"validation"),
           "--K", "4", "--Nv", "2", "--refined-K", "6", "--refined-Nv", "3",
           "--kfund", str(config["kfund"]), "--nu", str(config["nu"]),
           "--kappa", str(config["kappa"]), "--G", str(config["G"]),
           "--f", "1", "--dt", str(config["dt"]), "--duration", str(config["duration"]),
           "--base-runs", "4", "--refined-runs", "1", "--cpu-budget-seconds", "1800"]
subprocess.run(command, check=True)
assert (root/"validation/PASS.json").is_file()
manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob("*.py")}
(root/"source-manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
print("RELEASE_VALIDATION_PASS", flush=True)
