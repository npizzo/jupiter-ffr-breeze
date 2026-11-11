"""Convert the completed allocated-node audit's compact metadata to JSON."""
import json
import os
from pathlib import Path
import sys
import tomllib

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Run this publication conversion in the Unity audit allocation")
folder = Path(sys.argv[1])
report = tomllib.loads((folder / "audit.toml").read_text())
assert report["status"] == "COMPLETE" and report["parent_inputs_unchanged"]
assert report["model_steps_executed"] == 0
(folder / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps({"status": report["status"], "case": report["case"],
                  "worst_day": report["worst_probe_row"]["time_days"],
                  "worst_export_solver_relative_pv_difference": report["worst_probe_row"]["snapshot_pv_enstrophy_relative_error"],
                  "checkpoint_model_vs_independent_real_pv_relative_difference": report["checkpoint_state_diagnostics"]["model_vs_independent_real_pv_relative_difference"],
                  "checkpoint_model_vs_independent_spectral_pv_relative_difference": report["checkpoint_state_diagnostics"]["model_vs_independent_spectral_pv_relative_difference"],
                  "parent_inputs_unchanged": report["parent_inputs_unchanged"]}, indent=2))
