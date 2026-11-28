#!/bin/bash
#SBATCH --job-name=marcus-CV1-four-days
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=(a100|h100)&vram80
#SBATCH --exclude=uri-gpu006
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/marcus-fourday-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_fourday_v1/environment.sh
python scripts/marcus_cv1_fourday_v1/gate.py verify "$1"
python - "$2" <<'PY'
import sys,tomllib
from pathlib import Path
r=tomllib.loads((Path(sys.argv[1])/'GPU_PASS.toml').read_text())
assert r['status']=='PASS' and r['pilot_final_time']>=1800 and r['estimated_four_day_wall_hours']<12
print('RESOLVED_PILOT_RELEASE_PASS',flush=True)
PY
out="/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127/marcus-CV1-fourday-$SLURM_JOB_ID"
printf '%s\n' "$out" > "runs/marcus-cv1-fourday-v1/job-$SLURM_JOB_ID-output.txt"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/marcus_cv1_fourday_v1/driver.jl run "$out"
