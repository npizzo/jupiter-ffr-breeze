#!/bin/bash
#SBATCH --job-name=marcus-fine-pilot
#SBATCH --partition=uri-gpu,gpu-highvram,gpu-preempt-highvram
#SBATCH --constraint=(a100|h100)&vram80
#SBATCH --exclude=uri-gpu006
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/marcus-fine-pilot-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_extend_v1/environment.sh
python scripts/marcus_cv1_extend_v1/gate.py verify "$1"
out="/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127/marcus-fine-pilot-$SLURM_JOB_ID"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/marcus_cv1_extend_v1/driver.jl finepilot "$out"
python - "$out" "$1" <<'PY2'
import sys,tomllib,subprocess
from pathlib import Path
r=tomllib.loads((Path(sys.argv[1])/'FINE_RELEASE.toml').read_text());assert r['status']=='PASS'
subprocess.run(['python','scripts/marcus_cv1_extend_v1/chain.py','start','fine',sys.argv[2],r['checkpoint'],str(r['restart_time']),str(r['restart_iteration'])],check=True)
PY2
