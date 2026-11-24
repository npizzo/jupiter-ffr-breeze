#!/bin/bash
#SBATCH --job-name=ffr-weno9-lsv2000-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=logs/weno9-lsv2000-v5-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/weno9_lsv2000_v5/gate.py preflight "data/weno9-lsv2000-v5-cpu-$SLURM_JOB_ID"
python scripts/weno9_lsv2000_v5/release.py "data/weno9-lsv2000-v5-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
