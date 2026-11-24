#!/bin/bash
#SBATCH --job-name=ffr-budget-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-budget-resolution-v3-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_budget_resolution_v3/gate.py preflight "data/chan-budget-resolution-v3-cpu-$SLURM_JOB_ID"
python scripts/chan_budget_resolution_v3/release.py "data/chan-budget-resolution-v3-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
