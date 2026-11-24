#!/bin/bash
#SBATCH --job-name=ffr-budget-resolution
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --array=0-2%3
#SBATCH --output=logs/chan-budget-resolution-v2-run-%A_%a.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_budget_resolution_v2/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/chan_budget_resolution_v2/driver.jl run "data/chan-budget-resolution-v2-run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" "$SLURM_ARRAY_TASK_ID" "$(dirname "$1")/prepared"
