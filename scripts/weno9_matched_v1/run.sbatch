#!/bin/bash
#SBATCH --job-name=ffr-weno9-matched-long
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/weno9-matched-v1-run-%A_%a.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/weno9_matched_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/weno9_matched_v1/driver.jl run "data/weno9-matched-v1-run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" "$SLURM_ARRAY_TASK_ID" "data/weno9-matched-v1-pilot-${2}_${SLURM_ARRAY_TASK_ID}/parent"
