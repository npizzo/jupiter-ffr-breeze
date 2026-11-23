#!/bin/bash
#SBATCH --job-name=ffr-weno9-lsv100-long
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/weno9-lsv100-v3-run-%A_%a.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/weno9_lsv100_v3/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/weno9_lsv100_v3/driver.jl run "data/weno9-lsv100-v3-run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" "$SLURM_ARRAY_TASK_ID" "data/weno9-startup-v2-pilot-${2}_${SLURM_ARRAY_TASK_ID}/parent"
