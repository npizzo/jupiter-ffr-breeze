#!/bin/bash
#SBATCH --job-name=ffr-weno9-startup-pilot
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --array=0-3
#SBATCH --output=logs/weno9-startup-v2-pilot-%A_%a.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/weno9_startup_v2/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/weno9_startup_v2/driver.jl pilot "data/weno9-startup-v2-pilot-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" "$SLURM_ARRAY_TASK_ID"
