#!/bin/bash
#SBATCH --job-name=ffr-chan-clean-pair
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --array=0-1%2
#SBATCH --output=logs/chan-clean-pair-v1-run-%A_%a.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_clean_pair_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/chan_clean_pair_v1/driver.jl run "data/chan-clean-pair-v1-run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" "$SLURM_ARRAY_TASK_ID"
