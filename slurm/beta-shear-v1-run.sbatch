#!/bin/bash
#SBATCH --job-name=ffr-beta-shear-run
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/beta-shear-v1-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_v1/gate.py verify-stability "$1" "$3"
"$JULIA_BIN" --project=. scripts/beta_shear_v1/driver.jl run "data/beta-shear-v1-$2-$SLURM_JOB_ID" "$2"
