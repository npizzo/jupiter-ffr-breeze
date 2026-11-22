#!/bin/bash
#SBATCH --job-name=ffr-beta-shear-pilot
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/beta-shear-v1-pilot-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/beta_shear_v1/driver.jl run "data/beta-shear-v1-pilot-$SLURM_JOB_ID" pilot
python scripts/beta_shear_v1/gate.py release "$1" "data/beta-shear-v1-pilot-$SLURM_JOB_ID"
