#!/bin/bash
#SBATCH --job-name=ffr-repair-diagnostic
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=logs/beta-shear-profile-v1-gpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_profile_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/beta_shear_profile_v1/driver.jl "$2" "data/beta-shear-profile-v1-$2-$SLURM_JOB_ID"
