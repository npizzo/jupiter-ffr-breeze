#!/bin/bash
#SBATCH --job-name=ffr-jet-linear-four
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/beta-shear-linear-v1-gpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_linear_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/beta_shear_linear_v1/driver.jl run "data/beta-shear-linear-v1-$SLURM_JOB_ID"
