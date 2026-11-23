#!/bin/bash
#SBATCH --job-name=ffr-shear-physics-pair
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/beta-shear-physics-v1-gpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_physics_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/beta_shear_physics_v1/driver.jl run "data/beta-shear-physics-v1-run-$SLURM_JOB_ID" "${2:-complete}" "$1"
