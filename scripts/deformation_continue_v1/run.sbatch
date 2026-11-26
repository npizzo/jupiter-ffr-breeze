#!/bin/bash
#SBATCH --job-name=ffr-deformation100
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deformation-continue-v1-run-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/deformation_continue_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/deformation_continue_v1/driver.jl run "data/deformation-continue-v1-$3-$SLURM_JOB_ID" "$2"
