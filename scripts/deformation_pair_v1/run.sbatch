#!/bin/bash
#SBATCH --job-name=ffr-deformation-pair
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/deformation-pair-v1-run-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/deformation_pair_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/deformation_pair_v1/driver.jl run "data/deformation-pair-v1-$3-$SLURM_JOB_ID" "$2"
