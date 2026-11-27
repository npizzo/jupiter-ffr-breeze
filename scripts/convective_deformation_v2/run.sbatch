#!/bin/bash
#SBATCH --job-name=ffr-convective-deformation
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/convective-deformation-v2-run-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/convective_deformation_v2/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/convective_deformation_v2/driver.jl run "data/convective-deformation-v2-$3-$SLURM_JOB_ID" "$2"
