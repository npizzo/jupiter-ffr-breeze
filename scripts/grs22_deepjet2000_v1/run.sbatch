#!/bin/bash
#SBATCH --job-name=ffr-grs22-deepjet2000-long
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/grs22-deepjet2000-v1-run-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/grs22_deepjet2000_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/grs22_deepjet2000_v1/driver.jl run "data/grs22-deepjet2000-v1-run-${SLURM_JOB_ID}"
