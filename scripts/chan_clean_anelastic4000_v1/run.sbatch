#!/bin/bash
#SBATCH --job-name=ffr-chan-anel4000
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=logs/chan-clean-anelastic4000-v1-run-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_clean_anelastic4000_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/chan_clean_anelastic4000_v1/driver.jl run "data/chan-clean-anelastic4000-v1-run-${SLURM_JOB_ID}"
