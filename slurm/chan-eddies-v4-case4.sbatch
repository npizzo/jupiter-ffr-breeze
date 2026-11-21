#!/bin/bash
#SBATCH --job-name=ffr-eddies-case4
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --output=logs/chan-eddies-v4-case4-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_eddies_v4/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/chan_eddies_v4/driver.jl case4 "data/chan-eddies-v4-case4-$SLURM_JOB_ID" "$2"
