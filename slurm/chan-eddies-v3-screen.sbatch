#!/bin/bash
#SBATCH --job-name=ffr-eddies-screen
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --output=logs/chan-eddies-v3-screen-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_eddies_v3/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/chan_eddies_v3/driver.jl screen "data/chan-eddies-v3-screen-$SLURM_JOB_ID"
