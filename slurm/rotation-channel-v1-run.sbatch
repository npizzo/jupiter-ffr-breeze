#!/bin/bash
#SBATCH --job-name=ffr-rotation-channel
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=07:30:00
#SBATCH --output=logs/rotation-channel-v1-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/rotation_channel_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/rotation_channel_v1/driver.jl run "data/rotation-channel-v1-$2-$SLURM_JOB_ID" "$2"
