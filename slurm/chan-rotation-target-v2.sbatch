#!/bin/bash
#SBATCH --job-name=ffr-rotation-target
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/chan-rotation-target-v2-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_rotation_target_v2/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/chan_rotation_target_v2/driver.jl target "data/chan-rotation-target-v2-$SLURM_JOB_ID"
