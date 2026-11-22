#!/bin/bash
#SBATCH --job-name=ffr-response-ramp
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:40:00
#SBATCH --output=logs/chan-rotation-response-v1-ramp-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_rotation_response_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/chan_rotation_response_v1/driver.jl ramp "data/chan-rotation-response-v1-ramp-$SLURM_JOB_ID" "$2"
