#!/bin/bash
#SBATCH --job-name=ffr-cuda122-probe
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/cuda122-probe-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
nvidia-smi --query-gpu=name,driver_version --format=csv
"$JULIA_BIN" --project=. scripts/cuda_repair/probe.jl "$1"
