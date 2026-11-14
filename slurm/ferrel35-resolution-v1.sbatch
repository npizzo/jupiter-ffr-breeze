#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-resolution
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel35-resolution-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --project=. scripts/ferrel35_resolution_v1/cuda_launch.jl "$@"
