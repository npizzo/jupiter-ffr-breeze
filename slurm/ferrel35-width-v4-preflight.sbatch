#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-width-check
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ferrel35-width-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel35_width_v4/cuda_launch.jl scripts/ferrel35_width_v4/preflight.jl "$1"
