#!/bin/bash
#SBATCH --job-name=ffr-unforced55-check
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/beta55-unforced-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/cuda_repair/launch.jl scripts/beta55_unforced100/preflight.jl "$1"
