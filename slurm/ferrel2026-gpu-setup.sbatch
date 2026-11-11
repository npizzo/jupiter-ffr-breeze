#!/bin/bash
#SBATCH --job-name=ffr-ferrel-runtime
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ferrel2026-gpu-setup-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-ferrel2026-gpu:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
nvidia-smi --query-gpu=name,driver_version,uuid --format=csv,noheader
"$JULIA_BIN" --startup-file=no --project=. scripts/ferrel2026_gpu/setup.jl "$1"
