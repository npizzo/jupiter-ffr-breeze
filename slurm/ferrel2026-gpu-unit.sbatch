#!/bin/bash
#SBATCH --job-name=ffr-ferrgpu-unit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/ferrel2026-gpu-unit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-ferrel2026-gpu:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --startup-file=no --project=environments/ferrel2026_gpu scripts/ferrel2026_gpu/launch.jl "$1" CPU scripts/ferrel2026_v3/check_initial.jl "$1"
