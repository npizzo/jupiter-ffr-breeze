#!/bin/bash
#SBATCH --job-name=ffr-gammav3-pilot
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/gamma2026_v3-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-gamma2026-gpu:$FFR_ROOT/depot-gamma2026"
export JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
case_name="$2"
grid="${3:?grid}"
"$JULIA_BIN" --project=environments/gamma2026_gpu scripts/gamma2026_v3/run.jl "configs/gamma2026_v2/$case_name.toml" "$1/$case_name" "$1/preflight/gpu-$grid.toml"
