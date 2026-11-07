#!/bin/bash
#SBATCH --job-name=ffr-pilot
#SBATCH --partition=uri-gpu,gpu
#SBATCH --gpus=1
#SBATCH --constraint=sm_80
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/pilot-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/run.jl configs/smoke.toml data/gpu-smoke GPU
"$JULIA_BIN" --project=. scripts/check_output.jl data/gpu-smoke
"$JULIA_BIN" --project=. scripts/compare_architectures.jl data/cpu-smoke data/gpu-smoke
"$JULIA_BIN" --project=. scripts/run.jl configs/pilot.toml data/pilot GPU
"$JULIA_BIN" --project=. scripts/check_output.jl data/pilot
