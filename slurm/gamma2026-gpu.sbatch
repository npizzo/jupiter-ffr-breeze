#!/bin/bash
#SBATCH --job-name=ffr-gamma-gate
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/gamma2026-gpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-gamma2026"
export JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=environments/gamma2026 scripts/gamma2026/validate.jl gpu configs/gamma2026/polar_seed11.toml "$1/preflight"
