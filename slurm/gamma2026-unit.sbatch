#!/bin/bash
#SBATCH --job-name=ffr-gamma-unit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/gamma2026-unit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-gamma2026"
export JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=environments/gamma2026 scripts/gamma2026/validate.jl unit configs/gamma2026/polar_seed11.toml "$1/preflight"
