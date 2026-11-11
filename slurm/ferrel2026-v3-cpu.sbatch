#!/bin/bash
#SBATCH --job-name=ffr-ferr26v3-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ferrel26v3-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=8
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel2026_v3/cpu_smoke.jl "$1/cpu_smoke" "$1"
