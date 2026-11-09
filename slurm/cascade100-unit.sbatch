#!/bin/bash
#SBATCH --job-name=ffr-cascade-unit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/cascade-unit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/cascade100/check_initial.jl "$1"
