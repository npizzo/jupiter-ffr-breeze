#!/bin/bash
#SBATCH --job-name=ffr-ferr26v2-unit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/ferrel26v2-unit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel2026_v2/check_initial.jl "$1"
