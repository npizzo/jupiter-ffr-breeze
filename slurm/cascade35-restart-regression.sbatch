#!/bin/bash
#SBATCH --job-name=ffr-beta35-restart-regression
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/cascade35-restart-regression-%j.out
set -euo pipefail
source local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-cascade-cuda122:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/cascade35_restart_repair/regression.jl
