#!/bin/bash
#SBATCH --job-name=ffr-ferr26-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/ferrel26-restart-audit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel2026_restart_audit/audit.jl "$1" "$2"
