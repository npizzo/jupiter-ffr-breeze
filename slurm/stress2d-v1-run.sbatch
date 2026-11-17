#!/bin/bash
#SBATCH --job-name=ffr-stress-run
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/stress2d-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
CHECK=$1
CASE=$2
OUT=$3
TARGET=${4:-400}
"$JULIA_BIN" --project=. scripts/stress2d_v1/run.jl "$OUT" "$CHECK/initial.jld2" "$TARGET" "$CHECK/PASS.toml" "$CASE"
