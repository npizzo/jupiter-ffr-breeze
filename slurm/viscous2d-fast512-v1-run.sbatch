#!/bin/bash
#SBATCH --job-name=ffr-512-fast
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/viscous2d-fast512-%j.out
#SBATCH --no-requeue
set -euo pipefail
test -f "$1/PASS.toml"
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/viscous2d_fast512_v1/run.jl "$2" "$3" "$4" 400 "$1/PASS.toml"
