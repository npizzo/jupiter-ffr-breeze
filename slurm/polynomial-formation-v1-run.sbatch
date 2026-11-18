#!/bin/bash
#SBATCH --job-name=ffr-poly-noadv
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:36:00
#SBATCH --output=logs/polynomial-formation-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --project=. scripts/polynomial_formation_v1/run.jl "$1" "$2" "$3" "$4" no_advection
