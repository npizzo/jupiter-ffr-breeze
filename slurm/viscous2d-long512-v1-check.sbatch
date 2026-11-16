#!/bin/bash
#SBATCH --job-name=ffr-512-restart-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/viscous2d-long512-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/viscous2d_long512_v1/check.jl data/viscous2d-refine-v1-1848240/initial.jld2 "data/viscous2d-long512-check-${SLURM_JOB_ID}"
