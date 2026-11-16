#!/bin/bash
#SBATCH --job-name=ffr-2d-long-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/instability2d-long-v1-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/instability2d_long_v1/check.jl "data/instability2d-long-v1-check-${SLURM_JOB_ID}" data/instability2d-light-v1-1810727/checkpoint.jld2
