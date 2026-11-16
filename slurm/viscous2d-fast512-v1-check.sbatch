#!/bin/bash
#SBATCH --job-name=ffr-512-fast-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/viscous2d-fast512-check-%j.out
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/viscous2d_fast512_v1/check.jl data/viscous2d-long512-v1-1848831 "data/viscous2d-fast512-check-${SLURM_JOB_ID}"
