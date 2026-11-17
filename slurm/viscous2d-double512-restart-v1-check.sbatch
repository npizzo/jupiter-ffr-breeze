#!/bin/bash
#SBATCH --job-name=ffr-double-restart-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/viscous2d-double512-restart-check-%j.out
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/viscous2d_double512_restart_v1/check.jl data/viscous2d-double512-v1-1885121/to_200/checkpoint.jld2 data/viscous2d-double512-v1-1885121/to_200/zeta_00064000.f32 "data/viscous2d-double512-restart-check-${SLURM_JOB_ID}"
