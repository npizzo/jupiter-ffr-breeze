#!/bin/bash
#SBATCH --job-name=ffr-2d-light
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:50:00
#SBATCH --output=logs/instability2d-light-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/instability2d_light_v1/check.jl "data/instability2d-light-v1-check-${SLURM_JOB_ID}"
test -f "data/instability2d-light-v1-check-${SLURM_JOB_ID}/PASS.toml"
"$JULIA_BIN" --project=. scripts/instability2d_light_v1/run.jl "data/instability2d-light-v1-${SLURM_JOB_ID}"
