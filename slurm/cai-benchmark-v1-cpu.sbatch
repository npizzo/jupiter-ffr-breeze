#!/bin/bash
#SBATCH --job-name=ffr-cai-cpu-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:08:00
#SBATCH --output=logs/cai-benchmark-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1
exec "$JULIA_BIN" --project=. scripts/cai_benchmark_v1/cpu_check.jl "data/cai-benchmark-v1-cpu-$SLURM_JOB_ID"
