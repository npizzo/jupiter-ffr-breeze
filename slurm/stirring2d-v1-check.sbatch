#!/bin/bash
#SBATCH --job-name=ffr-2d-stir-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/stirring2d-v1-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --startup-file=no -e 'for p in ARGS; Meta.parseall(read(p,String)); end; println("SOURCE_PARSE_PASS")' scripts/stirring2d_v1/{forcing,run,check}.jl
"$JULIA_BIN" --project=. scripts/stirring2d_v1/check.jl "data/stirring2d-v1-check-${SLURM_JOB_ID}" data/instability2d-light-v1-1810727/checkpoint.jld2
