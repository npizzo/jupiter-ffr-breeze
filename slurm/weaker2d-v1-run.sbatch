#!/bin/bash
#SBATCH --job-name=ffr-2d-weaker
#SBATCH --partition=uri-cpu,cpu
#SBATCH --array=0-2%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/weaker2d-v1-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
test -f "$1/PASS.json"
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cases=(forced random_sign cyclonic_sign)
case_name=${cases[$SLURM_ARRAY_TASK_ID]}
"$JULIA_BIN" --project=. scripts/weaker2d_v1/run.jl "$2/$case_name" "$case_name" data/instability2d-light-v1-1810727/checkpoint.jld2
