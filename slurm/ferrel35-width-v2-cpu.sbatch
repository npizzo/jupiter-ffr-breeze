#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-width-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ferrel35-width-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel35_width_v2/preflight.jl "$1/cpu" CPU
