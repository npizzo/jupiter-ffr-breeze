#!/bin/bash
#SBATCH --job-name=ffr-convection120
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --constraint=sm_80
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --array=1-4%1
#SBATCH --output=logs/convection120-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/run_convection120.jl "$SLURM_ARRAY_TASK_ID" "${1:-data/convection120-v1}" GPU
