#!/bin/bash
#SBATCH --job-name=ffr-strong48
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --constraint=sm_80
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/strong48-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/run_strong_campaign.jl "${1:-data/strong48-v1}" GPU
