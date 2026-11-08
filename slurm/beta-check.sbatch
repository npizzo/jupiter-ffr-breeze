#!/bin/bash
#SBATCH --job-name=ffr-beta-check
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --constraint=sm_80
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/beta-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/check_beta_campaign.jl "${1:-data/regional-check-v1}"
