#!/bin/bash
#SBATCH --job-name=ffr-large40
#SBATCH --partition=uri-gpu
#SBATCH --gpus=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/large40-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/run_large_campaign.jl "$1" "${2:-data/large40-v1}" "${3:-data/large40-preflight-v1}"
