#!/bin/bash
#SBATCH --job-name=ffr-convection48
#SBATCH --partition=uri-gpu,gpu
#SBATCH --gpus=1
#SBATCH --constraint=sm_80
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/convection48-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
nvidia-smi
"$JULIA_BIN" --project=. scripts/run_convection_campaign.jl data/convection48-v1 GPU
