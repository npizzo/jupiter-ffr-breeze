#!/bin/bash
#SBATCH --job-name=ffr-cascade35
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/cascade35-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/cascade35/campaign.jl "$1" "$2" "$3"
tools/field-export-env/bin/python scripts/cascade35/transfer.py "$2/cascade35_heated/segment_$1"
