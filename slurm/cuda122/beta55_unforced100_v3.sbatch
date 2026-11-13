#!/bin/bash
#SBATCH --job-name=ffr-beta55_unforced100
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/beta55_unforced100_v3-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/beta55_unforced100_v3/cuda_launch.jl scripts/beta55_unforced100_v3/campaign.jl "$1" "$2" "$3"
tools/field-export-env/bin/python scripts/beta55_unforced100_v3/transfer.py "$2/beta55_unforced100/segment_$1"
