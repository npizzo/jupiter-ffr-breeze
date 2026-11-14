#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-width
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ferrel35-width-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/ferrel35_width_v3/cuda_launch.jl scripts/ferrel35_width_v3/campaign.jl "$1" "$2" "$3" "$4"
tools/field-export-env/bin/python scripts/ferrel35_width_v3/transfer.py "$3/$1/segment_$2"
