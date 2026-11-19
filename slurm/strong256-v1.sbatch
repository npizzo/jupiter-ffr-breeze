#!/bin/bash
#SBATCH --job-name=ffr-strong256-50
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint="a100|h100|l40s"
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/strong256-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/strong256_v1/launch.jl "data/strong256-v1-$SLURM_JOB_ID" "$1"
