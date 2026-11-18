#!/bin/bash
#SBATCH --job-name=ffr-strong-seed-50
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=l40s
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/strong-seed-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/strong_seed_v1/launch.jl "data/strong-seed-v1-$SLURM_JOB_ID"
