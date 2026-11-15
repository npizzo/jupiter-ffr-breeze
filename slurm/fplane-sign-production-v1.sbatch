#!/bin/bash
#SBATCH --job-name=ffr-sign-sixhour
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=logs/fplane-sign-production-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cases=(plus minus)
"$JULIA_BIN" --project=. scripts/fplane_sign_v1/cuda_launch.jl campaign data/fplane-sign-v1 "${cases[$SLURM_ARRAY_TASK_ID]}"
