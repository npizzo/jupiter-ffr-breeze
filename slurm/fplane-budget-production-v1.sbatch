#!/bin/bash
#SBATCH --job-name=ffr-dense-vorticity
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram80
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=logs/fplane-budget-production-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cases=(mature512 mature1024)
"$JULIA_BIN" --project=. scripts/fplane_budget_v1/cuda_launch.jl campaign data/fplane-budget-v1 "${cases[$SLURM_ARRAY_TASK_ID]}"
