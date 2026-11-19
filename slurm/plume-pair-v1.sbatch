#!/bin/bash
#SBATCH --job-name=ffr-plume-pair
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=a100
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/plume-pair-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/plume_pair_v1/run.jl "$1" "data/plume-pair-v1-$SLURM_JOB_ID"
