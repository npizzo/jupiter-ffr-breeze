#!/bin/bash
#SBATCH --job-name=ffr-cai-a100
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=a100
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cai-a100-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export FFR_JOB_START_EPOCH=$(date +%s)
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec tools/point4-analysis-env/bin/python scripts/cai_a100_v1/progress.py "$1" "data/cai-a100-v1-$SLURM_JOB_ID"
