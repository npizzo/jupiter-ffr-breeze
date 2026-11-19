#!/bin/bash
#SBATCH --job-name=ffr-cai-validate
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint="a100|h100"
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=logs/cai-benchmark-v1-validate-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "$JULIA_BIN" --project=. scripts/cai_benchmark_v1/validate.jl "data/cai-benchmark-v1-validation-$SLURM_JOB_ID"
