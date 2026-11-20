#!/bin/bash
#SBATCH --job-name=ffr-chan-validate
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=a100
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan2013-v1-validate-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/chan2013_v1/validate.jl "data/chan2013-v1-validation-$SLURM_JOB_ID"
