#!/bin/bash
#SBATCH --job-name=ffr-strong256-resume
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint="a100|h100"
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05:50:00
#SBATCH --output=logs/strong256-resume-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export FFR_JOB_START_EPOCH=$(date +%s)
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/strong256_resume_v1/launch.jl "data/strong256-resume-v1-$SLURM_JOB_ID" "$1" "$2" "$3"
