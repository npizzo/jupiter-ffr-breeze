#!/bin/bash
#SBATCH --job-name=ffr-chan-pair
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=a100
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/chan2013-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/chan2013_v1/launch.jl "$1" "data/chan2013-v1-$SLURM_JOB_ID" "$2"
