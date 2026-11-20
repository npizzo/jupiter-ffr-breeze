#!/bin/bash
#SBATCH --job-name=ffr-anelastic-parent
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/anelastic-parent-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/anelastic_parent_v1/pilot.jl "data/anelastic-parent-v1-$SLURM_JOB_ID" "$1"
