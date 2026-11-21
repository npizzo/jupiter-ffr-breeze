#!/bin/bash
#SBATCH --job-name=ffr-vortex-stability
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=logs/formed-stability-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/formed_stability_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/formed_stability_v1/driver.jl "$2" "data/formed-stability-v1-$2-$SLURM_JOB_ID" "${3:-none}"
