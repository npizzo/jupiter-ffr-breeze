#!/bin/bash
#SBATCH --job-name=ffr-response
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=logs/anelastic-response-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/anelastic_response_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/anelastic_response_v1/driver.jl "$2" "data/anelastic-response-v1-$2-$SLURM_JOB_ID" "${3:-none}"
