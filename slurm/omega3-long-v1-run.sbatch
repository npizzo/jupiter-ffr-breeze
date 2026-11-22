#!/bin/bash
#SBATCH --job-name=ffr-omega3-long-run
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:30:00
#SBATCH --output=logs/omega3-long-v1-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega3_long_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/omega3_long_v1/driver.jl run "data/omega3-long-v1-run-$SLURM_JOB_ID"
