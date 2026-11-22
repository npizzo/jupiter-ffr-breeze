#!/bin/bash
#SBATCH --job-name=ffr-omega15-long-run
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/omega15-long-v1-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega15_long_v1/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/omega15_long_v1/driver.jl run "data/omega15-long-v1-run-$SLURM_JOB_ID"
