#!/bin/bash
#SBATCH --job-name=ffr-fixed45-run
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/chan-fixed45-v2-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_fixed45_v2/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/chan_fixed45_v2/driver.jl run "data/chan-fixed45-v2-run-$SLURM_JOB_ID"
