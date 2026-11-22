#!/bin/bash
#SBATCH --job-name=ffr-omega3-anelastic-run
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=logs/omega3-anelastic-v3-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega3_anelastic_v3/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/omega3_anelastic_v3/driver.jl run "data/omega3-anelastic-v3-run-$SLURM_JOB_ID"
