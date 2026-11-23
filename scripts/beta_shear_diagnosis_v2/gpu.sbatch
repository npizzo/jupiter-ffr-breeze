#!/bin/bash
#SBATCH --job-name=ffr-jet-diag-four
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=logs/beta-shear-diag-v2-gpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_diagnosis_v2/gate.py verify "$1"
"$JULIA_BIN" --project=. scripts/beta_shear_diagnosis_v2/driver.jl run "data/beta-shear-diag-v2-$SLURM_JOB_ID"
