#!/bin/bash
#SBATCH --job-name=ffr-totalE-validation
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=00:20:00
#SBATCH --output=logs/chan-total-energy-v1-gpu-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy_v1/environment.sh
python scripts/chan_total_energy_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/chan_total_energy_v1/driver.jl gpu_validate "data/chan-total-energy-v1-gpu-${SLURM_JOB_ID}"
