#!/bin/bash
#SBATCH --job-name=ffr-totalE-long
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/chan-total-energy-long-v1-run-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy_long_v1/environment.sh
python scripts/chan_total_energy_long_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/chan_total_energy_long_v1/driver.jl run "data/chan-total-energy-long-v1-run-${SLURM_JOB_ID}"
