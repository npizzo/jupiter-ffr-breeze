#!/bin/bash
#SBATCH --job-name=chan-total-energy4000-v2
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/chan-total-energy4000-v2-run-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy4000_v2/environment.sh
python scripts/chan_total_energy4000_v2/gate.py verify "$1"
tools/point4-analysis-env/bin/python scripts/chan_total_energy4000_v2/health.py verify "$2"
export HEALTH_TOML="${2%.json}.toml"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/chan_total_energy4000_v2/driver.jl run "data/chan-total-energy4000-v2-run-${SLURM_JOB_ID}"
