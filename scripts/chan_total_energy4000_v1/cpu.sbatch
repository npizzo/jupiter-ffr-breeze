#!/bin/bash
#SBATCH --job-name=chan-total-energy4000-v1-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-total-energy4000-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy4000_v1/environment.sh
python scripts/chan_total_energy4000_v1/gate.py preflight "data/chan-total-energy4000-v1-cpu-$SLURM_JOB_ID"
