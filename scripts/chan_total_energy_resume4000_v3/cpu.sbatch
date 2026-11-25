#!/bin/bash
#SBATCH --job-name=chan-total-energy-resume4000-v3-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-total-energy-resume4000-v3-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy_resume4000_v3/environment.sh
python scripts/chan_total_energy_resume4000_v3/gate.py preflight "data/chan-total-energy-resume4000-v3-cpu-$SLURM_JOB_ID"

tools/point4-analysis-env/bin/python scripts/chan_total_energy_resume4000_v3/health.py check runs/chan-total-energy-resume4000-v3/HEALTH_PASS.json
python scripts/chan_total_energy_resume4000_v3/release.py "data/chan-total-energy-resume4000-v3-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json" runs/chan-total-energy-resume4000-v3/HEALTH_PASS.json
