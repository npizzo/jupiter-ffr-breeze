#!/bin/bash
#SBATCH --job-name=chan-clean-anelastic6000-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-clean-anelastic10000-v2-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_clean_anelastic10000_v2/gate.py preflight "data/chan-clean-anelastic10000-v2-cpu-$SLURM_JOB_ID"

tools/point4-analysis-env/bin/python scripts/chan_clean_anelastic10000_v2/health.py check runs/chan-clean-anelastic10000-v2/HEALTH_PASS.json
python scripts/chan_clean_anelastic10000_v2/release.py "data/chan-clean-anelastic10000-v2-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json" runs/chan-clean-anelastic10000-v2/HEALTH_PASS.json
