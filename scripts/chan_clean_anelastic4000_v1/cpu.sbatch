#!/bin/bash
#SBATCH --job-name=ffr-chan-anel4000-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-clean-anelastic4000-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_clean_anelastic4000_v1/gate.py preflight "data/chan-clean-anelastic4000-v1-cpu-$SLURM_JOB_ID"
python scripts/chan_clean_anelastic4000_v1/release.py "data/chan-clean-anelastic4000-v1-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
