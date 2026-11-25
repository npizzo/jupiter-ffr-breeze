#!/bin/bash
#SBATCH --job-name=chan-clean-anelastic60-gate
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/chan-clean-anelastic6000-v1-release-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
tools/point4-analysis-env/bin/python scripts/chan_clean_anelastic6000_v1/health.py check "runs/chan-clean-anelastic6000-v1/HEALTH_PASS.json"
python scripts/chan_clean_anelastic6000_v1/release.py "$1" "runs/chan-clean-anelastic6000-v1/HEALTH_PASS.json"
