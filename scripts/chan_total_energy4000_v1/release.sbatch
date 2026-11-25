#!/bin/bash
#SBATCH --job-name=chan-total-energy4000--gate
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/chan-total-energy4000-v1-release-%j.out
set -euo pipefail
source ./local.env
source scripts/chan_total_energy4000_v1/environment.sh
tools/point4-analysis-env/bin/python scripts/chan_total_energy4000_v1/health.py check "runs/chan-total-energy4000-v1/HEALTH_PASS.json"
python scripts/chan_total_energy4000_v1/release.py "$1" "runs/chan-total-energy4000-v1/HEALTH_PASS.json"
