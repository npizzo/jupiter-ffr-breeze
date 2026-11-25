#!/bin/bash
#SBATCH --job-name=chan-theta-diagnostic4-gate
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/chan-theta-diagnostic4000-v1-release-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
tools/point4-analysis-env/bin/python scripts/chan_theta_diagnostic4000_v1/health.py check "runs/chan-theta-diagnostic4000-v1/HEALTH_PASS.json"
python scripts/chan_theta_diagnostic4000_v1/release.py "$1" "runs/chan-theta-diagnostic4000-v1/HEALTH_PASS.json"
