#!/bin/bash
#SBATCH --job-name=ffr-eddies-release
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:03:00
#SBATCH --output=logs/chan-eddies-v3-release-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_eddies_v3/gate.py verify "$1"
printf '%s  %s\n' 6fad1b2e63d427244fcc8640705b1c17e7f4186c020ce58861c1cfd15aaf845b work/chan_eddies_v3_screen_review.py | sha256sum --check --status
python work/chan_eddies_v3_screen_review.py "$2" "$1"
