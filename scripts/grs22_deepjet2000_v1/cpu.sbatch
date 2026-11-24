#!/bin/bash
#SBATCH --job-name=ffr-grs22-deepjet2000-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=logs/grs22-deepjet2000-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/grs22_deepjet2000_v1/gate.py preflight "data/grs22-deepjet2000-v1-cpu-$SLURM_JOB_ID"
python scripts/grs22_deepjet2000_v1/release.py "data/grs22-deepjet2000-v1-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
