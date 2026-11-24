#!/bin/bash
#SBATCH --job-name=ffr-chan-clean-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-clean-pair-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_clean_pair_v1/gate.py preflight "data/chan-clean-pair-v1-cpu-$SLURM_JOB_ID"
python scripts/chan_clean_pair_v1/release.py "data/chan-clean-pair-v1-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"
