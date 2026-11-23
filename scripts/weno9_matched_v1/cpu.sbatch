#!/bin/bash
#SBATCH --job-name=ffr-weno9-matched-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=logs/weno9-matched-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/weno9_matched_v1/gate.py preflight "data/weno9-matched-v1-cpu-$SLURM_JOB_ID"
