#!/bin/bash
#SBATCH --job-name=chan-theta-diagnostic400-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-theta-diagnostic4000-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_theta_diagnostic4000_v1/gate.py preflight "data/chan-theta-diagnostic4000-v1-cpu-$SLURM_JOB_ID"
