#!/bin/bash
#SBATCH --job-name=ffr-fixed45-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/chan-fixed45-v2-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_fixed45_v2/gate.py preflight "data/chan-fixed45-v2-cpu-$SLURM_JOB_ID"
