#!/bin/bash
#SBATCH --job-name=ffr-omega45-long-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/omega45-long-v2-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega45_long_v2/gate.py preflight "data/omega45-long-v2-cpu-$SLURM_JOB_ID"
