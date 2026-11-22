#!/bin/bash
#SBATCH --job-name=ffr-omega15-long-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/omega15-long-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega15_long_v1/gate.py preflight "data/omega15-long-v1-cpu-$SLURM_JOB_ID"
