#!/bin/bash
#SBATCH --job-name=ffr-rotation-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=logs/chan-rotation-v2-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_rotation_v2/gate.py preflight "data/chan-rotation-v2-cpu-$SLURM_JOB_ID"
