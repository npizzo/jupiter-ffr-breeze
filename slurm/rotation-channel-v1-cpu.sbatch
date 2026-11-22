#!/bin/bash
#SBATCH --job-name=ffr-channel-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=logs/rotation-channel-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/rotation_channel_v1/gate.py preflight "data/rotation-channel-v1-cpu-$SLURM_JOB_ID"
