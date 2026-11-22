#!/bin/bash
#SBATCH --job-name=ffr-response-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/chan-rotation-response-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/chan_rotation_response_v1/gate.py preflight "data/chan-rotation-response-v1-cpu-$SLURM_JOB_ID"
