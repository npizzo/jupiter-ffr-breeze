#!/bin/bash
#SBATCH --job-name=ffr-omega3-anelastic-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/omega3-anelastic-v3-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/omega3_anelastic_v3/gate.py preflight "data/omega3-anelastic-v3-cpu-$SLURM_JOB_ID"
