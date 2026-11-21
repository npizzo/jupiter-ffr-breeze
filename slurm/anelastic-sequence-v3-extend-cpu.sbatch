#!/bin/bash
#SBATCH --job-name=ffr-exact-launch-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --output=logs/anelastic-sequence-v3-extend-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec python scripts/anelastic_sequence_v3_extend/launch_gate.py preflight "data/anelastic-sequence-v3-extend-cpu-$SLURM_JOB_ID" "$1"
