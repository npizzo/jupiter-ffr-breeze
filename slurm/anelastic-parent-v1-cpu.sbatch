#!/bin/bash
#SBATCH --job-name=ffr-anelastic-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=logs/anelastic-parent-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/anelastic_parent_v1/preflight_cpu.jl "data/anelastic-parent-v1-cpu-$SLURM_JOB_ID"
