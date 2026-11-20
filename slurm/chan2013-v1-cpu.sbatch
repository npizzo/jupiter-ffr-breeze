#!/bin/bash
#SBATCH --job-name=ffr-chan-cpu-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=logs/chan2013-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/chan2013_v1/preflight_cpu.jl "data/chan2013-v1-cpu-$SLURM_JOB_ID"
