#!/bin/bash
#SBATCH --job-name=ffr-kapyla-cpu-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=logs/kapyla-pair-v1-cpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
exec "$JULIA_BIN" --project=. scripts/kapyla_pair_v1/preflight_cpu.jl "data/kapyla-pair-v1-cpu-$SLURM_JOB_ID"
