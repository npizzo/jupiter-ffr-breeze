#!/bin/bash
#SBATCH --job-name=ffr-plume-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --output=logs/plume-pair-v1-validate-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/plume_pair_v1/validate.jl "data/plume-pair-v1-validation-$SLURM_JOB_ID"
