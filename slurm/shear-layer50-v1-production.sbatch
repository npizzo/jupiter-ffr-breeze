#!/bin/bash
#SBATCH --job-name=ffr-shear-layer-to50
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=l40s
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:30:00
#SBATCH --output=logs/shear-layer50-v1-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cases=(stable_cyclonic stable_none stable_anticyclonic neutral_cyclonic neutral_none neutral_anticyclonic)
"$JULIA_BIN" --project=. scripts/shear_layer50_v1/continue.jl "$1" "$2" "${cases[$SLURM_ARRAY_TASK_ID]}"
