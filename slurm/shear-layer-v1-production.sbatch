#!/bin/bash
#SBATCH --job-name=ffr-shear-layer
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/shear-layer-v1-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cases=(stable_cyclonic stable_none stable_anticyclonic neutral_cyclonic neutral_none neutral_anticyclonic)
"$JULIA_BIN" --project=. scripts/shear_layer_v1/cuda_launch.jl campaign "$1" "${cases[$SLURM_ARRAY_TASK_ID]}"
