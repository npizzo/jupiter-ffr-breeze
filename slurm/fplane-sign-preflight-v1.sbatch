#!/bin/bash
#SBATCH --job-name=ffr-sign-budget-gate
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=logs/fplane-sign-preflight-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/fplane_sign_v1/cuda_launch.jl preflight data/fplane-sign-v1
tools/point4-analysis-env/bin/python scripts/fplane_sign_v1/gate.py data/fplane-sign-v1
