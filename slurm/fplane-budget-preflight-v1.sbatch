#!/bin/bash
#SBATCH --job-name=ffr-vorticity-budget-gate
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram80
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fplane-budget-preflight-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/fplane_budget_v1/cuda_launch.jl preflight data/fplane-budget-v1
tools/point4-analysis-env/bin/python scripts/fplane_budget_v1/gate.py data/fplane-budget-v1
