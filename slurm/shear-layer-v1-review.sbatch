#!/bin/bash
#SBATCH --job-name=ffr-shear-layer-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/shear-layer-v1-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/shear_layer_publication_v1/prepare.jl "$1" "$2"
tools/point4-analysis-env/bin/python scripts/shear_layer_publication_v1/review.py "$1" "$2"
