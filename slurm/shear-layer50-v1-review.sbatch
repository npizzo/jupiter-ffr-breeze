#!/bin/bash
#SBATCH --job-name=ffr-shear-layer-50movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/shear-layer50-v1-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHON=tools/point4-analysis-env/bin/python
"$PYTHON" scripts/shear_layer50_publication_v1/assemble.py "$1" "$2" "$3/joined"
"$JULIA_BIN" --project=. scripts/shear_layer50_publication_v1/prepare.jl "$3/joined" "$3"
"$PYTHON" scripts/shear_layer50_publication_v1/review.py "$3/joined" "$3"
