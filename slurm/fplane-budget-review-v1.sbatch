#!/bin/bash
#SBATCH --job-name=ffr-dense-budget-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/fplane-dense-review-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/point4-analysis-env/bin/python scripts/fplane_budget_v1/review.py data/fplane-budget-v1 rendered/fplane-budget-v1
