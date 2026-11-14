#!/bin/bash
#SBATCH --job-name=ffr-point4-analysis-smoke
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fplane-point4-smoke-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
test -f data/fplane-point4-v2/preflight/PASS.toml
tools/point4-analysis-env/bin/python scripts/fplane_point4_v2/review.py data/fplane-point4-v2/preflight late_heated rendered/fplane-point4-analysis-validation-v3
