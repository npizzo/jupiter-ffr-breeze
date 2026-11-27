#!/bin/bash
#SBATCH --job-name=convective-cap-v2-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/convective-cap-v2-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-convective-cap-v2"
tools/point4-analysis-env/bin/python scripts/convective_deformation_results_v2/render.py
