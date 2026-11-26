#!/bin/bash
#SBATCH --job-name=deformation-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/deformation-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-deformation-movies"
tools/point4-analysis-env/bin/python scripts/deformation_pair_results_v1/render.py
