#!/bin/bash
#SBATCH --job-name=ffr-clean-comp-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-clean-compressible2000-render-v2-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-clean-comp2000"
tools/point4-analysis-env/bin/python scripts/chan_clean_compressible_publication_v2/render.py
