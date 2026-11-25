#!/bin/bash
#SBATCH --job-name=chan-three-case-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-three-case2000-render-v1-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-chan-three-case2000"
tools/point4-analysis-env/bin/python scripts/chan_three_case2000_publication_v1/render.py
