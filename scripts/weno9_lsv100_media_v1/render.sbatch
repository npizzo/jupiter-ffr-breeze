#!/bin/bash
#SBATCH --job-name=ffr-weno9-lsv100-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/weno9-lsv100-media-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m py_compile scripts/weno9_lsv100_media_v1/render.py
tools/point4-analysis-env/bin/python scripts/weno9_lsv100_media_v1/render.py
