#!/bin/bash
#SBATCH --job-name=ffr-lsv1000-layout
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/weno9-lsv1000-layout-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/weno9_lsv1000_media_v1/fix_anelastic_layout.py
