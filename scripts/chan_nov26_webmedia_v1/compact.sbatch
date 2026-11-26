#!/bin/bash
#SBATCH --job-name=chan-web-media
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=logs/chan-nov26-web-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/chan_nov26_webmedia_v1/compact.py
