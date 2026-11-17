#!/bin/bash
#SBATCH --job-name=ffr-double310-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/viscous2d-double310-movie-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/viscous2d_double310_v1/render.py "$1" "$2" "$3"
