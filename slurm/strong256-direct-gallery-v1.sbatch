#!/bin/bash
#SBATCH --job-name=ffr-strong256-direct-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --output=logs/strong256-direct-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHON=tools/point4-analysis-env/bin/python
OUT="rendered/strong256-direct-gallery-$SLURM_JOB_ID"
"$PYTHON" scripts/strong_shear_gallery_v1/assemble.py "$1" "$2" "$3" "$OUT/joined"
"$PYTHON" scripts/strong_shear_gallery_v1/render.py "$OUT"
