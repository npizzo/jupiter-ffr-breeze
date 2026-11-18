#!/bin/bash
#SBATCH --job-name=ffr-stable50-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/shear-stable50-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHON=tools/point4-analysis-env/bin/python
OUT="rendered/shear-stable50-gallery-$SLURM_JOB_ID"
"$PYTHON" scripts/shear_stable50_gallery_v1/assemble.py data/shear-layer-v1-1962402 data/shear-layer50-v1-1962402 "$OUT/joined"
"$PYTHON" scripts/shear_stable50_gallery_v1/render.py "$OUT"
