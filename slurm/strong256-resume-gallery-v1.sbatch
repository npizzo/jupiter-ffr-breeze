#!/bin/bash
#SBATCH --job-name=ffr-strong256-resume-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --output=logs/strong256-resume-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHON=tools/point4-analysis-env/bin/python
OUT="rendered/strong256-resume-gallery-$SLURM_JOB_ID"
"$PYTHON" scripts/strong256_resume_gallery_v1/assemble.py "$1" "$2" "$3" "$OUT/joined"
"$PYTHON" scripts/strong_shear_gallery_v1/render.py "$OUT"
