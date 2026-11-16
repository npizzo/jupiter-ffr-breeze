#!/bin/bash
#SBATCH --job-name=ffr-LSV-instability-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/fplane-centrifugal-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=tools/point4-analysis-env/bin/python
OUT=rendered/fplane-centrifugal-v1
$PY scripts/fplane_centrifugal_v1/analyze.py data/fplane-point4-v2/late_heated "$OUT"
$PY scripts/fplane_centrifugal_v1/plot.py "$OUT"
