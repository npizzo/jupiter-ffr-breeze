#!/bin/bash
#SBATCH --job-name=ffr-point4-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fplane-point4-review-v3-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
if [ "$1" = compare ]; then
 tools/point4-analysis-env/bin/python scripts/fplane_point4_v2/compare.py "$2" "$3"
else
 mapfile -t cases < scripts/fplane_point4_v2/case-order.txt
 tools/point4-analysis-env/bin/python scripts/fplane_point4_v2/review.py "$1" "${cases[$SLURM_ARRAY_TASK_ID]}" "$2"
fi
