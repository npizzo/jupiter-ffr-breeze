#!/bin/bash
#SBATCH --job-name=ffr-point4-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fplane-point4-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
if [ "$1" = compare ]; then
 tools/gallery-env/bin/python scripts/fplane_point4_v1/compare.py "$2" "$3"
else
 mapfile -t cases < scripts/fplane_point4_v1/case-order.txt
 tools/gallery-env/bin/python scripts/fplane_point4_v1/review.py "$1" "${cases[$SLURM_ARRAY_TASK_ID]}" "$2"
fi
