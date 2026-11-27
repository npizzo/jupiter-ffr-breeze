#!/bin/bash
#SBATCH --job-name=convective-cap-v2-assessment
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/convective-cap-v2-assessment-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/convective_deformation_v2/assess.py "$@"
