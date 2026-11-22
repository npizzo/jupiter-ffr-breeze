#!/bin/bash
#SBATCH --job-name=ffr-vortex-thermal
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/vortex-thermal-v2-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="/tmp/ffr-mpl-$SLURM_JOB_ID"
exec tools/point4-analysis-env/bin/python scripts/vortex_thermal_v2/analyze.py "data/vortex-thermal-v2-$SLURM_JOB_ID/publication"
