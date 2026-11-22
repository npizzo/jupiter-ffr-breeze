#!/bin/bash
#SBATCH --job-name=ffr-vortex-flux-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/vortex-flux-audit-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec tools/point4-analysis-env/bin/python scripts/vortex_thermal_v2/flux_audit.py data/vortex-thermal-v2-2108995/publication
