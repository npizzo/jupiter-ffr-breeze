#!/bin/bash
#SBATCH --job-name=ffr-plume-interface-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --output=logs/plume-interface-audit-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/ffr-mpl-$SLURM_JOB_ID
tools/point4-analysis-env/bin/python scripts/plume_interface_audit_v1/analyze.py rendered/plume-interface-audit-v1-$SLURM_JOB_ID
