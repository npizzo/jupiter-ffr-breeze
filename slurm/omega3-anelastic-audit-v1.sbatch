#!/bin/bash
#SBATCH --job-name=ffr-anelastic-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=logs/omega3-anelastic-audit-v1-%j.out
set -euo pipefail
tools/point4-analysis-env/bin/python scripts/omega3_anelastic_audit_v1/review.py "data/omega3-anelastic-audit-v1-$SLURM_JOB_ID" "$1"
