#!/bin/bash
#SBATCH --job-name=ffr-tracer100-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/tracer100-audit-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
tools/gallery-env/bin/python scripts/tracer100_contrast_audit.py
