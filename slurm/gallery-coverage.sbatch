#!/bin/bash
#SBATCH --job-name=ffr-gallery-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/gallery-coverage-%j.out
set -euo pipefail
tools/gallery-env/bin/python scripts/audit_gallery_coverage.py
