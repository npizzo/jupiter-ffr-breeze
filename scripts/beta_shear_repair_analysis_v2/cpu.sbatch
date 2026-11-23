#!/bin/bash
#SBATCH --job-name=ffr-repair-analysis-render
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=logs/beta-shear-repair-analysis-v2-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
out="data/beta-shear-repair-analysis-v2-$SLURM_JOB_ID"
tools/point4-analysis-env/bin/python scripts/beta_shear_repair_analysis_v2/analyze.py "$out" data/beta-shear-repair-v5-run-2140541 data/beta-shear-repair-v5-forensic-2140542 data/beta-shear-profile-v1-run-2140644 data/beta-shear-spectrum-v3-2140690 data/heated-energy-audit-v1-2140833
tools/point4-analysis-env/bin/python scripts/beta_shear_repair_analysis_v2/render.py "$out" data/beta-shear-repair-v5-run-2140541 data/beta-shear-profile-v1-run-2140644
