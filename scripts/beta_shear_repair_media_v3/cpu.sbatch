#!/bin/bash
#SBATCH --job-name=ffr-repair-readable-colorbars
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=logs/beta-shear-repair-media-v3-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/beta_shear_repair_media_v3/render.py "data/beta-shear-repair-media-v3-$SLURM_JOB_ID" data/beta-shear-repair-v5-run-2140541 data/beta-shear-profile-v1-run-2140644
