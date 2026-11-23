#!/bin/bash
#SBATCH --job-name=ffr-linear-y-refinement
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --output=logs/beta-shear-spectrum-v3-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export OPENBLAS_NUM_THREADS=1
out="data/beta-shear-spectrum-v3-$SLURM_JOB_ID"
"$JULIA_BIN" --project=. scripts/beta_shear_spectrum_v3/export.jl "$out"
tools/point4-analysis-env/bin/python scripts/beta_shear_spectrum_v3/spectrum.py "$out" "$out"
