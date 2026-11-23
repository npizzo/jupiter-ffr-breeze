#!/bin/bash
#SBATCH --job-name=ffr-matched-geometry
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/beta-shear-heated-v3-screen-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/beta_shear_heated_v3/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project=. scripts/beta_shear_heated_v3/driver.jl screen "data/beta-shear-heated-v3-screen-$SLURM_JOB_ID"
