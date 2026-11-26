#!/bin/bash
#SBATCH --job-name=ffr-deformation-init
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=logs/deformation-pair-v1-phase-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --pkgimages=no --project=. scripts/deformation_pair_v1/choose_phase.jl "data/deformation-pair-v1-phase-$SLURM_JOB_ID"
