#!/bin/bash
#SBATCH --job-name=ffr-existing-heated-energy
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=logs/heated-energy-audit-v1-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --project=. scripts/heated_energy_audit_v1/driver.jl "data/heated-energy-audit-v1-$SLURM_JOB_ID"
