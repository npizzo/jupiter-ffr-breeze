#!/bin/bash
#SBATCH --job-name=chan-restart-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=00:20:00
#SBATCH --output=logs/chan-restart-audit-v1-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
tools/point4-analysis-env/bin/python scripts/chan_clean_anelastic6000_v1/health.py verify runs/chan-clean-anelastic6000-v1/HEALTH_PASS.json
tools/point4-analysis-env/bin/python scripts/chan_theta_diagnostic4000_v1/health.py verify runs/chan-theta-diagnostic4000-v1/HEALTH_PASS.json
"$JULIA_BIN" --pkgimages=no --project=. scripts/chan_restart_audit_v1/audit.jl "data/chan-restart-audit-v1-${SLURM_JOB_ID}"
