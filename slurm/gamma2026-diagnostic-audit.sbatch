#!/bin/bash
#SBATCH --job-name=ffr-gamma-pv-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/gamma2026-diagnostic-audit-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-gamma2026-gpu:$FFR_ROOT/depot-gamma2026"
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
"$JULIA_BIN" --startup-file=no --project=environments/gamma2026_gpu scripts/gamma2026_diagnostic_audit/audit.jl "$@"
tools/gallery-env/bin/python scripts/gamma2026_diagnostic_audit/to_json.py "$4"
