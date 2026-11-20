#!/bin/bash
#SBATCH --job-name=ffr-compressible-pilots
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=01:02:00
#SBATCH --output=logs/compressible-pilots-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=16
export JULIA_NUM_PRECOMPILE_TASKS=4
root="data/compressible-pilots-v1-$SLURM_JOB_ID"
mkdir -p "$root"
status=0
for family in kapyla chan; do
  code=0
  timeout --signal=TERM --kill-after=30s 1800 "$JULIA_BIN" --project=. scripts/compressible_pilots_v1/run.jl "$family" "$root/$family" > "$root/$family.log" 2>&1 || code=$?
  printf '%s %s\n' "$family" "$code" >> "$root/exit-codes.txt"
  if [ "$code" -ne 0 ]; then status=1; fi
  echo "Pilot $family exited $code; production remains held."
done
exit "$status"
