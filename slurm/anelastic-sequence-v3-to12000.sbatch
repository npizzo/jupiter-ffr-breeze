#!/bin/bash
#SBATCH --job-name=ffr-anelastic-to12000
#SBATCH --partition=uri-gpu,gpu,gpu-highvram
#SBATCH --qos=short
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/anelastic-sequence-v3-to12000-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
root="data/anelastic-sequence-v3-to12000-$SLURM_JOB_ID"
trap 'code=$?; python scripts/anelastic_sequence_v3_to12000/launch_gate.py outcome "$root" "$code"; exit "$code"' EXIT
python scripts/anelastic_sequence_v3_to12000/launch_gate.py verify "$1" "$2"
printf 'GPU_PROCESS_START job=%s\n' "$SLURM_JOB_ID"
"$JULIA_BIN" --project=. scripts/anelastic_sequence_v3_to12000/pilot.jl "$root" "$2"
