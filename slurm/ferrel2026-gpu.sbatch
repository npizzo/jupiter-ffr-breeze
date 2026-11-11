#!/bin/bash
#SBATCH --job-name=ffr-ferrgpu
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel2026-gpu-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026_gpu/gpu_job.sh "$3" scripts/ferrel2026_v3/campaign.jl "$1" "$2" "$3"
tools/field-export-env/bin/python scripts/ferrel2026_v3/transfer.py "$2/ferrel2026_$1"
