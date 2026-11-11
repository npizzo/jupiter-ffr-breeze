#!/bin/bash
#SBATCH --job-name=ffr-ferr26v2
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel26v2-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/ferrel2026_v2/gpu_job.sh scripts/ferrel2026_v2/campaign.jl "$1" "$2" "$3"
tools/field-export-env/bin/python scripts/ferrel2026_v2/transfer.py "$2/ferrel2026_$1"
