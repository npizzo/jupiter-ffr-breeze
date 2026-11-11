#!/bin/bash
#SBATCH --job-name=ffr-conv26
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/convection26-%j.out
#SBATCH --no-requeue
set -euo pipefail
bash scripts/convection2026/gpu_job.sh scripts/convection2026/campaign.jl "$1" "$2" "$3" "$4"
tools/field-export-env/bin/python scripts/convection2026/transfer.py "$3/convection2026_f$1/segment_$2"
