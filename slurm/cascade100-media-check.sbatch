#!/bin/bash
#SBATCH --job-name=ffr-gallery100-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/cascade100-media-check-%j.out
set -euo pipefail
tools/gallery-env/bin/python scripts/verify_cascade100_media.py
