#!/bin/bash
#SBATCH --job-name=ffr-gallery100-web
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/cascade100-web-%j.out
set -euo pipefail
tools/gallery-env/bin/python -u scripts/cascade100_web_encode.py
