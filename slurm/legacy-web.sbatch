#!/bin/bash
#SBATCH --job-name=ffr-legacy-web
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/legacy-web-%j.out
set -euo pipefail
tools/gallery-env/bin/python -u scripts/legacy_web_encode.py
