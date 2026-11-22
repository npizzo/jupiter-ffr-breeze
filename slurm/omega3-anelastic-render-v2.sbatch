#!/bin/bash
#SBATCH --job-name=ffr-anelastic-render
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=logs/omega3-anelastic-render-v2-%j.out
set -euo pipefail
tools/point4-analysis-env/bin/python scripts/omega3_anelastic_render_v2/render.py "data/omega3-anelastic-render-v2-$SLURM_JOB_ID"
