#!/bin/bash
#SBATCH --job-name=ffr-poly-form-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/polynomial-formation-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --project=. scripts/polynomial_formation_v1/check.jl data/viscous2d-refine-v1-1848240/initial.jld2 "data/polynomial-formation-check-${SLURM_JOB_ID}"
"$JULIA_BIN" --project=. scripts/polynomial_formation_v1/run.jl "data/polynomial-formation-check-${SLURM_JOB_ID}/driver" data/viscous2d-refine-v1-1848240/initial.jld2 .5 "data/polynomial-formation-check-${SLURM_JOB_ID}/PASS.toml" no_advection
tools/point4-analysis-env/bin/python - "$SLURM_JOB_ID" <<'PY'
import csv,json,sys,tomllib
from pathlib import Path
p=Path('data/polynomial-formation-check-'+sys.argv[1]);d=p/'driver'
assert tomllib.loads((d/'status.toml').read_text())['status']=='COMPLETED'
assert [float(r['t']) for r in csv.DictReader((d/'scalars.csv').open())]==[0,.25,.5]
assert [float(r['t']) for r in csv.DictReader((d/'frames.csv').open())]==[0,.5]
(p/'DRIVER_PASS.json').write_text(json.dumps({'status':'PASS','target':.5}))
PY
