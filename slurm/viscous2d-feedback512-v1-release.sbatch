#!/bin/bash
#SBATCH --job-name=ffr-feedback-release
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/viscous2d-feedback512-release-%j.out
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/viscous2d_feedback512_v1/run.jl "$1/driver_smoke" "$2" 1 "$1/PASS.toml" biased_plus
tools/point4-analysis-env/bin/python - "$1/driver_smoke" <<'PY'
from pathlib import Path
import csv,math,sys,tomllib
p=Path(sys.argv[1]);s=tomllib.loads((p/'status.toml').read_text());assert s['status']=='COMPLETED' and abs(s['t']-1)<1e-8
frames=list(csv.DictReader((p/'frames.csv').open()));assert len(frames)==3
assert all(abs(float(row['t'])-i*.5)<1e-8 for i,row in enumerate(frames))
for row in frames:assert (p/row['file']).stat().st_size==512*512*4
scalars=list(csv.DictReader((p/'scalars.csv').open()));assert len(scalars)==3
assert all(math.isfinite(float(v)) for row in scalars for v in row.values())
assert (p/'checkpoint.jld2').stat().st_size>0
(p.parent/'DRIVER_PASS.txt').write_text('Production driver, time-aligned fields, scalar diagnostics and checkpoint: PASS\n')
print('DRIVER_PASS',flush=True)
PY
