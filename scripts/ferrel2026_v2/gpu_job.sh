#!/bin/bash
# Run only inside a Unity allocation. Keep per-process GPU memory observations.
set -euo pipefail
source ./local.env
export JULIA_DEPOT_PATH="$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4
export JULIA_NUM_PRECOMPILE_TASKS=4
export OPENBLAS_NUM_THREADS=1
memory_log="logs/ferrel2026-v2-memory-${SLURM_JOB_ID}.csv"
hardware_log="logs/ferrel2026-v2-hardware-${SLURM_JOB_ID}.csv"
nvidia-smi --query-gpu=uuid,name,memory.total --format=csv,noheader,nounits > "$hardware_log"
printf 'epoch_seconds,pid,gpu_uuid,gpu_memory_MiB\n' > "$memory_log"
"$JULIA_BIN" --project=. "$@" &
model_pid=$!
(
  while kill -0 "$model_pid" 2>/dev/null; do
    stamp=$(date +%s)
    nvidia-smi --query-compute-apps=pid,gpu_uuid,used_gpu_memory --format=csv,noheader,nounits 2>/dev/null | awk -F, -v wanted="$model_pid" -v stamp="$stamp" '$1+0==wanted {gsub(/ /,"",$2);gsub(/ /,"",$3);print stamp "," wanted "," $2 "," $3}' >> "$memory_log" || true
    sleep 15
  done
) &
monitor_pid=$!
trap 'kill "$monitor_pid" 2>/dev/null || true' EXIT
wait "$model_pid"
python3 - "$memory_log" "$hardware_log" <<'PY'
import csv,sys
with open(sys.argv[2]) as f:
    capacity={r[0].strip():float(r[2]) for r in csv.reader(f)}
with open(sys.argv[1]) as f:
    rows=list(csv.DictReader(f))
if not rows:
    raise SystemExit('No per-process GPU memory samples: allocation validation incomplete')
ratios=[float(r['gpu_memory_MiB'])/capacity[r['gpu_uuid']] for r in rows]
peak=max(ratios)
print(f'Sampled peak process GPU memory fraction: {peak:.4f}',flush=True)
if peak>=0.95:
    raise SystemExit('Insufficient sampled GPU-memory headroom (<5%): block dependent production')
PY
