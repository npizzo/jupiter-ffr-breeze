#!/bin/bash
# Run from repository root on Unity: bash scripts/gamma2026/submit.sh ROOT SETUP_JOB
set -euo pipefail
root="${1:?output root}";setup="${2:?successful setup job dependency}"
mkdir -p "$root" logs
record="$root/submission.tsv"
test ! -e "$record" || { echo "Submission record exists; refusing duplicate chain";exit 1; }
printf 'case\tjob\tafterok\n' > "$record"
unit=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$setup" slurm/gamma2026-unit.sbatch "$root")
printf 'unit\t%s\t%s\n' "$unit" "$setup" >> "$record"
gpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/gamma2026-gpu.sbatch "$root")
printf 'gpu_gate\t%s\t%s\n' "$gpu" "$unit" >> "$record"
previous="$gpu"
for seed in 11 22 33; do
  for kind in polar flat; do
    case_name="${kind}_seed${seed}"
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026-run.sbatch "$root" "$case_name")
    printf '%s\t%s\t%s\n' "$case_name" "$job" "$previous" >> "$record"
    previous="$job"
  done
done
cat "$record"
