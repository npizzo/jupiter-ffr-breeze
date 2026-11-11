#!/bin/bash
# Unity only; all numerical scripts/configurations remain the frozen SI-informed v2.
# Without an optional SETUP_JOB, first submit the isolated GPU-runtime setup.
set -euo pipefail
root="${1:?output root}";setup="${2:-}"
mkdir -p "$root" logs
record="$root/submission.tsv"
test ! -e "$record" || { echo "Submission record exists; refusing duplicate chain";exit 1; }
printf 'case\tjob\tafterok\n' > "$record"
if [[ -z "$setup" ]]; then
  setup=$(sbatch --parsable --kill-on-invalid-dep=yes slurm/gamma2026_gpu-setup.sbatch "$root")
fi
printf 'runtime_setup\t%s\t-\n' "$setup" >> "$record"
unit=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$setup" slurm/gamma2026_gpu-unit.sbatch "$root")
printf 'unit\t%s\t%s\n' "$unit" "$setup" >> "$record"
gpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/gamma2026_gpu-gpu.sbatch "$root" 1024)
printf 'gpu_1024\t%s\t%s\n' "$gpu" "$unit" >> "$record"
previous="$gpu"
for seed in 11 22 33; do
  for kind in polar flat; do
    case_name="${kind}_n1024_seed${seed}"
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_gpu-run.sbatch "$root" "$case_name" 1024)
    printf '%s\t%s\t%s\n' "$case_name" "$job" "$previous" >> "$record"
    previous="$job"
  done
done
gpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_gpu-gpu.sbatch "$root" 2048)
printf 'gpu_2048\t%s\t%s\n' "$gpu" "$previous" >> "$record"
previous="$gpu"
for kind in polar flat; do
  case_name="${kind}_n2048_seed11"
  job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_gpu-run.sbatch "$root" "$case_name" 2048)
  printf '%s\t%s\t%s\n' "$case_name" "$job" "$previous" >> "$record"
  previous="$job"
done
cat "$record"
