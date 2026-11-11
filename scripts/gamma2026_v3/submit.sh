#!/bin/bash
# Unity only; uses the verified gamma2026_gpu runtime environment unchanged.
# No dependency on an aged-out environment setup job. All scientific gates rerun.
set -euo pipefail
root="${1:?output root}"
mkdir -p "$root" logs
record="$root/submission.tsv"
test ! -e "$record" || { echo "Submission record exists; refusing duplicate chain";exit 1; }
printf 'case\tjob\tafterok\n' > "$record"
unit=$(sbatch --parsable --kill-on-invalid-dep=yes slurm/gamma2026_v3-unit.sbatch "$root")
printf 'unit\t%s\tverified_runtime\n' "$unit" >> "$record"
gpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$unit" slurm/gamma2026_v3-gpu.sbatch "$root" 1024)
printf 'gpu_1024\t%s\t%s\n' "$gpu" "$unit" >> "$record"
previous="$gpu"
for seed in 11 22 33; do
  for kind in polar flat; do
    case_name="${kind}_n1024_seed${seed}"
    job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_v3-run.sbatch "$root" "$case_name" 1024)
    printf '%s\t%s\t%s\n' "$case_name" "$job" "$previous" >> "$record"
    previous="$job"
  done
done
gpu=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_v3-gpu.sbatch "$root" 2048)
printf 'gpu_2048\t%s\t%s\n' "$gpu" "$previous" >> "$record"
previous="$gpu"
for kind in polar flat; do
  case_name="${kind}_n2048_seed11"
  job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$previous" slurm/gamma2026_v3-run.sbatch "$root" "$case_name" 2048)
  printf '%s\t%s\t%s\n' "$case_name" "$job" "$previous" >> "$record"
  previous="$job"
done
cat "$record"
