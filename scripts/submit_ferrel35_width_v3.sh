#!/bin/bash
# Run from the Unity repository root after deploying this isolated campaign.
# Usage: bash scripts/submit_ferrel35_width_v3.sh CPU_VALIDATION_JOB
set -euo pipefail
cpu="$1"
root=data/ferrel35-width-v3
preflight=data/ferrel35-width-preflight-v3
record=campaigns/ferrel35-width-v3/jobs.tsv
mkdir -p campaigns/ferrel35-width-v3
test ! -e "$record"
printf 'case\tsegment\tjob_id\tdependency\n' > "$record"
submit() { sbatch --parsable --kill-on-invalid-dep=yes "$@"; }
gpu=$(submit --dependency="afterok:$cpu" slurm/ferrel35-width-v3-preflight.sbatch "$preflight")
printf 'gpu_validation\t0\t%s\t%s\n' "$gpu" "$cpu" >> "$record"
diag=$(submit --job-name=ffr-ferrel35-budget210 --dependency="afterok:$gpu" slurm/ferrel35-width-v3.sbatch diagnostic 1 "$root" "$preflight")
printf 'diagnostic\t1\t%s\t%s\n' "$diag" "$gpu" >> "$record"
movie=$(submit --job-name=ffr-ferrel35-budget-movie --dependency="afterok:$diag" slurm/ferrel35-width-gallery.sbatch diagnostic "$root" rendered/ferrel35-budget200-210-v3)
printf 'diagnostic_movie\t0\t%s\t%s\n' "$movie" "$diag" >> "$record"
previous="$diag"
for segment in $(seq 1 20); do
 job=$(submit --job-name="ffr-ferrel35-wide-d$((segment*10))" --dependency="afterok:$previous" slurm/ferrel35-width-v3.sbatch wide "$segment" "$root" "$preflight")
 printf 'wide\t%s\t%s\t%s\n' "$segment" "$job" "$previous" >> "$record"
 previous="$job"
done
movie=$(submit --job-name=ffr-ferrel35-wide-movie --dependency="afterok:$previous" slurm/ferrel35-width-gallery.sbatch wide "$root" rendered/ferrel35-wide0-200-v3)
printf 'wide_movie\t0\t%s\t%s\n' "$movie" "$previous" >> "$record"
cat "$record"
