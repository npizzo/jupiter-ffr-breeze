# Gamma publication v2: separate vorticity validation and total-PV audit

This reporting-only revision retains the original diagnostic tolerances
(`rtol=5e-6`, `atol=1e-18`). It never changes a simulation or replaces a recorded
solver diagnostic.

The first flat-trap report failed the total-PV comparison at day 1, although
kinetic energy and relative-vorticity enstrophy agreed within the original
tolerance. All snapshot times, native dimensions, finite fields, paired initial
conditions, kinetic energy, and relative enstrophy remain hard requirements for
publication. A failed total-PV comparison is recorded as a failed audit, with
per-frame pass flags, failed days, maximum discrepancy, and an explicit
`COMPLETE_WITH_PV_DIAGNOSTIC_DISCREPANCY` publication status. The aggregate
`snapshot_crosscheck_passed` remains false if any selected frame fails.

The figures retain both solver and independently reconstructed total-PV series.
The README flags failures prominently. Publication of verified relative-vorticity
movies and descriptive statistics is not certification of the total-PV budget.
The cause of any mismatch remains unresolved unless a separate native audit
establishes it; it is not automatically attributed to Float32 rounding.

Run the probe on Unity first, then pass its numeric CSV to the unchanged native
audit to compare the exact completed checkpoint, reconstructed initial state,
and the worst saved-time mismatch:

```sh
sbatch slurm/gamma2026-publication-probe-v2.sbatch \
  data/gamma2026-v3 rendered/gamma2026-flat-probe-v2 \
  --case flat_n1024_seed11

# After the probe completes:
sbatch slurm/gamma2026-diagnostic-audit.sbatch \
  configs/gamma2026_v2/flat_n1024_seed11.toml \
  data/gamma2026-v3/flat_n1024_seed11 \
  rendered/gamma2026-flat-probe-v2/flat_n1024_seed11_diagnostics.csv \
  rendered/gamma2026-flat-native-audit-v2
```

Review both reports before publishing. A successful Slurm exit for the probe
means it completed the audit; read `probe.json` for the scientific check status.
The v2 movie renderer continues to use actual daily relative-vorticity fields,
shared fixed physical scales, full-frame decoding, hashes and frame-count checks.
