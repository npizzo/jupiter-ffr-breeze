# Fixed-region gamma observer

This separate observer reads completed frozen 75-day gamma2026-v3 cases and never evolves or modifies a simulation. All numerical tests, diagnostics and plotting run inside an allocated Unity CPU job. It accepts the completed 1024² seed 11/22/33 pairs and the later 2048² seed 11 pair. Every requested pair must be complete; partial publication is explicit through selecting fewer pairs, not dropping missing cases silently.

From the Unity repository root, use a fresh output directory:

```sh
sbatch slurm/gamma2026-regions.sbatch data/gamma2026-v3 rendered/gamma2026-regions-selftest-v1 --self-test-only
sbatch slurm/gamma2026-regions.sbatch data/gamma2026-v3 rendered/gamma2026-regions-v1 --pairs 1024:11 1024:22 1024:33
```

The optional first job tests only constructed fields and writes an explicitly analytical-only marker. The full observer also reruns those checks before reading model fields.

After the fine pair finishes, a separate fresh output can include `2048:11` along with the three coarse pairs. Existing output directories are rejected. The allocation is 4 CPUs / 24 GB / two hours; no job was submitted during source preparation.

The physical masks are fixed before analysis: core r < Lγ, annuli 1–3, 3–4, 4–5 Lγ, exterior r ≥ 5 Lγ, and cumulative interiors r < 3, 4, 5 Lγ. Boundaries are half-open. Exterior means the portion within the periodic square. Per region/time, the observer records area, signed circulation, separate positive/negative circulation magnitudes, weighted radii and Cartesian centroids, and relative-enstrophy integral/fraction. Fractions use whole-domain denominators; absent-sign/empty-mask moments remain null/blank. These are Eulerian field moments, not detected vortex cores or tracked objects. No mass, regional-energy attribution or inverse-cascade flux is inferred.

Before analysis, allocated numerical tests verify exact boundary inclusion, constructed constant/sign annuli, a two-dimensional lattice annulus, rotation/sign reversal, zero/absent-sign handling and partition/cumulative integral and moment identities. The same finite/partition checks run on every actual saved frame. The frozen publication_v2 energy/relative-enstrophy tests and paired initial-field checks are reused at their original thresholds. Unresolved total-PV failures remain failures in the report. All 76 native daily frames, original configurations and source revisions are checked, and input hashes are verified again after observation.

Outputs are per-case long-form CSVs, `regions.json`, `analytical_checks.json`, one fixed-region diagnostic figure, a readable report and the reused native diagnostic audit. No raw fields are exported. The report links the already verified synchronized 2D movies on the [GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/#gamma75); no new movie is generated.

Local preparation is limited to source editing, Python AST parsing, hashes and shell syntax. Actual analytical tests, data gates and plots remain pending Unity execution until the parent submits the wrapper.
