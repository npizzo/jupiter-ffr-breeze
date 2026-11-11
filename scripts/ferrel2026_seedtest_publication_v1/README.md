# Read-only publication for the four horizontal-seed tests

This new wrapper preserves every frozen scientific driver, existing observer and
renderer. It imports the hash-pinned Ferrel resolution primitives and v2 native
collector/movie encoder. Only in-memory presentation labels change. All reading,
diagnostics, rendering, analytical checks and media validation run on allocated
Unity CPUs. Nothing is submitted, cancelled, uploaded or modified in model data.

## Required completion and provenance

All four new cases must be complete: `zero512_b0`, `zero512_b162`,
`modes512_b162`, `modes1024_b162`. Each must have its own matching `CASE_PASS.toml`
and a successful enclosing Slurm job (`COMPLETED`, exit0:0). This is important:
the model's `COMPLETE` flag precedes final checks, and the GPU-memory-headroom
check runs after `CASE_PASS` is written. The initialization job must also have
completed successfully. Slurm IDs come from the gates, not arbitrary supplied
job numbers.

`guards.py` verifies complete source/configuration inventories, the recorded
initialization gate and common mode table, both prior v3 grid gates, the repaired
runtime Project/Manifest/preferences and package identities, the runtime guard
sources, the inheritance audit, executed configurations and metadata. The
per-run config may differ from its frozen template only by the verified mode
reference path and SHA. Model inputs and observer/helper source fingerprints are
rechecked before completion. The observer refuses an existing output directory.

## Diagnostics and precision

Each new run supplies five Float64 centered-w / temperature-anomaly planes at
**5.625,9.375,19.375,59.375,99.375km** above the model bottom. Their saved integer
indices, actual heights, shapes and overlap with the three video planes are
verified explicitly. The lowest video's absolute Float64 temperature is checked
against the corresponding anomaly. All five plane covariances and vertical
velocity powers must reproduce the native Float64 heat-flux and vertical-KE
profiles within derived arithmetic allowances. The anomaly is already relative
to the model's plane mean; that mean is not subtracted twice.

The spectra reuse the tested full complex FFT on periodic x and even-reflected
y. All Fourier modes are retained. Wavelength2π/k and inverse wavenumber1/k remain
separate. Both plane-mean-removed and zonal-mean-removed results are saved.
Signed heat cospectra and absolute mode-throughput tails are reported at4/8/16
cells, using radial and component criteria. Fractions apply to individual hourly
snapshots and each grid's own cell width; the code does not average fractions or
pretend their physical cutoffs are identical. It time-averages spectra and
profiles over0–6h only. No results are described as developed convection.

Float64 anomaly storage avoids the earlier Float32 absolute-temperature loss.
It does not remove Float64 precision from the dynamics. The spectral screen
therefore separately labels a one-rounding representation-sensitivity allowance:
half-ULP of absolute Float64 theta plus its anomaly, and half-ULP of centered w.
It is **not** an output quantization error or a bound on accumulated solver
error. Constant plane-mean errors are removed by the orthogonal projections.
Stored-profile identities use the saved Float64 arrays directly, with no added
format-quantization term. The inherited FFT arithmetic model assumes a stable
Float64 implementation; it is not certified interval arithmetic. Weak or zero
signals remain indeterminate where appropriate.

Allocated analytical checks retain the old known-wavelength/sign/Parseval and
Float32 weak-signal tests, and add a weak Float64 anomaly that must remain
detectable, a zero anomaly that must remain indeterminate, and a direct
anomaly/profile check that rejects corrupted heat flux. These are not physical
model convergence tests.

## Movies, comparisons and small publication assets

Four individual movies plus three paired movies use one set of fixed scales
across **all four cases and all seven saved frames**. Each has three actual depth
columns, with vorticity above tracer. The video planes remain
19.375/59.375/99.375km; the added lower diagnostic planes are not passed off as
extra movie depths. The frozen encoder preserves asinh widths0.1 inζ/f₀ and0.01
in tracer, full shared extrema, exact frame counts, browser-safe H.264/yuv420p,
decode validation, dimensions and hashes. Every asset must be below24MB.

The three comparisons have different interpretations:

| Pair | Allowed interpretation |
|---|---|
| Zero-seed heated minus zero-seed unheated | Heating response without an imposed thermal perturbation, including adjustment/heating interactions |
| Heated modal512 minus heated zero-seed512 | Sensitivity to the thermal seed |
| Heated modal1024 minus heated modal512 | Horizontal-grid sensitivity of the shared thermal seed; velocity discretization also changes |

There is no unheated modal case, so modal changes cannot be called a heating
effect. A band-limited seed may delay shorter unstable modes; a clean result
does not establish resolution of intrinsic convection. Existing cellwise-noise
pilots remain separately published references. This wrapper does not relabel
their three planes as five, manufacture their missing Float64 fields, or use
them as a missing matched control.

Outputs include seven movies and posters, per-movie metadata, four comparison
figures, five compact CSV tables per case, native mean-flow CSVs, an analytical
check report, a provenance-rich `report.json`, and a hashed `manifest.json`.
Only the manifest's `READY_FOR_REVIEW` result is eligible for parent review and
GitHub publication. The actual gallery remains
[the repository's GitHub Pages site](https://npizzo.github.io/jupiter-ffr-breeze/).

## Parent execution

First run an optional allocated reader self-check, before model cases finish:

```sh
sbatch slurm/ferrel2026-seedtest-publication-v1.sbatch \
  --self-test-only --output rendered/ferrel2026-seed-reader-qa-v1
```

After all four model jobs succeed, submit once with dependencies on those jobs:

```sh
sbatch slurm/ferrel2026-seedtest-publication-v1.sbatch \
  --production-root data/ferrel2026-seedtest-gpu-v1 \
  --gate-root data/ferrel2026-seedtest-preflight-v1 \
  --prior-root data/ferrel2026-preflight-gpu-v1 \
  --output rendered/ferrel2026-seedtest-v1
```

The wrapper independently checks the recorded scheduler states; a dependency is
not a substitute for successful completion. A failed reader uses a fresh output
root on retry. No deployment, submission or scientific execution was performed
by the publication agent while preparing this version.
