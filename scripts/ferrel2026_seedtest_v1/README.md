# Six-hour horizontal thermal-seed tests

The completed512 pilot has weak lower-level vertical motion whose spectrum peaks
near two horizontal cells, with more than96% of its variance below four-cell
radial wavelengths. Both heated and unheated cases show this. Conservation checks
passed, but this is not resolved convection. These four new tests isolate the
thermal seed before another forcing amplitude or long extension is considered.
All existing v3 drivers, configurations, gates and jobs remain unchanged.
Watch the published pilot's synchronized three-depth vorticity/tracer movies on
the [GitHub gallery landing page](https://npizzo.github.io/jupiter-ffr-breeze/).

| Case | Grid | Bottom input | Initial thermal perturbation |
|---|---:|---:|---|
| `zero512_b0` |512²×96|0 W/m²|zero|
| `zero512_b162` |512²×96|162 W/m²|zero|
| `modes512_b162` |512²×96|162 W/m²|fixed horizontal analytic modes|
| `modes1024_b162` |1024²×96|162 W/m²|the same physical modes and per-layer variance|

Each case starts fresh, stops at6 Earth hours, and writes7 hourly frames. No
checkpoint pickup or extension is accepted by this driver. The existing velocity
seed, jets, reference state, rotation, viscosity, diffusivity, thermal feedback,
cap treatment, timestep policy and numerical operators are unchanged.

## Seed definition and what it tests

`thermal_seed.jl` contains an explicit eight-mode coefficient/phase specification:
periodic zonal cosines multiplied by meridional cosines satisfying zero normal
derivative at the channel walls. The admitted band is375–750km; actual mode
wavelengths, coefficients and every layer's phase numerator are recorded in each
case's `thermal_seed.toml`. Actual shortest wavelength is about412km. Modes are
evaluated at each grid's physical cell centers. Coarse centers do not coincide
with every second fine center.

The original `perturbation=1e-5K` multiplies uniform noise on[-0.5,0.5], so it was
not an RMS amplitude of1e-5K. On the allocated node, replay the original512²×96
seed42 draw and horizontal demeaning, then compute its realized variance in every
layer. Preserve those exact pre-addition variances and the original vertical
envelope `exp(-(z/(0.35H))^2)` in both modal grids. This is a fixed reference,
independent of the target horizontal resolution. Production reads the exact
hashed CPU-gate table; it does not recompute the reference reductions on another
node. The per-run configuration records that reference path and SHA. The CPU
gate checks successful reuse and rejection of a changed hash. The zero controls intentionally
omit that variance. Adding the reference absolute Float64 temperature causes
ordinary roundoff; construction and actual output checks report separate
arithmetic allowances.

These are **horizontally band-limited**, not vertically smooth perturbations.
Predetermined layer-dependent phases retain short vertical correlations at
unchanged96 levels. This avoids silently imposing a coherent vertical column,
but does not establish vertical or bottom-boundary-layer resolution. One modal
realization is a controlled seed test, not an ensemble or unbiased instability
wavelength measurement.

The zero-noise heated/unheated pair tests whether vertical motion develops
without the original cellwise thermal seed. Any remaining motion can include
nonlinear/beta adjustment, numerical error and their interaction with heating;
zero thermal noise does not identify their individual causes. The two modal
heated runs compare horizontal resolution at the same seed. A matched modal
unheated control would still be required before attributing its later mean-flow
response to heating.

## Validation and source inheritance

Only Unity allocations may execute any scientific initialization, model,
diagnostic or rendering work. Source editing and text/hash/shell-syntax checks are
the only local work performed in preparation.

`driver_inheritance.json` specifies exact textual changes to the frozen v3
driver. `verify_inheritance.py` checks the original hash, every allowlisted change
and byte-identical copied physics components. Configuration validation permits
only the case name, horizontal grid,6h stop time, seed settings and added
precision-output levels to differ from the original matching heating case.
Inherited v3 physics checks are identified separately; old gates do not silently
validate a new source hash.

The CPU initialization gate requires the existing successful v3 `PASS512.toml`.
Both CPU and GPU launchers use the repaired `environments/ferrel2026_gpu` project,
the `depot-ferrel2026-gpu:depot` ordering, and `--startup-file=no`. A new narrowly
scoped launcher calls the frozen repaired-runtime verifier against the prior
`runtime.toml`, checks the active project, and permits only the two new entry
points. The old runtime launcher and environment files remain unchanged.
It verifies periodicity, wall parity, zero mean, the original variance reference,
zero controls and absolute-temperature addition. It samples analytic modes on
64/128 grids without evolving a model; these are mathematical sampling checks,
not physical-resolution claims. It separately checks every one of96 production
planes on both actual512/1024 grids against the common analytic field and RMS.
The gate pins all new sources, configurations and the recorded mode reference.

The GPU campaign requires that new gate plus unchanged v3 `PASS512.toml` and
`PASS1024.toml`. Parent review of the old fine pilot precedes submission. Each new
case reruns the existing native momentum, mass circulation, tracer and thermal
checks, then verifies the new Float64 plane indices, initial seed, native profile
agreement and exact0–6h frame coverage. `CASE_PASS.toml` is written only after
these checks and the existing transfer observer succeed. Publication must also
require the enclosing Slurm job to finish successfully, because the wrapper's
GPU-memory-headroom check runs after the Julia process exits. A `COMPLETE`
metadata flag alone is insufficient.

## Additional output and comparison criteria

The usual synchronized vorticity/tracer movies still use19.375,59.375,99.375km.
Only new runs add Float64 centered-w and horizontally demeaned-theta arrays at
those exact planes and at5.625/9.375km, where the original heated pilot's
heat-transport/vertical-motion profiles peaked. Each frame stores
`precision_wcenter_slices`, `precision_theta_anomaly_slices`, `precision_heights`
and `precision_indices_julia_1based`. Existing Float64 mean-temperature profiles
are preserved. This prevents Float32 absolute-temperature quantization from
erasing weak thermal structure. It does not imply the dynamics itself is
converged. Existing observers must be explicitly extended to use these new
arrays rather than silently continuing to read the old Float32 temperature
slices at the upper levels.

Use native budgets, profiles and spectra from these levels—not movie morphology
alone. At6h compare the peak's physical wavelength and cells per wavelength,
4/8/16-cell variance and absolute heat-throughput tails, signed heat transport,
and the heated/zero-control difference. If the peak follows grid spacing or
transport remains concentrated near a few cells, do not extend these as a
validated convective experiment. If zero-noise motion stays near the grid, test
the adjustment/timestep pathway separately. If the modal response is resolved
but still weak or immature, a later24h test may be justified after review, with
an unheated modal counterpart. No additional24h runs are included here.

## Parent submission instructions

From the Unity repository root, with fresh paths and logs directory available:

```sh
sbatch slurm/ferrel2026-seedtest-v1-unit.sbatch \
  data/ferrel2026-seedtest-preflight-v1 data/ferrel2026-preflight-gpu-v1
```

After that allocated gate passes and the existing fine pilot is reviewed, submit
the three512 cases using `slurm/ferrel2026-seedtest-v1-512.sbatch`, and the fine
case using `slurm/ferrel2026-seedtest-v1-1024.sbatch`. Both accept:

```text
CASE data/ferrel2026-seedtest-gpu-v1 data/ferrel2026-seedtest-preflight-v1 data/ferrel2026-preflight-gpu-v1
```

Use scheduler dependencies on the successful initialization and fine-pilot jobs.
Submit each case once. Do not overwrite or cancel any existing run. No jobs have
been submitted by the publication agent. Actual runtime results remain pending.
