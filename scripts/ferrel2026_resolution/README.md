# Ferrel horizontal-resolution and heat-transport observer

Allocated Unity CPU analysis only; no model, checkpoint or publisher source is
modified. This is an under-resolution screen, not a convergence certificate.
The observer requires successful Slurm jobs, passed source/configuration gates,
COMPLETE metadata and every expected frame. It refuses nonempty output roots.

Run the bounded analytical checks before real-field interpretation:

```bash
sbatch slurm/ferrel2026-resolution.sbatch --self-test-only \
  --output rendered/ferrel2026-resolution-selfcheck-v1
```

The same checks run automatically at the start of a full observation. They verify
a known cosine mode in an even-reflected channel, its variance, physical
wavelength and inverse wavenumber, positive/negative heat covariance, and a weak
absolute-Float32-temperature signal that must be labeled indeterminate. A
deliberately inconsistent heat profile must be rejected by the derived precision
bound. No analytical checks are intended to run locally.

Example complete pilot observation:

```bash
sbatch --kill-on-invalid-dep=yes --dependency=afterok:1409341 \
  slurm/ferrel2026-resolution.sbatch --stage pilot \
  --completed-job 1409341 --output rendered/ferrel2026-resolution-pilot-v1
```

Example complete production plus independent horizontal refinement:

```bash
sbatch --kill-on-invalid-dep=yes --dependency=afterok:1409345:1409342 \
  slurm/ferrel2026-resolution.sbatch --stage production --include-fine \
  --completed-job 1409341 --completed-job 1409342 \
  --completed-job 1409343 --completed-job 1409344 --completed-job 1409345 \
  --output rendered/ferrel2026-resolution-production-v1
```

Omit `--include-fine` and its job/dependency for an earlier coarse report. Every
relevant simulation job must be listed; model metadata becomes COMPLETE before
post-run validation, so metadata alone is insufficient. Default roots match the
existing Ferrel campaign and can be overridden explicitly.

## Fields and geometry

The observer reads hourly saved `wcenter_slices` and `theta_slices` at the three
actual video heights, plus Float64 `theta`, vertical KE, heat-flux, diffusivity,
temperature, density, pressure and mean-N² profiles. It verifies the exact plane
indices and that casting the first Float64 theta plane reproduces the first
Float32 theta slice. Only that matching first plane uses Float64 temperature.
The other two absolute-temperature planes have explicit Float32 uncertainty.

Full complex FFTs use `norm=forward`, periodic x and an even y reflection to
2Ly. Therefore ky is spaced by π/Ly. No modes are dropped or filtered. Separate
plane-mean-removed and zonal-mean-removed spectra answer different questions;
only the former corresponds to total resolved horizontal-mean heat transport.
Wavelength is 2π/k; inverse wavenumber is 1/k. Both are reported. The peak is the
power-weighted k within the maximum-power radial shell, with shell edges and
candidate-count uncertainty. Tail fractions below 4/8/16 grid spacings are
reported separately for radial wavelengths and directional component wavelengths.
Those thresholds are descriptive screens, not universal resolution criteria.

## Precision accounting

Each saved Float32 value has a half-ULP bound from its adjacent representable
neighbors, converted to Float64 before subtraction. This handles zero,
subnormals and power-of-two binade boundaries. Mean removal is an orthogonal
projection, so the raw RMS error bound also bounds the demeaned field.

For observed demeaned fields W,T, RMS input bounds eW,eT give the Fourier L1
bounds `2*rms(W)*eW + eW²` for variance and
`eW*rms(T) + eT*rms(W) + eW*eT` for heat cospectra. The same bounds cover any
subset or shell. Arithmetic and quantization contributions are stored separately;
their combined envelope includes mixed products. Mean-removal arithmetic uses
the standard `gamma_n = n*eps/(1-n*eps)` summation/subtraction envelope. The
conservative FFT identity allowance uses gamma_(32*N_extended), assuming a
numerically stable Float64 FFT and exceeding a direct complex DFT's operation
depth plus products/reductions. It is an arithmetic-model allowance, not certified
interval arithmetic for arbitrary FFT implementations or twiddle generation.
Observed Parseval and cospectral identity residuals must fit that separate bound;
these identities do not independently establish modewise forward-error bounds.

Profile agreement compares the raw saved-plane variance against Float64 vertical
KE/rho and uses the *stored Float64 theta-profile mean* when comparing heat flux.
This reproduces the model's subtraction reference rather than assuming Julia and
NumPy means agree bitwise. The profile arithmetic envelope covers both independent
N-value reductions and products, with throughput enlarged by quantization error.
Pointwise quantization bounds, reduction error and the absolute-scale subtraction
error are reported separately. No arbitrary relaxed
profile tolerance is substituted for those bounds.

Heat cospectra are labeled indeterminate when anomaly or absolute-cospectral
throughput cannot be separated from the derived precision floor. Tail fractions
have conservative lower/upper intervals; when the denominator can include zero,
the estimate is omitted. A bounded whole spectrum does not certify each mode.
Positive/negative signed covariance and absolute modal throughput are distinct.

## Time windows and interpretation

Every complete case has a 0–6-hour mean; 48-hour cases also have a 12–48-hour mean.
The excluded first12h are a convention, not a demonstrated equilibration time.
Six-hour refinement/controls cannot validate developed48h dynamics. Time averages
use the actual saved times with trapezoidal weights; no temporal spectrum is
computed. Hourly snapshots undersample the reference cap's buoyancy-wave period.

Full-depth aggregate vertical-KE, heat-flux and N² profiles identify activity
outside the sampled planes. The three planes cannot establish bottom-layer or
vertical resolution, and 1024²×96 changes horizontal resolution only. Nominal
Pr/Ek/flux-Rayleigh/flux-Rossby analogues explicitly use the configured neutral
layer depth and lowest-cell density; they are not calibrated RBC criticality
criteria. Sampled actual SSPRK3 step lengths are reported separately from the
configured30s cap, with no assertion about the active limiter.

Outputs are compact scalar CSVs, time-mean shell/profile CSVs, two diagnostic
figures per case, `report.json` and checksummed `manifest.json`. No raw fields or
checkpoints are exported. The parent must review precision flags and figures,
combine these with native budgets and matched controls, and decide the next
experiment from the actual completed results. Existing synchronized movies remain
the required accompaniment to any new scientific result presentation.
