# Long, coarse f-plane experiment: convection and upscale transfer

The aim is to test how bottom-driven convection modifies upscale transfer in initially large, rotating eddies. This replaces morphology-only inference with measured transfers. It is a screening experiment, not a converged simulation of Jovian folded filaments.

| Setting | Both cases |
|---|---|
| Domain | 6,000 × 6,000 × 120 km, periodic zonal direction, free-slip meridional walls |
| Grid | 256 × 256 × 64; 23.4375-km horizontal and 1.875-km vertical spacing |
| Rotation | f plane, f = 2.88095774 × 10⁻⁴ s⁻¹ (Jupiter at 55°N); beta = 0 |
| Initial random flow | Full-depth coherent, wavelengths 960–1600 km centered on 1200 km; seed 31415 |
| Initial random speed | 5.5022 m/s RMS; Ro = U/(f/k₀) = 0.1, k₀ = 2π/(1200 km) |
| Imposed shear/jet | None |
| Heating | Matched 162 W/m² and zero-bottom-flux cases |
| Thermodynamics | Existing dry anelastic hydrogen–helium model; stable cap; existing enthalpy-feedback controller |
| Duration | 100 Earth days per case; ten checkpoint segments of ten days |
| Saved times | Three hours, including real vorticity/tracer planes near 20, 60 and 100 km |
| Allocation | One L40S, 4 CPUs, 32 GB RAM, six-hour limit per segment; measured throughput gate before production |

The large random seeds are resolved by roughly 41–68 horizontal cells per wavelength. The same physical initialization is used in both cases. The reference time (1/k₀)/U is about 0.40 Earth days; a wavelength crossing time is about 2.52 days. One hundred days is a deliberately longer window, not a prediction of when an inverse cascade or condensate must appear. The box contains only about five initial central wavelengths, so its range of accessible larger scales is limited.

There are one-sixteenth as many cells as in the 1024² × 64 experiments, but no factor-of-sixteen runtime improvement is promised: pressure solves, timesteps, startup, diagnostics and hardware differ. Vertical resolution and prescribed diffusivities are retained. Horizontal convection and numerical dissipation are less well resolved; a positive or negative cascade result needs follow-up at higher resolution. Strong forcing does not by itself imply a convective flow with small Rossby number.

## What is actually compared

The heating-off run is decaying random turbulence in the same stable-cap/controller model, not a statistically stationary mechanically forced turbulence reference. There is no sustained random momentum forcing in either case. Thermal relaxation, numerical/explicit viscosity and the sponge remain present. The two runs isolate the imposed bottom heat input conditional on those model choices; nonlinear trajectories will diverge.

The earlier four-case `eddy20` campaign retains imposed initial cyclonic shear and only a twenty-day target. Those jobs remain unchanged. This new pair removes that shear confound and extends time. Beta is deliberately deferred until the f-plane pair establishes a useful baseline; there are no new beta jobs in this campaign.

## Transfer diagnostics

Each segment computes diagnostics on Unity from saved density-weighted depth-mean velocities (`ubar`, `vbar`) and depth-dependent horizontal Reynolds stresses (`Rxx`, `Rxy`, `Ryy`). No fields are downloaded for analysis or rendering.

1. Reflect the free-slip channel fields with their correct even/odd parity to a periodic domain of meridional length 2 Ly.
2. Fourier transform, retain the strict two-thirds band, and project the interpolated depth-mean velocity onto its solenoidal component. Record the energy removed by truncation and projection.
3. Calculate dealiased nonlinear advection of that resolved mean flow. Bin its energy transfer in horizontal wavenumber and accumulate the flux: **negative flux means upscale redistribution**.
4. Separately calculate work on the mean flow by minus the divergence of the depth-dependent Reynolds stresses. Positive values indicate energy entering the depth-mean flow. This is an exchange term, not a conservative scale flux; not every depth-dependent motion is convective.
5. Track energy at wavelengths larger than the seeded band (greater than 1600 km), energy-weighted inverse wavenumber, zonal energy fraction, and time-dependent spectra. Compare these with vertical motion, cap stability and heat budgets in heated/control runs.

This is a continuum spectral diagnostic of the projected resolved flow, **not an exact discrete WENO energy budget**. It does not close numerical dissipation, sponge, Coriolis transfer or all discretization residuals. A large projection loss, loss of resolved energy, or growth only near the grid cutoff limits interpretation. A negative mean flux over a finite interval supports upscale transfer; it does not by itself establish a broad inertial-range cascade, nor identify every contributing interaction as local in scale.

Rotating convection can feed depth-invariant structures through depth-dependent Reynolds stresses, so that route must be distinguished from redistribution within the depth-mean flow: [Guervilly, Hughes & Jones (2014)](https://arxiv.org/abs/1403.7442), [Maffei et al. (2021)](https://arxiv.org/abs/2003.01669).

## Gates and reproducibility

Source is isolated in `scripts/cascade100`; running/queued older campaigns are not edited. Before production:

- Verify initial speed/Rossby number, discrete divergence, walls, geostrophic pressure estimate, vertical coherence and identical initial velocities in the pair.
- Test the spectral diagnostic against known-mode energy, Parseval, conservative advective transfer, velocity-reversal and scaling identities, and stress-work integration by parts.
- Compare uninterrupted/checkpointed short runs; run both actual 256² × 64 configurations for one day; check finite outputs, tracer conservation, protected cap, thermal accounting, matching initial states and transfer outputs.
- Require estimated ten-day runtime below four wall hours, within six-hour production requests. Freeze/check source and configuration hashes after validation.

Jobs alternate heating on/off to day 10, then day 20, through day 100. Every job requires its predecessor to succeed. Production also requires the saved passing validation record. Each successful segment computes transfer outputs on Unity. No new results are claimed by submitting these jobs.

```sh
bash scripts/cascade100/submit.sh data/cascade100-v1 data/cascade100-preflight-v1
```

Movies will continue to show three depths simultaneously, vorticity above tracer, with fixed labeled scales. Rendering and diagnostics stay on Unity. A depth-mean-flow movie and spectra can supplement those required three-depth views so grid-scale vorticity does not obscure the evolving large-scale flow.

## Submission — November 9, 2026

Unit gate **66992547**, GPU restart/one-day/throughput gate **66992548**, followed by twenty production segments **66992549–66992568**. The order alternates heating on/off at day 10, then day 20, up to day 100. The machine-readable [launch record](cascade100-launch.json) contains every dependency. Production is submitted but still requires both validation gates; do not read submission as a passing GPU validation result.

Ten initial-condition checks passed on Unity. The transfer self-tests also passed on Unity after correcting the normalization for analytically zero transfer. The earlier development test failure and repeat job are retained in the launch record. No existing production job was cancelled or modified.

The submitted unit gate 66992547 subsequently **passed** (14 seconds elapsed), including all ten initialization checks and the spectral-transfer tests. GPU gate 66992548 is pending scheduling priority; all twenty production segments remain dependent on its success.


## Fifty-day publication — November 10, 2026

The validation gates passed and both matched cases completed days 0–50, including all five segment diagnostics. Heated segment five was job 66992557 and unheated segment five was 66992558. Dependent rendering job 67036876 completed on Unity.

[Watch the comparison and separate-case movies](https://npizzo.github.io/jupiter-ffr-breeze/#cascade50). Each contains all 401 real saved times, three depth columns, vorticity above passive tracer, and shared fixed color scales. [Publication records](../gallery/cascade100-first50/README.md) include scale definitions, source/configuration hashes and segment metadata. Scientific computation, rendering and media validation stayed on Unity.

This publishes the first fifty days of the hundred-day campaign; it does not report later queue progress or infer inverse transfer from appearance. Earlier submission and queue observations above are retained as dated history.


## Full 100-day publication — November 11, 2026

Both matched cases completed all ten segments, through day 100. Final jobs 66992567 and 66992568 completed successfully, and all segment transfer diagnostics completed. [Watch the complete comparison and separate movies](https://npizzo.github.io/jupiter-ffr-breeze/#cascade100). Each contains all 801 three-hour outputs at three actual depths. [Publication and encoding provenance](../gallery/cascade100-full100/README.md). Earlier queue observations above are dated history, not current status.
