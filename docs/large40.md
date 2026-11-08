# Larger channel with a maintained stable cap

This campaign tests whether bottom-driven convection develops larger vortices and beta-dependent organization when the domain is wide enough and the cap is maintained. It replaces the fixed upper heat extraction that eroded the base of the cap in the 20-day experiments. These remain dry, idealized anelastic calculations, not a calibrated reconstruction of Jupiter or a demonstrated FFR simulation.

## Configuration

| Quantity | Setting |
|---|---|
| Channel | 6,000 × 6,000 × 120 km; periodic x, free-slip impermeable y and z boundaries |
| Grid | 1,024 × 1,024 × 64; Float64 dynamics |
| Spacing | 5.859 km horizontally; 1.875 km vertically |
| Aspect ratio | 50:50:1 for total depth; 75:75:1 for the initially neutral 80-km layer |
| Rotation | Jupiter Ω = 1.7585 × 10⁻⁴ s⁻¹, latitude 55°N; f₀ = 2.88095774 × 10⁻⁴ s⁻¹ |
| Matched cases | f plane, and beta plane with β = 2.8216700143 × 10⁻¹² m⁻¹ s⁻¹ |
| Beta variation | f changes by 5.88% of f₀ across the full channel |
| Bottom heat | 162 W/m², the previous 30× sensitivity case |
| Initial turbulence | Same random seed in both cases, wavelength band around 300 km, eddy Rossby number 0.1 |
| Cap | Initial N² = 2.5 × 10⁻⁵ s⁻² above 80 km; local potential-temperature relaxation described below |
| Duration | 40 Earth days per case, four restart segments of 10 days |
| Fields | Three vorticity planes every 3 Earth hours, plus density-weighted depth-averaged velocity and stresses |
| Monitoring | Hourly bulk quantities, full mean N² profile in saved frames, protected-cap gradient at every diagnostic time |
| Checkpoints | Every 2 Earth days and at segment ends; only the newest checkpoint within each segment is retained |
| Compute | One 80-GB H100, 4 CPU cores and 128 GB host RAM per production job |

The horizontal spacing is 25% coarser than the previous 256², 1,200-km runs. This is a large-domain pilot, not a resolution-converged filament experiment. The new box accommodates roughly one to two wavelengths of the previously estimated 3,300–4,600 km Rhines wavelength; it does not guarantee a large statistical population of jets. The conventional Rhines inverse-wavenumber length and 2π times that length must not be confused.

## Thermal reservoir

The cap tendency is

```
∂θ/∂t = −M(z) [θ − θ_reference(z)] / (1 Earth day),
M = s²(3 − 2s),  s = clamp((z − 80 km)/(10 km), 0, 1).
```

Thus it vanishes below 80 km and reaches full strength above 90 km. The reference profile is the same initial stable profile. There is no simultaneous fixed upper heat sink. The actual upper heat exchange is measured from the evolving state and can be positive or negative; it is not prescribed to balance the bottom flux. The lower-layer heat budget can continue to evolve. This is an idealized thermal reservoir, not a radiative-transfer calculation. Local relaxation also damps cap temperature perturbations and waves.

The initial mean stratification is neutral below 80 km. The existing diffusivity taper, upper vertical-velocity sponge, dry hydrogen–helium constants, and bottom heat-flux conversion are retained. A runtime check stops a run if the horizontally averaged N² anywhere in the protected 90–115 km cap core falls below 20% of its initial value. The 80–90 km transition is allowed to adjust. This is a diagnostic stop, not a reset or clipping of the temperature field.

Breeze's pinned [radiative convection example](https://github.com/NumericalEarth/Breeze.jl/blob/0479d9c8ce006adcb1e3314c324b84980ad4871b/examples/radiative_convection.jl) also uses upper-atmosphere Newtonian relaxation, though with a different reference state, variable and timescale. Our cap parameters are numerical experiment choices.

## Validation and execution

1. Compare the faster, separated horizontal/vertical turbulence construction against the original initializer; check mass-flux divergence, velocity amplitudes and impermeable walls.
2. Compare an uninterrupted integration against a checkpoint continuation, including the three vorticity slices and depth-averaged fields. A checkpointed callback explicitly preserves the adaptive time-step setting, which is omitted by the pinned Simulation checkpoint implementation.
3. Run 40 days in an inexpensive 600 × 600 km, 64² × 64 channel. Require the protected cap core to retain at least 50% of its initial mean N², and check the integrated rho-theta source budget using the measured time-dependent cap exchange. This validates the thermal treatment in a small box, not large-domain morphology.
4. Run a two-hour probe on the actual 1,024² × 64 grid and H100. Check memory, outputs, checkpoints and initial throughput. Require a projected ten-day segment below 8 wall-clock hours, leaving margin within the 12-hour job limit. Mature turbulence may run more slowly.
5. Run the f/beta pair in alternating ten-day segments, with scheduler dependencies preventing continuation after a failed segment. Physics and configuration must agree with the previous segment. The production runtime cap check remains active because a small-box test cannot guarantee large-box stability.

Launch the thermal test with `slurm/large40-thermal.sbatch`, then `slurm/large40-preflight.sbatch` after successful completion. `slurm/large40.sbatch INDEX ROOT PREFLIGHT_ROOT` runs indices 1–8 in order: f10, beta10, f20, beta20, f30, beta30, f40, beta40. Each has a separate output directory and an exact restart from the same case's preceding segment. Existing campaigns use separate files and are unaffected.

The code explicitly refuses CPU fallback for GPU runs and refuses nonempty output directories. Source hashes, configuration hashes, actual GPU name, initial turbulence diagnostics and restart provenance are recorded in metadata. Source-budget checks allow 1% residual relative to absolute source throughput because the integral uses sampled rather than exact discrete tendencies; this is not an exact kinetic-energy budget closure.

## Three-depth movies required for results

Every new results presentation includes actual vorticity videos at all three levels. Requested heights above bottom are **20, 60 and 100 km**; on this grid the nearest centers are **19.6875, 59.0625 and 100.3125 km**. These sample the lower convection, upper convection and stable cap. Saved vorticity slices use Float32 for storage; dynamics, bulk diagnostics and depth-averaged fields remain Float64.

```sh
python scripts/vorticity_depth_videos.py \
  data/large40-v1/large40_channel_f data/large40-v1/large40_channel_beta \
  --labels 'f plane' 'Jupiter beta plane' --output gallery/large40
```

Use `--preview` to label an incomplete but contiguous time series. The renderer requires all saved times, three distinct real planes and synchronized times for paired comparisons. It creates three MP4 files, final-frame posters and a JSON record. Each depth has a fixed symmetric color scale shared across all times and both cases, with no clipping or temporal interpolation. Older runs that saved only one horizontal level cannot supply the other two depths retroactively.

Assess inverse transfer with depth-averaged transfer diagnostics alongside the movies, spectra, correlation lengths and zonal mean flow. Growth of coherent structures or a zonal band alone does not establish an inverse cascade or persistent jets. Cap heat exchange and mean stratification must accompany any physical interpretation.

## Submission record — November 8, 2026

Thermal/restart validation: job **66946536**, running when submitted. Full-size H100 probe: **66946541**, dependent on validation. Production segments: **66946542–66946549**, ordered as in the table below. They are submitted with `afterok` dependencies and cancellation of invalid dependencies; no production segment can start before its predecessors succeed. The validation and full-size probe have now passed; production is waiting on scheduling priority. See the update below.

| End day | f-plane job | beta-plane job |
|---|---|---|
| 10 | 66946542 | 66946543 |
| 20 | 66946544 | 66946545 |
| 30 | 66946546 | 66946547 |
| 40 | 66946548 | 66946549 |

[Machine-readable launch record](large40-launch.json). The initial validation submission was replaced before production to preserve the adaptive time step explicitly across checkpoints. Raw records remain on Unity.

For a fresh campaign, `bash scripts/submit_large.sh NEW_PRODUCTION_ROOT NEW_PREFLIGHT_ROOT` submits the complete dependency chain. This does not imply the scientific outcome is known: failed stability, restart, budget or throughput tests block production.

## Validation update — November 8, 2026

Both gates completed: thermal job 66946536 took 18 minutes 48 seconds including startup and restart tests; full-grid H100 probe 66946541 took 5 minutes 34 seconds. The 40-day small-box test retained at least 99.90% of initial protected-core mean-profile N². See [all three vorticity movies and the heat-budget interpretation](../gallery/large40-validation/README.md).

The H100 probe's measured post-startup rate projects 2.34 wall-clock hours per ten simulated days. All eight production allocations were shortened to 12 hours (about fivefold margin), and the future-launch throughput gate tightened accordingly. This changes the scheduling request, not the physical integrations. At this update the first production segment is still pending on priority; no mature large-domain morphology or jets have yet been measured.
