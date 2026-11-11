# Polar gamma-plane reproduction pilot

**Superseded before science runs began.** The CPU gate failed an overly strict
bitwise FFT equality check and all dependent jobs were cancelled automatically.
Use the [SI-informed v2 campaign](gamma2026_v2.md) and
[launch record](new-campaigns-launch.json). The description below preserves the
provisional design; it is not the configuration of the submitted v2 science runs.

This isolated first batch reconstructs the main-text barotropic QG setup of
[Siegelman, Young and Ingersoll (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170063/).
It is **not yet an exact numerical reproduction**: the corrected supplementary
initial taper, numerical grid and filter settings have not been independently
verified. No convection is present in this baseline.

The model evolves relative vorticity with the published polar or flat PV trap:

```
q = zeta + eta,    Dq/Dt = 0
eta_polar = -gamma*r^2/2 inside r_trap; 0 outside
eta_flat  = -gamma*r_trap^2/2 inside r_trap; 0 outside
```

The solver is official GeophysicalFlows.jl SingleLayerQG, pinned to commit
`51ca98a0d471ce4ef8bbba8b73daf2bad54d53fb`, in its own Julia environment and
`depot-gamma2026`. It does not modify the running Breeze environment.

## Initial batch

Six cases: polar and flat geometry for seeds 11, 22 and 33; 2048² dynamics,
75 Earth days per case. This tests early evolution and segregation, **not mature
crystal formation**, which occurred around 1.5 years in the published reference.
Reference parameters: gamma=7.8690e-20 m^-2 s^-1, domain RMS speed 80 m/s,
Lgamma=(U/gamma)^(1/3), domain width12 Lgamma, trap radius5 Lgamma,
initial inverse-wavenumber scale200km (wavelength about1257km).
The published November2022 correction to the 2pi convention is respected.

Explicit reconstructed numerical choices:

- Real Gaussian random noise projected onto a spectral annulus of half-width
  2.5% around k=1/(200km), then spatially tapered.
- The taper is unity to2.5 Lgamma, cosine-ramped to zero at3 Lgamma.
- Remove the mean vorticity and normalize the domain RMS speed after tapering.
- Float64, FilteredRK4, strict two-thirds dealiasing, adaptive advective CFL0.25,
  maximum60-second timestep, no explicit viscosity or linear drag.
- Same initial flow for each polar/flat seed pair.

The 4096² comparison and longer integrations are deferred until measured
throughput, numerical losses, initial-band fidelity and this first batch are
reviewed. Neither a successful gate nor a visually organized pattern proves
convergence, a particular vortex count, or reproduction of Jupiter.

## Validation and resource bounds

The setup job installs dependencies only. A CPU job checks an analytic steady
Fourier mode, energy normalization, finite polar-trap evolution and exact restart
continuation. A GPU job repeats those checks, compares CPU and GPU solutions,
then measures100 steps at the actual2048² grid.

Production requires the passing gate with the matching core-source hash and
grid. Estimated75-day cost must be below80% of the4.5-hour wall budget.
Every production job requests one L40S,4CPU cores,32GB RAM and6hours.
The pilot fails with a resumable checkpoint if energy exceeds102% or falls
below70% of its initial value, or its4.5-hour budget is reached. These broad
thresholds detect numerical failure; actual energy loss must still be compared
with the paper and a refined calculation.

The six cases execute sequentially with successful-predecessor dependencies.
No old simulation is changed or cancelled. Submission is idempotent through an
exclusive submission record:

```sh
bash scripts/gamma2026/submit.sh data/gamma2026-v1 SETUP_JOB_ID
```

## Output and interpretation

Full spectral checkpoints remain on Unity. Daily output stores512² sampled
vorticity and total PV maps for publication, with sampling stride and physical
coordinates. Dynamics use the full2048² grid. Estimated output is below2GB for
the six-case batch. Raw fields and checkpoints are not published.

This is genuinely2D barotropic dynamics: it has no three physical depth samples.
Its gallery must explicitly say so, rather than fabricating three-depth copies.
The parallel Breeze convection experiments retain actual three-depth
vorticity/tracer movies. Gamma movies do not represent cloud radiance.

Required comparisons before a scientific conclusion are energy/enstrophy loss,
radial cyclone/anticyclone segregation, PV profiles, crystal radius, vortex
counts/spacing and merger rates across seeds; refined-grid and trap sensitivity
remain necessary. A75-day absence of a mature crystal is not evidence against
the published mechanism.
