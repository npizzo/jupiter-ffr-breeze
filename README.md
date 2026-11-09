# Jupiter folded filamentary regions with Breeze.jl

**[Watch the movies on the gallery landing page](https://npizzo.github.io/jupiter-ffr-breeze/)**


Regional, bottom-heated, rotating convection experiments to determine whether Breeze.jl can reproduce the dynamics and morphology of Jupiter’s folded filamentary regions (FFRs).

**Status (November 9, afternoon):** Both original f/beta cases completed twenty days, and the heated/heating-off shear pair completed ten days. The new large-eddy validation passed all gates. [Watch the updated synchronized comparisons and validation movies](https://npizzo.github.io/jupiter-ffr-breeze/).

[Watch gallery](https://npizzo.github.io/jupiter-ffr-breeze/) · [Latest campaign](docs/large40.md) · [Initial validation](docs/results.md) · [Scientific plan](docs/science.md) · [Unity instructions](docs/unity.md)

## Initial configuration

Dry anelastic hydrogen–helium atmosphere on a periodic regional f-plane at 55°N; 1,200 × 1,200 × 120 km; lower weakly unstable layer, stable cap above 80 km, bottom enthalpy flux 5.4 W/m², relaxation aloft, and a vertical-velocity sponge in the upper 10%. Free-slip, impermeable vertical boundaries. Gravity and thermodynamics are explicit configuration parameters. Potential temperature is referenced to 1 bar; bottom pressure is 20 bar.

The weak initial instability and random thermal perturbations accelerate the initial convection test: its energy supply is not exclusively the imposed bottom flux. Radiative relaxation is idealized. The experiment is not a reconstruction of a measured Jovian sounding. The initial 1,200-km box tests convection and smaller-scale organization, not a complete several-thousand-kilometre FFR.

## Reproduce

Julia 1.12 is required. Breeze is pinned to commit `0479d9c8ce006adcb1e3314c324b84980ad4871b`; the generated Manifest records all resolved dependencies.

```sh
julia --project=. scripts/setup.jl
julia --project=. scripts/run.jl configs/smoke.toml data/smoke CPU
julia --project=. scripts/run.jl configs/pilot.toml data/pilot GPU
python -m pip install -r requirements.txt
python scripts/gallery.py data/pilot
```

The GPU option fails if CUDA is unavailable rather than silently running on CPU. Output directories must be empty. Raw output, credentials, local cluster paths, and environments are excluded from Git. Figures and sanitized parameters belong in `gallery/`.

## Cases

| Case | Grid | Simulated duration | Purpose |
|---|---|---|---|
| rest / heating | 16 × 16 × 24 | 60 s | hydrostatic rest and imposed heat-flux budget |
| smoke | 16 × 16 × 24 | 60 s | construction, stepping, and output check |
| pilot | 128 × 128 × 64 | 4 h | initial convection and measured resource use |
| regional | 256 × 256 × 96 | 2 Earth days | optional longer experiment; not submitted automatically |

`jet_speed=0` by default. A nonzero value adds a sinusoidal initial jet and a prescribed meridional body force balancing its Coriolis acceleration; it is a separate forced-shear experiment, not a self-generated jet. This option requires its own validation before scientific interpretation.

## Provenance

Built with [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl) and [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). Upstream libraries retain their licenses. See the scientific plan for references and limitations.

## Controlled convection experiment

Four 48-hour experiments compare bottom heating on/off at two initial temperature-perturbation amplitudes (0.2 and 0.02 K; actual seeded noise is smaller than these amplitude parameters). The grid, reference stratification, and seed are identical. These runs isolate dependence on the imposed bottom heat flux from the initial instability. See [experiment design](docs/convection48.md).

The initial convection experiments used a **traditional f-plane**: Ω = 1.7585 × 10⁻⁴ s⁻¹, latitude 55°N, f = 2.88096 × 10⁻⁴ s⁻¹, beta = 0. Jupiter’s rotation rate stays fixed. Rotational dominance must be diagnosed from the flow; the initial convective Rossby number is about 1.1.

## Revised forcing experiment

The [completed revised experiment](docs/controlled48.md) replaces reference-profile relaxation with a prescribed upper heat sink and reduces thermal mixing across the stable cap. It repeats the heating controls and adds an unforced case and a vertical-resolution check. The model and geometry remain dry anelastic on the same Jovian f-plane.

[Revised gallery](gallery/controlled48/README.md) · [Measured results and next diagnostic](docs/controlled48-results.md)

## Neutral initial-layer experiment

The [neutral-layer test](docs/neutral48.md) removes the initial mean instability and reduces seeded temperature noise to 10⁻⁵ K. It compares heated/cooled, cooled-only, and unforced cases, with a vertically refined heated case. All four runs completed; [snapshots at 6, 12, 24, and 48 hours](gallery/neutral48/README.md) and [results](docs/neutral48-results.md) are available.

## Seeded turbulence and convection

The [new comparison](docs/turbulence48.md) starts all three thermal-forcing cases from identical three-dimensional, band-limited velocity fields with initial eddy Rossby number 0.1. It tests how convection modifies decaying turbulence. All checks and three runs completed. See [results](docs/turbulence48-results.md) and [snapshots](gallery/turbulence48/README.md).

## Enhanced forcing campaign

The [stronger-convection comparison](docs/strong48.md) raises heating to 30× and 100× baseline, with matching cooling-only controls and a horizontally refined 100× pair. All six cases completed; see [results](docs/strong48-results.md) and [snapshots](gallery/strong48/README.md). These fluxes are idealized sensitivity choices, not inferred Jovian heat fluxes.

## Five-day 30× experiment

The [longer 30× campaign](docs/convection120.md) completed: matched heated/cooling-only pairs on 256 × 256 × 64 and 256 × 256 × 128 grids, for 120 Earth hours. See [results](docs/convection120-results.md) and [vorticity videos and snapshots](gallery/convection120/README.md).

## Longer runs and beta

The [20-day campaign](docs/regional20.md) adds a physical Jupiter beta plane, matched f-plane channels, a longer periodic reference, and short 512² horizontal-resolution comparisons. New depth-averaged transfer diagnostics distinguish upscale energy transfer from morphology. Production is gated on GPU preflight checks.

## Maintained cap and larger channel

The [large-domain campaign](docs/large40.md) uses matched 6,000-km f/beta channels, a 1,024² × 64 grid and local cap relaxation. Every new results entry must include vorticity videos near 20, 60 and 100 km above the model bottom.

## Shear, tracer, thermal feedback and moisture

[Shear/tracer campaign](docs/shear20.md): validation and full-grid checks passed; both heated and heating-off ten-day segments completed. No-shear and horizontally refined cases remain queued. A separate thermal-feedback control addresses net heat accumulation.

[Moist-convection validation](docs/moist-validation.md): thermodynamic checks, a resting atmosphere, and dry/moist plumes at two resolutions completed. [Watch the synchronized three-depth movies](https://npizzo.github.io/jupiter-ffr-breeze/#moist-validation). Local plume extremes are still resolution-sensitive; moisture is not yet enabled in the FFR production campaign.


New campaign: [larger random eddies on f and beta planes](docs/eddy20.md), four matched heating-on/off cases over twenty Earth days. [Queue record](docs/eddy20-launch.json). Initializer, short-run and production-grid validation passed; production is awaiting scheduling. Completed validation movies are included in the gallery.
