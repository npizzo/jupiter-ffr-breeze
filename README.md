# Jupiter folded filamentary regions with Breeze.jl

Regional, bottom-heated, rotating convection experiments to determine whether Breeze.jl can reproduce the dynamics and morphology of Jupiter’s folded filamentary regions (FFRs).

**Status:** CPU/GPU validation and four paired 48-hour convection runs are complete. The controls reveal initial-transient convection, unintended upper-layer heating, and strong diffusion of the stable cap. The next step is to correct those forcing/transport choices; sustained FFR dynamics are not yet established.

[Output gallery](gallery/README.md) · [Convection results](docs/convection48-results.md) · [Initial validation](docs/results.md) · [Scientific plan](docs/science.md) · [Unity instructions](docs/unity.md)

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

All current runs use a **traditional f-plane**, not a beta plane: Ω = 1.7585 × 10⁻⁴ s⁻¹, latitude 55°N, f = 2.88096 × 10⁻⁴ s⁻¹, beta = 0. Jupiter’s rotation rate stays fixed. Rotational dominance must be diagnosed from the flow; the initial convective Rossby number is about 1.1.
