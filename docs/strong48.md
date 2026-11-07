# Stronger convection interacting with seeded turbulence

The physical-flux comparison produced only a small bulk response to heating over 48 hours. This authorized sensitivity campaign increases bottom heating to **162 and 540 W/m²**, respectively **30 and 100 times** the 5.4 W/m² baseline. These are enhanced forcing levels, not estimates of Jupiter's actual flux.

The initial random velocity field remains normalized to eddy Rossby 0.1, with the same 300-km reference wavelength and seed. The initially neutral lower layer, stable cap, anelastic equations, Jovian f, viscosity and diffusivity taper remain unchanged. Beta and imposed shear remain zero. Each level has an upper sink of equal prescribed power and a matching cooling-only control, which isolates the effect of the bottom heat input from the strengthened upper cooling.

| Cases | Bottom input [W/m²] | Upper input [W/m²] | Grid |
|---|---:|---:|---|
| f030 heated / cooled | 162 / 0 | −162 | 128 × 128 × 64 |
| f100 heated / cooled | 540 / 0 | −540 | 128 × 128 × 64 |
| f100 heated / cooled, h256 | 540 / 0 | −540 | 256 × 256 × 64 |

Every case runs for 48 Earth hours. The higher-resolution pair halves the horizontal spacing from 9.375 to 4.6875 km. The analytical initial velocity modes and phases are shared across grids, but their discrete derivatives and RMS normalization differ slightly; the tiny thermal noise also has a different spatial realization on the changed grid. Within each pair the initial state is identical. This is a sensitivity check, not a claim of convergence.

## Literature and parameter conventions

[Aurnou, Horn & Julien (2020)](https://arxiv.org/abs/2009.03447), equations 8–12, motivate recording Pr = ν/κ, Ek = ν/(fH²), flux Rayleigh Ra_F = B₀H⁴/(νκ²), and modified flux Rayleigh Ra_F* = B₀/(f³H²). Here the reference H is 80 km, κ is the lower-layer value, and the imposed bottom buoyancy-flux scale is B₀ = gF/(ρ₀cₚT₀), evaluated in the hydrostatic bottom reference state. These are explicit comparison conventions for a stratified, capped atmosphere; the Boussinesq paper's onset thresholds and asymptotic scalings do not directly validate this setup.

| Heating factor | Reference Ra_F | Ra_F* | Flux Rossby (Ra_F*)^(1/3) |
|---|---:|---:|---:|
| 1 | 7.69 × 10⁵ | 1.23 × 10⁻⁴ | 0.0497 |
| 30 | 2.31 × 10⁷ | 3.68 × 10⁻³ | 0.1544 |
| 100 | 7.69 × 10⁷ | 1.23 × 10⁻² | 0.2306 |

All cases have reference Ek ≈ 5.42 × 10⁻⁴ and lower-layer Pr = 1. Exact values and conventions are saved in [the parameter record](../configs/strong48-parameters.json). The flux Rossby parameter differs from the initial velocity-based eddy Rossby 0.1 and from a temperature-based convective Rossby number inferred from the evolved stratification. Increasing forcing at fixed f can weaken rotational dominance; measured velocity, vorticity and spatial scales must establish the resulting regime. No asymptotically rapid-rotation or geostrophic-turbulence regime is presumed.

[Bouillaut et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8612362/) show that the distribution of heating and cooling affects transport and that boundary layers can restrict convection. Their internally heated experiment informs interpretation, not a change to this campaign's bottom-boundary forcing. We retain source geometry to isolate forcing strength, while monitoring cap changes and resolution sensitivity. A future distributed-heating test would be a separate experiment.

## Validation and output

Existing continuity, initial-amplitude, finite-field and thermal-inventory checks apply to every run. The cooling-only partners are essential because stronger upper cooling itself can alter the cap. Compare kinetic-energy evolution, vertical speed, buoyancy work, heat transport, vorticity, spatial spectra, and cap preservation against the controls and the completed baseline. Snapshots at 6, 12, 24 and 48 hours use common scales for matched pairs; difference images retain physical units.

Submit `sbatch slurm/strong48.sbatch`. The six cases run sequentially on one L40S, with four CPUs, 16 GB host memory, and a two-hour time limit. The larger host-memory request covers the 256² horizontal grids. Output is isolated in `data/strong48-v1`; completed baseline results are preserved.

Analysis entry point: `python scripts/analyze_strong.py data/strong48-v1 --output gallery/strong48`. The campaign must complete before the full comparison runs. It produces paired metrics, spectra, common-scale final sections and vorticity maps, and time-sequence snapshots.
