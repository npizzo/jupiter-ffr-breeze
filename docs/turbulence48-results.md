# Initial turbulence with Rossby 0.1: completed comparison

Three matched 48-hour runs completed with identical initial three-dimensional random velocity fields. At the retained physical heating rate, bottom heating changes second-day kinetic energy by **+0.117%** and vertical RMS speed by **+0.086%**, relative to the cooling-only control. The bulk evolution is still dominated by decay and adjustment of the seeded flow. This is a measurable difference in one paired realization, not a statistically established or large convective reorganization of turbulence.

[Snapshots, spectra and data](../gallery/turbulence48/README.md) · [Design](turbulence48.md) · [Summary](../gallery/turbulence48/summary.json)

## Initial condition and checks

The reference eddy wavelength is 300 km, with Fourier wavenumbers within 25% of its central wavenumber. These wavelengths span about 26–43 horizontal cells. Rossby normalization uses L = λ/(2π) = 47.75 km, lower-layer horizontal RMS velocity **1.3756 m/s**, and Jovian f = 2.88096 × 10⁻⁴ s⁻¹. Thus U/(fL) = 0.1. The lower-layer vertical RMS speed is **0.2751 m/s**. With viscosity 1,000 m²/s, UL/ν is approximately **65.7**. This is a modest-Reynolds-number experiment, not a high-Reynolds-number reconstruction of Jupiter.

Toroidal and poloidal potentials generate a discrete divergence-free mass flux consistent with the anelastic density profile. Unit checks passed for amplitude, boundary flow, continuity, repeated-seed identity and sensitivity to changing the seed. A short GPU test verified the actual solver output preserved the prescribed amplitudes. All long runs passed the same source-budget and finite-output checks. The analysis verified exact equality of the saved initial velocity maps and energy profiles across the three cases, together with matching physical parameters and source-code hashes. Initial relative mass-flux divergence was **2.9e-16**.

The initial mean thermal profile is neutral below 80 km, and the tiny temperature perturbation is the same in every case. The stable cap, prescribed diffusivity taper, box, anelastic model, Jovian rotation, beta = 0 and zero imposed shear are unchanged. The random field is not an equilibrated turbulent state; early evolution includes adjustment and wave motion. There is no continuous mechanical forcing, and no nonrotating control in this campaign.

## Comparison

Second-day averages cover hours 24–48. Column kinetic energy includes the model's face-weighted horizontal and vertical velocities. RMS speeds in the table are volume averages over the domain at cell centers, unlike the lower-layer normalization values above.

| Case | Second-day column KE [J/m²] | Second-day vertical RMS speed [m/s] | Final KE / initial KE | Final maximum resolved upward heat flux [W/m²] |
|---|---:|---:|---:|---:|
| Heated + cooled | 19558.57 | 0.121994 | 27.049% | 0.260872 |
| Cooling only | 19535.72 | 0.121890 | 26.985% | 0.256821 |
| Unforced | 19535.56 | 0.121888 | 26.991% | 0.257682 |

All runs start with **60077.9 J/m²** of column kinetic energy. Bottom heating changes final kinetic energy by **+0.235%** relative to the cooling-only control. The two controls differ by **+0.0008%** in second-day energy and **+0.0016%** in second-day vertical RMS speed, indicating how little the prescribed upper cooling affects these bulk quantities in this window.

The heated case's final maximum resolved upward heat flux is **0.2609 W/m²**, versus **0.2568 W/m²** in the cooling-only case. The latter is important: a nonzero heat flux can result from the seeded flow redistributing the existing stratification, even without bottom heating. The difference between these maxima is **0.00405 W/m²**; maxima need not occur at the same height. This diagnostic is not an exact boundary flux or complete heat budget.

The final heated and cooling-only vorticity maps near 20 km have correlation **0.999054**. Difference images use separate, smaller color scales and must not be read as the size of the change relative to the total flow. A single pair cannot separate a robust change in turbulent statistics from trajectory sensitivity.

## Spatial structure and snapshots

The gallery shows initial and evolving horizontal vorticity, final common-scale comparisons, heated-minus-control differences, and vertical sections at 6, 12, 24 and 48 hours. Horizontal maps are near 20 km, inside the initially turbulent layer. The time sequences show the evolution of a well-resolved imposed field; they do not establish FFR morphology.

Radial horizontal energy spectra use the saved near-20-km slice, with a checked Parseval identity. They are not full-depth spectra or spectral transfer budgets. The reported spectrum-derived length scales and low-mode energy fractions help compare cases; growth at low wavenumber alone does not demonstrate an inverse cascade. Vertical energy profiles are also shown. The separate neutral-layer tests identified short convection-driven wavelengths close to the mesh scale, and seeding larger eddies does not remove that concern.

## Interpretation and cost

This comparison implements the intended question: how convection changes an existing rotating flow. At 5.4 W/m² over 48 hours, the effect on the chosen bulk turbulence measures is small. Larger motions than in the neutral tests chiefly reflect the imposed initial kinetic energy, not stronger convection. Before claiming appreciable turbulent modification, useful follow-ups would be longer matched integrations, additional matched seeds and a horizontal-resolution check. Stronger heating would be a distinct nondimensional sensitivity experiment, not the same physical forcing.

Measured solver times were 3.24, 3.09 and 3.07 minutes. The [resource record](../gallery/turbulence48/resources.json) gives total job time including startup and validation. No additional runs, beta, imposed shear or model switch were submitted after this comparison.
