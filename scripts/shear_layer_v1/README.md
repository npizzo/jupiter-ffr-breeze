# Shear versus upper-layer anticyclone production — Breeze pilot

Six matched dry anelastic experiments distinguish ambient-shear selection from an upper-layer response to convection. This is a new, isolated campaign; it does not modify the original LSV or the minimal polynomial campaigns. Scientific computation and rendering run on Unity only.

## Design fixed before launch

| Upper-layer target | Cyclonic initial shear | No initial shear | Anticyclonic initial shear |
|---|---|---|---|
| Existing stable cap, N² = 2.5e-5 s⁻² | stable_cyclonic | stable_none | stable_anticyclonic |
| Neutral target, N² = 0 | neutral_cyclonic | neutral_none | neutral_anticyclonic |

All cases retain positive f = 2 Ω sin(55°), Ω = 1.7585e-4 s⁻¹, and beta = 0. Positive relative vorticity is cyclonic. This is a northern-hemisphere mechanism test, not a calibrated Great Red Spot simulation. The neutral case retains the upper relaxation region, its one-day relaxation time, and the same tapered diffusivity; it removes the target stratification, not the thermal controller. Each hydrostatic reference is built consistently with its target, so density and thermodynamic response above 80 km also differ. This comparison is not an intervention solely on one term of the vorticity equation.

- Grid: 128 × 128 × 64, 1500 × 1500 × 120 km; 11.71875 km horizontal spacing, 1.875 km vertical spacing. Horizontal-to-convecting-depth aspect ratio 18.75. This is a smaller-domain pilot, not a continuation of the 6000 km production experiment.
- Boundaries: periodic east–west, free-slip impermeable north–south and top/bottom; no horizontal momentum drag. Existing vertical-velocity sponge remains above 108 km.
- Uniform bottom heat flux 162 W m⁻² (the previous enhanced-heating value), existing enthalpy controller; dry hydrogen–helium thermodynamic constants, no moisture.
- Laplacian viscosity 1000 m² s⁻¹. Thermal diffusivity tapers from 1000 to 1 m² s⁻¹ across 60–80 km, identically in all cases. Tracer diffusivity 100 m² s⁻¹.
- Same full-depth random horizontal velocity at every height, nominal wavelength 400 km and Rossby number 0.1 using L = wavelength/(2π). No initial w. Only odd zonal Fourier indices: translation by half the zonal period reverses the random velocity, giving zero initial random-vorticity skewness. Same thermal perturbations and tracer patches across cases.
- Initial zonal flow U(y) = -S Ly/pi sin(pi y/Ly), S/f = +0.05, 0, or -0.05. Discrete geostrophic initial pressure; no continued shear forcing. Opposite shear cases have equal initial background energy. The initial zonal-flow vorticity is one-signed and entails wall circulation, explicitly tracked; it is not counted as a vortex.
- Ten Earth days per pilot, no automatic long extension. One realization cannot establish a robust preferred sign.

## Decisive comparisons

Let zeta_prime = zeta - mean_x(zeta), subtracting the instantaneous zonal mean separately at every depth. Report both the total field and fluctuations; skewness of total vorticity alone is not evidence for shear selection.

1. **Shear:** at fixed cap, does changing shear sign change coherent-vortex sign and nonzonal signed skewness? Track actual shear throughout. If shear disappears, a null outcome is not a maintained-shear test.
2. **Upper layer:** without shear, does a stable cap produce an anticyclonic upper-layer response relative to the neutral-target control while the lower/full-column flow remains cyclonic? Negative upper skewness alone is insufficient; require finite coherent nonzonal energy and persistent structures.
3. **Interaction:** does shear amplify, suppress, or reverse the cap response? A factorial difference of cap effects across opposite shears separates an interaction from the individual contrasts.
4. **Scale:** measure absolute low-mode nonzonal energy and spatial structure. Do not call a vortex an LSV from skewness alone or claim an inverse cascade without an energy-transfer calculation.

Layer averages are density weighted: lower 10–70 km, interface 75–90 km, upper 90–105 km (below the sponge), and the entire column. Include full domain and central 80% meridional statistics. Save profiles, actual zonal flow, resolved divergence/stretching/tilting, three movie planes near 20/60/100 km, and daily full staggered fields. The stretching and tilting estimates are centered diagnostic correlations, NOT a closed native WENO budget. They cannot alone establish the production mechanism. No silently discarded boundary or numerical remainder.

## Validation and limits

Preflight runs all six initializations on the actual grid; verifies sign-neutral random fields, matched velocities, shear signs, finite output, tracer conservation, thermal balance, stable-cap preservation, and restart agreement. It measures throughput including diagnostics; projected ten-day cost must be below 1800 seconds, leaving headroom below the internal 3000-second stop. Production is an array of six jobs, maximum one running at once, with a hard one-hour Slurm limit each. The validation allowance is 45 minutes; any wall-limit stop is explicit and does not authorize another run.

Resolution is a screening limitation, especially for convective plume widths and upper-layer penetration. A null result at ten days is not proof that either mechanism is absent. Any extension or refined comparison should be selected after examining convection strength, resolved scales, persistent shear, cap state, and signal relative to initial/transient variability. Reversing f is not part of this experiment.

## Literature motivation

- Heimpel, Gastine & Wicht (2016), *Simulation of deep-seated zonal jets and shallow vortices in gas giant atmospheres*: https://doi.org/10.1038/ngeo2601 — stable outer layer and convection-driven anticyclones, with shear-dependent large vortices.
- Marcus (1988), *Numerical simulation of Jupiter's Great Red Spot*: https://doi.org/10.1038/331693a0 — coherent-vortex survival tied to ambient shear sign.
- Yadav, Heimpel & Bloxham (2020), *Deep convection-driven vortex formation on Jupiter and Saturn*: https://doi.org/10.1126/sciadv.abb9298 — distinguishes organization in shear from plume-divergence-driven anticyclones. These global simulations are motivation, not a quantitative benchmark for this Cartesian pilot.
