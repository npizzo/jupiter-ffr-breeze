# Convection-to-anticyclone benchmark: first implementation batch

This is a **preliminary Breeze adaptation**, not a reproduced result and not yet a Great Red Spot model. Its purpose is to recover a published cyclone / cyclone-plus-anticyclone transition before attempting another Jupiter cap/jet experiment.

Reference: Cai, Chan & Chow (2022), [arXiv:2111.11003](https://arxiv.org/pdf/2111.11003), Appendix A and Tables 1–2. The first three cases correspond to B0, B5, B7 in Table 1. The paper frequently obtains both signs together; anticyclone formation is not proof of exclusive anticyclonic selection.

## Declared configuration

- Dry fully compressible Breeze, pinned commit `0479d9c8ce006adcb1e3314c324b84980ad4871b`; ideal gas R=1, cp=2.5, gamma=5/3. These are benchmark units, **not the earlier hydrogen–helium constants**.
- Domain 4 × 4 × 1; periodic x/y; impermeable free-slip top and bottom. No beta, channel walls, stable cap, temperature relaxation, imposed/maintained shear, moisture, or momentum drag.
- Initial polytrope T=5−4z, rho=T, p=T², g=8; top reference values T=rho=p=1. Initialize from the discretely hydrostatic Breeze reference and record any discrepancy from the analytic profile.
- Identical seed 31415, zero initial velocity, relative potential-temperature noise 1e−5 with horizontal mean removed and vertical sin(pi z) envelope. No pre-inserted LSV. A passive tracer blob is used solely for visualization.
- Bottom total heat flux F=0.001; constant **temperature conductivity** K=F/4=0.00025. Fixed top temperature T=1 is imposed through the upper finite-volume conductive flux. No additional theta-diffusion source is used.
- Rotation rates Omega=0, 0.64, 2.56, all at 90 degrees for this first benchmark. ConstantCartesianCoriolis retains all components when later used at a tilted latitude. No low-latitude case is included in the present batch.
- Local Smagorinsky viscosity nu=0.28² sqrt(dx dy) dz sqrt(2 S:S). The native coefficient is adjusted for Oceananigans' volume-based filter convention. SGS scalar Pr=Inf because heat conduction is provided separately. Viscous heating 2 rho nu S:S is included as an energy source.
- Grid 128 × 128 × 64 for this implementation batch versus the paper's 512 × 512 × 101. A negative low-resolution outcome does not refute the mechanism.
- WENO5 finite-volume advection, Float64, Breeze split-explicit acoustic stepping (including its default thermal divergence damping), CFL 0.25, outer dt ≤ 0.02. These replace the paper's horizontal spectral / vertical finite-difference semi-implicit algorithm. The full local SGS stress also differs from their simplification of fluctuating-viscosity terms. Potential-temperature evolution plus explicit heating is not a discretely conservative total-energy integrator; monitor total-energy drift rather than claim exact closure.
- Target t=5000, or five nominal 1/F thermal times. This is **not a claim of achieved thermal equilibrium**. Stop after at most four allocated GPU-hours per case, retain restart checkpoints, and report actual final time. No automatic retries or extensions.
- Three-depth vorticity/tracer movies at cell centers nearest z/H=0.2,0.5,0.8. Save frame every 10 time units and checkpoints every 100. Raw data and rendering remain on Unity.

## Checks and decisions

Before production: manufactured temperature-conduction solution, discrete integrated heat source equals boundary heat-flux difference, hydrostatic initialization, finite positive thermodynamics, actual GPU integration, and continuous-versus-serialized-restart agreement. Check the target-grid cost. These are implementation checks, not evidence of physical fidelity.

Interpretation requires all of:

1. Sufficient thermal equilibration and turbulent evolution; report top/bottom heat flux and energy trend.
2. A measured nonrotating reference velocity for Cai's Ro_c=U_nonrot/(2 Omega H). Do not substitute seeded-eddy Ro or a local rotating-flow Ro.
3. A resolved convective spectrum and a coherent vortex larger than the convective scale. Skewness alone is insufficient; examine persistent circulation, area, velocity, and its vertical structure. Compare both vortex signs relative to f.
4. A later resolution check before accepting a sign transition. Absolute low-mode energy is not proof of inverse transfer; any energy-supply claim requires a separately verified transfer diagnostic.

## Conditional next experiments, not submitted in this batch

Only after the benchmark produces a coherent sufficiently resolved anticyclone: test a tilted rotation vector at the GRS latitude, then one-at-a-time stratification and jet changes. For stratification, distinguish maintenance of the horizontal mean temperature profile from local damping of temperature anomalies; preserve equal mean heat removal when comparing them. For jets, measure the actual surviving shear, include a no-jet control, and distinguish freely evolving from explicitly maintained jets with their momentum/energy input recorded. Formation from small perturbations and survival of an inserted vortex answer different questions and must remain labeled separately.

No such follow-on may be released merely because skewness becomes negative. All changes and resource limits must be stated before submission.
