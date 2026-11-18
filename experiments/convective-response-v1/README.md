# Convective stress response around matched stable vortices

Authorized November 18, 2026. This is a bounded mechanism experiment, executed on Unity. The [equation derivation](DERIVATION.md) and [earlier adversarial review](../feedback-maintenance-v1/literature-review.md) explain the motivation. It does not replace or modify any Breeze campaign.

## Questions and declared approximation

Does convection produce a different stress response around otherwise matched cyclonic and anticyclonic vortices? Does the difference survive removing unstable thermal gradients, reflecting f, and changing the number of retained vertical modes?

The model linearizes finite-f, nonhydrostatic Boussinesq perturbation equations about a prescribed horizontal vortex and a frozen thermal profile. The vortex is held fixed, so stress work is a **virtual reinforcement tendency**, not measured growth of a freely evolving vortex. This experiment cannot demonstrate spontaneous LSV formation, saturated convection, or a sustained inverse cascade.

The reference is the published 64-level day-100 mean N² profile from the heated f-plane run. The physical height is 120 km, f=2.880957739764384e-4 s^-1, and the stable cap occupies approximately 80–120 km. The profile is linearly interpolated between cell centers and held at the nearest center value outside their range. It is not a reconstruction of the exact wall gradient.

Ordinary viscosity is 1000 m²/s. Thermal diffusivity is 1000 m²/s below 60 km, smoothly tapers to 1 m²/s by 80 km, and remains 1 in the cap. Thermal cap relaxation ramps from zero at 80 km to a one-day relaxation time at 90 km. Momentum boundaries are impermeable and stress-free; thermal perturbations have homogeneous Neumann conditions corresponding to unchanged spatially uniform flux boundaries. No beta, artificial momentum drag, polynomial reaction, or stochastic forcing is added.

This Boussinesq approximation omits the original anelastic density weighting, which varies by nearly a factor of ten across the source domain. It also omits the original vertical-velocity sponge and global enthalpy controller. Results therefore test a proposed mechanism; they are not a quantitative reduction of the full Breeze simulation.

The comparison background replaces negative N² values by zero and retains all positive stratification. It is named **gradient removed**. Turning off the bottom heat flux while leaving an identical frozen gradient would leave the linear equations unchanged; the comparison is not described as an instantaneous heat-flux switch-off.

## Vertical reference and mode selection

A one-dimensional reference solves the pressure-projected rotating perturbation equations over 33 predetermined horizontal wavenumbers, 0.5≤kH≤64. The initial 32-velocity/65-thermal-function calculation is compared with 64-velocity/97-thermal-function resolution. A separate analytic uniform-background dispersion benchmark and pressure/energy identities verify the operator.

The reference screen completed on Unity in 1.54 seconds. Its maximum growth rate changed from 0.20455f to 0.20622f on refinement; the selected kH is 29.9867. See [reference-screen.json](reference-screen.json). The fastest wavelength is approximately 25 km, comparable to the source Breeze horizontal cell width of 23.4 km. This further limits direct quantitative comparison with that simulation; the new horizontal response calculation resolves the selected wavelength separately.

One or two complete convective eigenmodes are then retained. Their vertical-velocity shapes are supplemented by the derivative and toroidal horizontal shapes required to retain finite-f coupling. Thus one physical mode can require more than one horizontal shape function. Independent cosine buoyancy modes satisfy the original thermal perturbation boundary condition. Actual function counts and selected thermal rank are recorded in the mode validation output.

Release requires convergence of the full vertical reference; growth of the intended convective branch; less than 20% growth-rate error across 0.8,1.0,1.2 times the selected wavenumber; agreement of one versus two retained modes; thermal-rank refinement; a non-growing gradient-removed control; and pressure, reflection, and energy checks. A failed check does not permit tuning the measured gradient until it passes.

## Matched response design

- Horizontal periodic grid: 128², square width 16 selected convective wavelengths.
- Prescribed compensated Gaussian vortex: radius 2 wavelengths, peak |ζ/f|=0.25. Opposite signs have identical amplitude and geometry. At f>0, positive ζ is cyclonic.
- Both backgrounds must pass the full axisymmetric inertial-discriminant check with a positive margin, plus a check on the small periodic-image correction.
- Perturbations: matched, small random initial disturbances near the selected convective scale; no subsequent stirring. Two fixed seeds are shared across signs and thermal backgrounds.
- One and two complete vertical modes, two vortex signs, two thermal backgrounds and two seeds give 16 primary cases. Two reflected-f cases verify symmetry over the same response window.
- Planned window: 0≤|f|t≤10, approximately 9.64 hours. Ordinary RK4 timestep at most 0.1/|f|, reduced if required by measured stability. No automatic extension.
- A production-driver pilot must pass output, temporal refinement, energy accounting, symmetry and cost checks. Projected total cost must be below 30 CPU-minutes; the complete production campaign has a one-CPU, 45-minute hard limit.

## Diagnostics and decision

Report perturbation amplification, actual stress work on the prescribed vortex and that work normalized by fluctuation energy. Independently evaluate stress-divergence and stress-strain forms. The perturbation kinetic budget must separate virtual barotropic exchange, buoyancy conversion and viscous loss. Report temporal-integration residuals explicitly. Buoyancy variance is tracked separately; with a sign-changing N² profile it is not positive available potential energy.

Save the actual reconstructed perturbation vorticity at 20,60 and100 km, and buoyancy perturbations at the same heights. These are modal reconstructions, not independent vertical grid levels or passive-tracer concentrations. The small perturbations do not substantially alter the prescribed vortex during this test.

A growing disturbance alone is not evidence that it feeds the vortex. A response that changes qualitatively with vertical rank or seed is insufficient to motivate a long simulation. If both vortices lose energy to the disturbances, this candidate does not explain their maintenance on its own. A robust difference while both vortices remain centrifugally stable would support investigating a feedback beyond selective centrifugal destruction.

## Performance gate history

The first complete response validation passed the physical and numerical checks but estimated 77 CPU-minutes for 18 cases, so production was withheld. The revised driver retains the equations and grid, tests a maximum timestep 0.1/|f| with independent refinement, and records scalars every 0.2/|f| and planes every 0.5/|f|. Production must use the exact timestep certified by the repeated accuracy and 30-CPU-minute cost gate.

## Launch record

Vertical reference 1956711, complete-mode check 1956981 and response validation 1957540 passed. The revised timestep changes the checked response by 1.8e-6 relative to its half-step reference. Measured projected cost is 779.6 CPU-seconds for all 18 cases. Production 1957606 was submitted on one Unity CPU, with analysis/rendering dependent on its successful completion. See [execution.json](execution.json), [mode validation](mode-validation.json) and [response validation](response-validation.json). This record identifies submission, not completed results.
