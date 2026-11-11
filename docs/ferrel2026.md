# Midlatitude circulation mechanism pilot

**Submitted on Unity:** initialization job 1401859 passed; quick closure/restart
gate 1401860 and full-grid gate 1401861 precede the 48-hour experiments
1401863–1401865. Fine check 1401862 branches after the full-grid gate.
See the [launch record](new-campaigns-launch.json) for the latest recorded status.
Submitting this chain is not a claim that its circulation has physical fidelity.

This isolated dry anelastic Breeze campaign asks whether convection changes the eddy torque and overturning associated with weak midlatitude jets. It is a **48-hour transient mechanism and diagnostic test**, not a steady Ferrel-cell simulation, a reproduction of Jupiter's deep circulation, or evidence that its stacked cells have been explained. The preceding 85°N convection pilots establish numerical readiness; they do not establish midlatitude fidelity.

## Configuration and matched experiments

The channel is 12,000 × 12,000 × 120 km, centered at 55°N, with periodic zonal boundaries and free-slip, impermeable meridional and vertical boundaries. Jupiter's rotation rate and the configured radius give approximately `f0 = 2.881e-4 s^-1` and `beta = 2.822e-12 m^-1 s^-1`. The traditional beta-plane approximation is used; deep spherical geometry and the horizontal component of planetary rotation are outside this pilot's scope.

The initial barotropic cosine jets have amplitude **2 m/s** and wavelength **6,000 km**, giving two periods across the channel. Their maximum curvature is about 0.78 times beta, so the initial zonal-mean barotropic potential-vorticity gradient does not reverse sign. This avoids deliberately releasing a strongly unstable, energetic jet as the dominant energy source. The jets are initialized only: no jet nudging, maintained mechanical forcing, or additional drag acts subsequently.

Matched full-depth random eddies have central wavelength 1,200 km and Rossby number 0.1 defined using the inverse angular wavenumber. Their initial RMS speed is approximately 5.5 m/s. Initial jet kinetic energy is therefore smaller than the random-eddy kinetic energy. The jet's own shear Rossby number is approximately 0.007; the eddy and jet Rossby numbers are distinct.

| Experiment | Heating | Initial eddies | Duration | Purpose |
|---|---:|---|---:|---|
| Axisymmetric unheated control | 0 W/m² | None | 6 hours | Jet balance, numerical adjustment and diagnostic baseline |
| Axisymmetric heated control | 162 W/m² | None | 6 hours | Direct thermal and boundary response without imposed nonaxisymmetric seeds |
| Unheated turbulent control | 0 W/m² | Matched random realization | 48 hours | Eddy–jet evolution without bottom heat input |
| Heated turbulent case | 162 W/m² | Same realization | 48 hours | Convection-related change relative to the unheated case |
| Combined boundary sensitivity | 162 W/m² | Same realization | 48 hours | Double the cap relaxation time and disable the vertical-velocity sponge |
| Fine-grid early check | Selected heated case | Same physical initialization | 6 hours | Construction, cost and early resolution sensitivity |

Screening uses 512 × 512 × 96 cells, with horizontal spacing about 23.4 km. The fine check uses 1024 × 1024 × 96. Quick validation uses 96 × 96 × 48 and a 3,200-km seed wavelength to retain the initializer's minimum cells-per-wavelength requirement; it is not the scientific configuration. Six hours at fine resolution cannot establish convergence of a 48-hour circulation or a mature convective state.

Axisymmetric controls remove both the random velocity field and all three-dimensional potential-temperature perturbations. They test the imposed symmetric initial state, not the physical stability of a perfectly symmetric atmosphere. The jet pressure is constructed by discrete integration of the native meridional Coriolis balance. With height-independent zonal velocity and traditional rotation, no imposed meridional thermal-wind temperature gradient is required. The random eddies retain a leading-order f0-geostrophic initialization, so beta-related and nonlinear adjustment remain and must be separated from later responses.

## Thermal treatment and limitations

The lower layer is initially neutral. A stable cap begins above 80 km, with thermal relaxation transitioning to its full strength above 90 km. The baseline includes a vertical-velocity sponge above 108 km and the existing globally uniform fractional potential-temperature feedback that balances net imposed enthalpy input. Heating at 162 W/m² is intentionally enhanced; it is not Jupiter's calibrated intrinsic heat flux. The model is dry and contains neither cloud microphysics nor realistic radiative transfer.

The cap, sponge, thermal feedback, numerical diffusion and rigid boundaries can influence the circulation. The combined boundary sensitivity is a stress test: changing two settings together cannot identify which one caused a difference. Separate one-factor tests, alternative heat removal, viscosity sensitivity, additional seeds and domain/depth sensitivity are needed before a robust physical attribution. A 120-km channel cannot represent the closure of deep Jovian cells thousands of kilometers below the clouds.

## Diagnostics and interpretation

The driver samples native density-weighted zonal momentum tendencies across complete SSPRK3 steps, initially for the first eight steps and then every 250 steps. Closure must use the actual stage weights and state changes, including the effects of the pressure correction. This is a sampled numerical closure test, not an assertion that every step or every equation has been audited.

Native meridional and vertical face mass fluxes provide the mean mass streamfunction and its continuity check. Additional zonal means and centered covariance fluxes describe resolved horizontal and vertical momentum, heat and passive-tracer transport. The centered eddy-flux decomposition and mean–eddy energy-transfer estimate are approximate; their remainder against the native advection budget must remain visible. A closed total momentum budget does not by itself prove a closed eddy decomposition.

For traditional rotation the continuum reference balance is

\[
\partial_t\bar u+\bar v\,\partial_y\bar u+\bar w\,\partial_z\bar u-f\bar v
=-\partial_y\overline{u'v'}-\rho_0^{-1}\partial_z(\rho_0\overline{u'w'})
+\text{viscous and imposed tendencies}.
\]

Both horizontal and vertical eddy momentum fluxes are retained. The resolved conversion estimate

\[
C_{\mathrm{eddy}\rightarrow\mathrm{mean}}
=\int\rho_0\left(\overline{u'v'}\,\partial_y\bar u
+\overline{u'w'}\,\partial_z\bar u\right)\,dy\,dz
\]

is positive for transfer into the zonal-mean flow. Negative conversion indicates extraction of mean-flow energy; it must not be described as eddies maintaining the jets. Convection modulation is assessed from heated-minus-unheated differences, not from an isolated circulation image.

Eulerian overturning and tracer transport are reported separately. In the nearly neutral lower layer, standard transformed-Eulerian-mean formulas involving division by the mean potential-temperature gradient are ill-conditioned. Use directly diagnosed mean and eddy tracer fluxes there; do not manufacture a residual circulation with an arbitrary denominator floor. TEM may be evaluated only where its stratification and scale assumptions are satisfied.

## Release and continuation criteria

Production follows successful initialization, restart, finite-state, thermal/cap, tracer, GPU-resource, throughput and native momentum/continuity checks on Unity. Passing these checks demonstrates executable numerical behavior, not Jovian fidelity. Diagnostic maxima across the axisymmetric controls matter, not only their final values, because a near-inertial oscillation can return close to its initial state.

The first 48 hours establish whether there is a measurable signal worth extending. A candidate response must exceed the axisymmetric adjustment envelope, survive a comparison with the unheated turbulent case, and be explainable by the measured momentum terms. Early random-eddy adjustment must be identified explicitly. Boundary/viscous terms or decomposition error must not dominate a response attributed to resolved eddies.

Extensions over tens of days require the budget, transfer sign and early resolution checks to be satisfactory. Statistical evidence then requires multiple time blocks and eventually additional seeds; the first pair cannot establish statistical robustness. A useful provisional continuation target is a circulation signal larger than three block-estimated standard errors and a budget/decomposition residual below roughly 20% of the leading terms. These are project screening criteria, not universal physical thresholds or substitutes for refinement.

All simulations, diagnostics and rendering run on Unity. Scientific result entries retain the synchronized three-depth vorticity and tracer movies required by the repository. Status-only reports need not include new media. Existing campaigns and frozen source files are preserved.

## Scientific context

[Duer et al. (2021)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL095651) connect Jovian circulation inferences with eddy momentum transport. [Duer et al. (2023)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023AV000908) show deep anelastic examples with important vertical momentum transport and strong boundary sensitivity. Those simulations motivate mechanism tests; they do not validate this shallow channel.

[Mendonça et al.](https://arxiv.org/html/2509.02759v1) distinguish Eulerian from residual transport and demonstrate why their vertical-motion signs can differ. [Lee et al. (2026)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025GL120433) provide an Earth-based proof of concept for stacked residual cells involving eddy heat transport and baroclinic dynamics. These are competing mechanisms to test, not evidence that convection alone must produce stacked Jovian cells. [Boissinot et al. (2024)](https://doi.org/10.1051/0004-6361/202245220) also show why successful dry/moist jet simulations do not automatically identify the associated Ferrel circulation.
