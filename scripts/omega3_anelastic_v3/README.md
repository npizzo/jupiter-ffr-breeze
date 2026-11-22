**Superseded by `omega3_anelastic_v4`.** CPU preflight passed, but GPU validation caught a diagnostic reduction using a host-only `mean` on a GPU view. No full-resolution integration was performed with v3. V4 corrects that readout; physics is unchanged.

# Matched Ω=3 anelastic comparison

Purpose: test whether the finite-amplitude anticyclonic preference in the original compressible Ω=3 periodic box also occurs with Breeze's anelastic dynamics. This is not a beta/channel experiment, and not a new cap experiment.

The original `omega3_long_v1` model and initializer are imported unchanged to generate the initial face velocities, potential temperature (specific entropy) and tracer. The anelastic reference pressure, temperature, density and potential temperature are checked against that original reference. The original 256×256×80 grid, 4×4×1 periodic box, tilted rotation vector, gas constants, bottom flux, radiative conduction, cap cooling, Smagorinsky viscosity, deviatoric stress and viscous heating, entropy diffusion, and WENO5 advection are retained.

## What necessarily changes

The continuity density becomes fixed rho0(z); the pressure perturbation is diagnosed, not independently evolved. The thermodynamic pressure is the hydrostatic reference pressure. Matching entropy and velocity therefore cannot also match the compressible temperature, density or dynamical pressure perturbations. The initial mapping records these differences rather than claiming identical thermodynamic states.

Density factors in the anelastic stresses, thermal flux conductances and cooling use rho0, consistent with its mass equation. The entropy gradient must use `cp log(T) - log(pref)` (equivalently `cp log(theta)`), not `cv log(T) - log(rho0)`. The latter would manufacture a spurious entropy gradient. A CPU test compares the constitutive flux to the compressible formula on a common EOS state and checks that it telescopes.

The native volumetric heat-source conversion is checked against `Q theta/(cp T)` in the prognostic rho0*theta equation, with zero velocities and SGS viscosity. This is separate from a claim of total-energy closure, which is not made. The reference-at-rest test, weighted-divergence check, wall impermeability, serialized restart and output readback all run in the actual driver.

## Release sequence

1. Fingerprinted CPU preflight using the actual driver, imports, configuration, outputs and restart.
2. GPU restart/readback validation within the allocation.
3. A full-resolution pilot to t=5 with a 30-minute model wall limit (45-minute total allocation). Save three-depth vorticity/tracer/temperature, full fields and the initialization mapping.
4. Review the actual approximation departures, transient adjustment, pressure constraint, stability and cost before a long continuation. The scientific target is t=1000, matching the existing compressible experiment; this driver does not release it automatically.

Anelastic dynamics retains density stratification and divergence associated with vertical motion. An anticyclone here would show that acoustic waves and prognostic background-density evolution are unnecessary in this configuration. A different sign would not, by itself, prove an acoustic mechanism: initial thermodynamic mapping, anelastic pressure/buoyancy dynamics, and numerical differences would need to be separated. Low Mach number alone is not an adequate validity test in rapidly rotating convection (Calkins, Julien & Marti 2015, https://doi.org/10.1098/rspa.2014.0689).

All simulations, diagnostics and rendering run on Unity. No shared package or active campaign is modified.

The anelastic reference fields are copied from the actual discrete Exner reference used by the compressible model, including halos. They are not independently reintegrated from the analytic profile. Reference potential temperature is derived from these pressure/density fields; the scalar `ReferenceState.potential_temperature` metadata is only the bottom value.
