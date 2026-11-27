# Convection beneath a stable layer: matched deformation pilot

Authorized November 26, 2026. Two dry anelastic Breeze experiments test whether convection below a stable layer can supply persistent anticyclonic structures in that layer. These are exploratory pilots, not reproductions of Chan or predictions of the Great Red Spot.

## Fixed parameters

- 256×256×80 grid in an 8×8×1 nondimensional box; horizontally periodic, impermeable stress-free lids and no scalar boundary flux.
- Vertical f=0.1, beta=0, no imposed jet. WENO9; ordinary constant viscosity and diffusivity, initial eddy Re=100, Pr=1. These transport choices are the same between the new cases, but differ from the unforced deformation pair (Re=2000, zero explicit scalar diffusivity).
- Depth-spanning, identical mixed-sign velocity fields, initial eddy Ro=0.1. This barotropic initialization avoids the large initial baroclinic eddy-energy reservoir of the earlier pair. Pressure includes geostrophic and nonlinear horizontal balance. Thermal perturbations confined below the cap are smooth resolved horizontal modes (6,5), amplitude 1e-7 relative theta; no grid-scale random thermal noise.
- Nearly adiabatic lower 60%; smoothly increasing background N² across z/H=0.55–0.65 and constant N² above. No local thermal restoration or imposed vortex sign.
- Identical horizontally uniform enthalpy input in the bottom 10%, and equal integrated removal in the top 10%. Layers are normalized by their actual discrete thickness. No net prescribed energy input across the whole box. Cooling nevertheless changes entropy/potential temperature; theta mass is not conserved and is not used as a conservation guard.
- Bottom buoyancy-flux scale B0=U³/0.6; corresponding enthalpy flux F0=B0*rho_bottom*cp*T_bottom/g. This targets a dimensional-analysis convective velocity comparable to the initial eddy velocity; actual convection is measured, not assumed.
- Passive tracer starts as C=0.1+0.05 cos(2πx/Lx)cos(2πy/Ly). It has the same ordinary diffusivity as the thermal scalar.

## Controlled difference

Calibrate the initial stable-layer N² to L/Ld=0.5 or 2, where L is the discrete velocity-gradient length and Ld=c1/f. The first internal hydrostatic mode uses the full reference column, including the neutral lower region. Its displacement eigenproblem is

    -(a_z/rho)_z = N² a/(rho c1²),  a=0 at both lids.

This formulation remains regular when N² vanishes below the stable layer. Both initial velocities match to roundoff. Reference density and stable-layer heat capacity are allowed to differ as required by hydrostatic balance. Altering stratification also changes plume penetration, mixing and transmission of energy: the pair is a test of the combined convection/stable-layer route, not an isolated proof of Polvani's mechanism.

## Duration, checks and scope

Initial pilot: 10 initial eddy turnover times, approximately native t=1000 (tau=L/U=100). Each GPU job has a one-hour allocation and a shorter runtime cutoff to leave checkpoint time. No automatic long extension.

CPU preflight parses the entire actual driver, constructs both full production initial states, checks velocity matching and initial adjustment, then runs all enabled diagnostic/output/checkpoint paths and serialized restart equivalence on a smaller grid. It checks equal integrated prescribed heating/cooling. A passing receipt binds source, environment, pinned Breeze and Julia. GPU restart checks repeat within each production allocation.

The pilot must demonstrate resolved, active convection, distinct retained stable-layer stratifications, acceptable numerical checks and measured affordable cost before a longer run. Changes to N², convective heat flux, horizontal vertical-velocity spectra and sign statistics must be reported. The viscosity and potential-temperature formulation mean an exact total-energy conservation claim is inappropriate; prescribed heat source/sink powers and thermal drift are recorded separately.

Output: diagnostics and three-depth planes every 0.1 turnover at z/H nearest 0.3, 0.7, 0.85; volumes/checkpoints every 2.5 turnovers. Exact final output is saved. The upper cooling layer is reported separately from the interior cap; the entire column remains subject to Mach, finite-state and positive-temperature guards. Scalar conservation checks apply to the passive tracer, not the heated theta variable.

The implementation is frozen while jobs run. Neither new case changes any existing simulation.
