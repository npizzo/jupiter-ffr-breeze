# Compressible f-plane sign-selection pair

Authorized November 19, 2026. Target: Käpylä, Mantere & Hackman (2011),
ApJ 742, 34, https://doi.org/10.1088/0004-637X/742/1/34, cases A9 and A1.
This is a documented reconstruction in Breeze, **not an exact replication**.
No beta, imposed shear, magnetic field, moisture, inserted vortex, drag,
hyperviscosity or SGS viscosity is present.

## Equations and inputs

Dry compressible ideal gas, gamma=5/3. Units d=g=cp=rho(z=0)=1; R=0.4,
cv=0.6. Constant kinematic viscosity, with stress
`2 rho nu [sym(grad u) - div(u) I/3]`, and matching viscous heating.
The native density-weighted ScalarDiffusivity supplies the symmetric stress;
three explicit momentum forcings supply only the missing isotropic trace
subtraction. A manufactured compressive flow verifies the combined 4/3 factor.

Domain 4d by 4d, z/d=-0.85 to 1.15. Periodic horizontally; impermeable,
stress-free upper/lower walls. Piecewise polytropic reference indices 3 below
z=0 and 1 above, with xi=R T_top/(gd)=1/3. The hydrostatic reference is built
by Breeze from the analytic potential temperature and pressure at z=0.
Identical mean-free random thermal perturbations (relative amplitude 1e-5,
seed 314159) restricted to 0<z<1. Initially zero velocity. A horizontally
varying passive dye is dynamically inert. No pre-relaxed progenitor is available;
both cases include their own thermal spin-up, which must be excluded from claims
of statistically stationary sign selection.

Both cases set reported Ra=2e6 and Pr=0.24 using the analytic unperturbed
mid-convection-zone entropy gradient. Specifically, `chi=sqrt((-ds/dz)/(Ra Pr))`,
`nu=Pr chi`, `K_cz=chi rho_mid cp`. K below the convection zone is twice K_cz,
so the piecewise initial temperature slopes carry equal conductive flux.
The bottom flux is `K_cz/(2R)`. Conductivity across the layer junction uses
the series-resistance formula. The top boundary temperature is fixed through
the finite-volume conduction operator, not a second prognostic-temperature
boundary source. Upper-layer cooling adds `-rho cv (T-T_top)/tau(z)`.

The paper does not specify a numerical cooling time or transition width in
its main description. Here the cooling time is 0.1 sqrt(d/g), with a smoothstep
transition over z/d=0.975..1.025. These are explicit implementation choices.
The quoted density contrast (~23) and normalized input flux (1.7e-5) are **not
independently enforced or presumed recovered** by reconstructing the printed
polytropic/xi/Ra/Pr definitions. This ambiguity is a physical fidelity limit:
record initial and evolving profiles and actual heat flux, and do not label a
different achieved regime an exact reproduction or use a negative result to
refute the published result. Further parameter changes require disclosure.

A9 target: Ta=3.6e7. A1 target: Ta=4e8. In each, f=nu sqrt(Ta), so f_A1/f_A9=10/3.
The only difference in the two input configurations is Ta (and the case label).
Reported Co=3.4/15.5 and Re=45/33 are **outcomes**, not imposed values.
Diagnostics use `u_conv=sqrt(3 mean(w^2))` in 0<z<1, `Co=|f|/(2pi u_conv)`,
`Re=u_conv/(2pi nu)`. The paper's precise averaging convention should be
considered when comparing these outputs.

## Numerical method and resource limits

256 x 256 x 128 Float64, WENO5 finite-volume advection; differs from Pencil's
high-order finite differences. Breeze split-explicit compressible stepping,
with native acoustic off-centering but **NoDivergenceDamping**, no sponge.
The outer step is limited by advection, explicit diffusion, cooling and rotation.
The two first production segments each have an eight-hour A100 allocation,
with simulation stop/checkpoint after 7.75 hours. Long-time target is t=30000
sqrt(d/g); reaching the wall limit is not completing that target. No automatic
continuation, retry, or resolution reduction. The initial target-grid cost check
requires at least t=500 estimated progress per eight-hour allocation, using a
50% timing margin. Early-time timing does not guarantee turbulent runtime.

The current thermodynamic variable is potential temperature, so this is not a
discretely conservative total-energy formulation. Record heat input/removal and
total energy. It does **not** yet provide an exactly closed split-explicit
vertical-vorticity budget. This first pair tests the existence of the sign
transition. Mechanistic claims require a separately validated native term budget;
saved sparse fields or vorticity skewness alone cannot supply that proof.

## Release and output

Before release: parameter identities; conductive equilibrium and integrated
boundary-flux balance; manufactured deviatoric stress; source-free hydrostatic
rest; heated transient; continuous vs serialized-checkpoint restart; actual
production-grid GPU smoke and timing. Nonfinite/negative states, mass drift,
or Mach>0.7 stop the run. The validation receipt includes a model hash and the
launcher rejects changed source. All computations stay on Unity.

Movie fields every 25 time units at z/d closest to 0.2, 0.5, 0.8: vertical
vorticity, passive concentration, temperature. Save physical heights and f.
Future rendering uses synchronized three-column vorticity/tracer panels and
fixed, labeled color maps. Checkpoints every 250 and at wall stop, retaining
the latest. Three-dimensional velocity, density and temperature every 500.
Analyze persistent structure, scale separation, actual convection and rotation
regime, thermal balance, and sign relative to f. Low-mode energy or morphology
alone is not proof of inverse energy transfer. A later reflected -f control is
conditional on a credible sign-transition result; it is not submitted here.
