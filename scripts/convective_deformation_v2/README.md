# Convection/stable-layer pilot v2: mean-profile thermal repair

Authorized November 26. This is a new driver; v1 files and completed outputs remain unchanged.

The two dry anelastic cases retain the original 256×256×80 grid, 8×8×1 domain, vertical f=0.1, beta=0, no jet, identical depth-spanning mixed-sign eddies at Ro=0.1, WENO9, constant ordinary viscosity/diffusivity (eddy Re=100, Pr=1), bottom enthalpy flux and tracer. Initial deformation ratios L/Ld=0.5 and 2 are calibrated to the same full-column displacement eigenproblem. Periodic horizontal directions and impermeable stress-free, zero-scalar-flux lids remain.

## Only physical change

Remove the prescribed top-10% cooling slab. Add directly to rho*theta:

    S_rhotheta(z,t) = -rho0(z) m(z) [mean_xy(theta)(z,t)-theta_ref(z)] / 20.

The smooth m(z) rises from zero at z/H=.55 to one at .65, matching the background cap transition. The initial cap reference theta is fixed. This source is uniform horizontally and therefore has no direct theta-anomaly damping and no momentum or vorticity-sign forcing. It may heat or cool; it represents an externally maintained mean stratification, not a self-consistent radiative equilibrium. Source energy is recorded and must accompany interpretation. The lower convection zone is not restored. Bottom heating remains unchanged.

Mean theta is refreshed on every model update stage and after restart. The cap source is converted to enthalpy-equivalent power using the actual cp*Pi. Instantaneous signed, positive and negative powers are saved; cumulative net cap energy uses trapezoidal quadrature at completed timesteps and is checkpointed. This is source accounting, not an exact total-energy budget for the potential-temperature formulation.

## Release sequence

1. Full-file parsing; both full-grid initial conditions and balance checks; actual driver outputs/readback/restart on CPU; mean-only source sign and spatial-uniformity check.
2. Horizontally uniform 8×8×80 columns at production vertical resolution and the production-calibrated physical parameters through t=1000. No eddies or thermal noise, same bottom heating and cap controller. Require completion, negligible vertical velocity, and interior-cap minimum mean N² >80% of initial. Passing source/environment receipt is created only after these pass.
3. Repeat GPU restart/readback checks, then the two 256×256×80 pilots to t=1000 (10 initial turnover times), at most one GPU-hour allocation each. No automatic long extension. Three-depth planes every .1 turnover; volumes/checkpoints every2.5. Exact endpoint saved.

Afterward assess evolving N² and internal deformation speed, lower convective heat flux, vertical-velocity spectra and sign statistics. A resolved convective regime and retained cap are prerequisites for a long comparison. A targeted 512² check, if needed, is separate and must preserve physical parameters. Merely passing numerical stability does not establish LSV physics or sufficient resolution.

Column validation uses absolute mass-divergence <1e-12 and |w|<1e-10: a relative residual is singular for an exactly resting column. Both absolute and relative values are logged. The 3D relative constraint tolerance remains 1e-8. Restart validation also compares cumulative cap-source heat.
