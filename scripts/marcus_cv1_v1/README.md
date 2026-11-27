# Zhang–Marcus CV1 reproduction: implementation scaffold

Authorized 2026-11-27. This is **not a running CV1 reproduction**. The only executable modes are a CPU preflight and a small GPU implementation check using a manufactured atmosphere. The full-grid physical run is deliberately unavailable until the reference inputs below are settled. Existing LSV/Cai campaigns are unchanged.

## Target and missing inputs

Target: Zhang & Marcus (2024), JFM 984 A61, case CV1, DOI https://doi.org/10.1017/jfm.2024.132. The target configuration is recorded in target.toml. It is the beta=0 constant-shear, tapered, spin-aligned anticyclone, followed to 120 Earth days. No driven convection or additional jets.

The publisher page supplies the paper and three supplementary movies. We did not find numerical initialization tables or source. Before calling this a reproduction we need:

1. Smoothed background z, pressure, temperature, density, potential temperature and N², with g, R, cp and pressure normalization. A native initial-condition file and its format would be even better.
2. Exact vertical domain and boundary placement; the paper reports sponge thicknesses but not a complete unambiguous numerical specification.
3. Sponge amplitude/profile, hyperviscosity and hyperdiffusivity coefficients/orders, and the zonal-flow continuation near the periodic y seam.
4. Confirmation of the physical sign/magnitude of the CV1 shear and vortex. In east–north–up coordinates, the paper's Figure 3 shows dU/dy<0 and its Southern Hemisphere vortex is anticyclonic. But the printed U=-sigma*y with sigma=-9.1e-6 s^-1 gives the opposite shear; equation 5.6 would then give negative omega0 at epsilon=1.5. Figure 5 instead quotes positive omega0. We will not silently copy the inconsistent signs into a production run.

An approximate profile digitized from Figure 2 would be a documented reconstruction, not the original reference state. No such approximation has been used to launch a physical experiment.

## Equation audit against the pinned Breeze source

Pin: environments/chan-total-energy-8273be8; package source location read from SOURCE_DIR.txt. Shared source is not edited.

Breeze's continuity constraint is div(rho0*u)=0 and its pressure acceleration is -grad(p'/rho0), matching the paper's structure. In dry fixed-reference-pressure thermodynamics, set delta=theta/theta0-1=T/T0-1. Native Breeze buoyancy is g*delta/(1+delta), whereas the paper uses g*delta. An isolated vertical momentum forcing adds rho0*g*delta²/(1+delta). It is evaluated at centres and interpolated with the same native vertical interpolation as buoyancy. This correction is tested against the actual native function at positive and negative thermal anomalies. It is not heating and does not alter any existing campaign.

Potential temperature is advected as a density-weighted prognostic quantity; with fixed rho0 and the anelastic constraint, its continuum equation is material conservation. The initial nonlinear horizontal pressure balance includes advection, Coriolis, and full-rotation vertical balance. The reference vorticity is computed from the paper's shielded velocity formula, then inverted for a discrete divergence-free horizontal velocity. A uniform vorticity correction from removing the finite-domain mean is recorded, never hidden. Its acceptability must be checked at the eventual production resolution.

## What this scaffold tests

A 24³ dry manufactured atmosphere at 23S with a known isothermal reference, not the Jovian profile. Domain 2000×2000×80 km, 2 s duration, fixed dt0.5 s. It exercises rest, rotation-aligned periodic shear, tapered vortex initialization, local boundary damping, WENO9 advection, three-depth outputs, diagnostics, and serialized restart. These sizes and times have no scientific interpretation as CV1 or GRS predictions. The trial boundary window/damping values are test fixtures, not claimed to be Marcus's settings.

The source is parsed fully before execution. A receipt fingerprints source/configuration/batch/environment/Julia and the pinned Breeze source. GPU submission requires a successful unchanged CPU receipt. Production is disabled independently of whether the fixture checks pass.

The diagnostic energy uses cell-centred kinetic energy plus rho0*cp*Pi*theta and records the signed sponge power. It is a diagnostic approximation; exact discrete energy closure is not claimed. Full production needs staggered-energy and time-integrated source accounting.

## Remaining implementation before production

- Import and verify the supplied physical reference profile; adapt initialization to Breeze's discrete hydrostatic reference without changing its stratification silently.
- Set the confirmed boundary, circulation, and dissipation conventions; quantify departures from the spectral solver.
- Test a stratified gravity-wave mode at the proposed timestep and density contrast.
- Resolve vortex collar and vertical taper; perform actual-grid memory/throughput validation on Unity before quoting elapsed time for 120 days.
- Match publication figures at 0, 40, 80, 120 days, including vorticity, temperature, vertical shape and circulation. Then run one resolution sensitivity check before adding beta/convection.

All numerical checks and analysis execute on Unity. Local activity is source editing and reading the publication. No CV1 long run or automatic extension is submitted by these scripts.
