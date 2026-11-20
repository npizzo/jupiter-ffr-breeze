# Chan–Mayr 2013 reconstruction: Cases 0, 3 and 4

Separate from the frozen `kapyla_pair_v1` campaign. Both campaigns test compressible-convection routes to vortex-sign selection. They do **not** isolate stratification, cooling, and closure effects from one another: multiple model ingredients differ between the papers. A sign reversal would motivate a subsequent matched intervention, not prove a unique mechanism.

Source: Chan & Mayr (2013), *Earth and Planetary Science Letters* 371–372, 212–219, [doi:10.1016/j.epsl.2013.03.046](https://doi.org/10.1016/j.epsl.2013.03.046). The user supplied the full PDF. Chan (2001), [doi:10.1086/318989](https://doi.org/10.1086/318989), describes radiative temperature diffusion and SGS entropy diffusion.

## Physical specification

- Dry three-dimensional compressible ideal gas, gamma=5/3, R=1, cp=2.5, cv=1.5, g=12. Time unit is box height divided by the top isothermal sound speed, **not an Earth day**.
- Box 4×4×1; periodic horizontal directions; impermeable, stress-free horizontal walls. x east, y north, z up. Positive vertical vorticity is cyclonic for the rotating northern-hemisphere cases.
- Interface z=0.78; lower polytropic index n=1.5, upper n=24. Ttop=rhotop=ptop=1 in the analytic reference atmosphere. Density, pressure and temperature continuous at the interface. Model initializes with its discretely hydrostatic reference density; audits record actual grid values.
- Uniform upward bottom heat flux F=0.03125. Radiative conductivity K=0.0004 F/(g/cp). Fixed top T=1 is applied through the half-cell boundary conductive flux. Conductivity is held constant throughout the box (explicit reconstruction choice).
- Upper cooling Q=-rho*(T-1)*S(z)/20, where S is cubic smoothstep from 0 at z=0.78 to 1 at z=1. There is **no cv multiplier** in this volumetric source, following the printed energy equation. The smoothstep profile is an explicit choice because the paper does not give its functional form.
- Native Smagorinsky nu=(0.16 Δ)^2 sqrt(2 S:S), Δ=(dx dy dz)^(1/3). Coefficient and cubic-volume filter are explicit reconstruction choices, not claimed recovered from the 2013 paper. Full compressible deviatoric stress is native symmetric stress plus the -2/3 rho nu div(u) diagonal correction; viscous dissipation is returned as heat.
- Turbulent Prandtl number 1/3. SGS heat flux q=-rho T (nu/Pr) grad(s), s=cv log(T)-R log(rho). Conservative face-flux divergence; zero SGS normal heat flux at both walls. Native thermodynamic/tracer diffusion is disabled to avoid double counting. Passive dye is advected only, with bounded WENO.
- No drag, imposed jets, beta, moisture, or artificial vortex-sign forcing. Same tiny potential-temperature noise (seed 314159, amplitude 1e-5), no initial velocity.
- Case 0: Omega=0. Case 3: Omega=1.5. Case 4: Omega=3. Rotating cases at 67.5N with **full** vector (0,Omega cos(latitude),Omega sin(latitude)).

The printed n, interface, gravity and top normalization yield about 6.2 pressure scale heights across the **whole** domain, not the lower layer alone; the text's lower-layer claim is inconsistent with those explicit parameters. `PROFILE.toml` records the calculation rather than altering gravity to force that sentence.

## Numerical differences and interpretation

Breeze 0.15.0, Oceananigans 0.117.2; Float64, WENO5, native split-explicit acoustics with no additional divergence damping. Breeze evolves potential-temperature density, **not the paper's strongly conservative total-energy variable**. The continuum heating/stress construction does not ensure exact discrete total-energy conservation. Total kinetic+internal+gravitational energy and the boundary/cooling power are recorded; significant imbalance stops a run for review. No closed discrete vorticity budget is claimed here.

Initial production grid is 256×256×80 uniform, below the paper's 540×540×80. SGS viscosity depends on resolution. Consequently numerical differences can change the achievable Reynolds number and vortex transition. Failure to produce an anticyclone at this grid is not a refutation of Chan–Mayr.

Diagnostics report V=sqrt(<u²+v²+w²>cz), Ro=V/(2 Omega d), Re=V d/<nu>cz, with d=0.78. Case 0's Ro field is a finite placeholder (division floor), not a meaningful rotation diagnostic. Once thermally equilibrated, measure V0 from Case 0 and evaluate RobOmega=V0/(2 Omega d), Robf=V0/(f d). The paper's .063/.031 Robf values are **comparison targets**, not enforced by imposing its reported V0. Examine coherent vortex persistence, depth dependence, core heat transport and temperature; neither skewness nor a single movie proves the mechanism.

The paper requires thermal relaxation exceeding 1e4 time units. Target t=20000, subject to demonstrated thermal stationarity. No early segment is to be called equilibrated or a mature sign-selection result.

## Gates and resource caps

1. CPU operator checks: hydrostatic profile identities, full Coriolis vector, conservative radiative/SGS heat flux, zero entropy diffusion for constant entropy, compressive stress correction, short rest/heated integrations.
2. One A100 validation allocation (maximum 30 minutes): same operators, hydrostatic rest, serialized restart equivalence, heated 32×32×80 transients at Omega=0 and 3, then timing on the actual production grid.
3. Only an `afterok` validation result with matching source/config hashes releases the three production jobs. Predicted progress, with a 1.5× timing penalty, must exceed t=1000 per eight-hour allocation. This early estimate does not guarantee turbulent throughput.
4. Each first production allocation: one A100, 4 CPU cores, 48 GB host memory, **maximum 8 hours**, graceful stop at 27900 seconds. No automatic repeated allocation. Review progress, heat balance and cost before continuation.
5. Stop for nonfinite/nonpositive fields, mass error >1e-7, max Mach>3, or interval total-energy imbalance >10% of bottom power after t=1. This energy tolerance is a feasibility screen, not precision budget closure. Endpoint quadrature can itself trigger a review; do not silently relax it.

The SGS diffusion limit is recomputed every step. Checkpoint every 250 units and on normal termination. Three vorticity/tracer planes nearest z=0.2,0.5,0.9 every 25 units; sparse full face-velocity, density and temperature fields every 1000. Diagnostic interval one unit. No raw fields go to GitHub; all scientific computations, validation and future rendering stay on Unity. Movies must use synchronized three-column vorticity/tracer panels and explicit fixed color scales.

These are bounded first segments, not a commitment to spend whatever it takes to complete t=20000. Restart validation tests the supported continuation path; new continuation jobs require output review.
