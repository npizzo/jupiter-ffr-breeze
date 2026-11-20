# Plan: a defensible, bounded test of anticyclonic LSV formation

**Status, 20 November 2026: design review. No replacement simulation has been submitted.**

The proposed tenfold reduction of the initial unstable temperature contrast has been prepared locally but is not approved here as a scientifically justified next experiment. The earlier long compressible campaigns remain held or cancelled. This document proposes a decision process, not an automatic launch sequence.

## Scientific question

Can stratified, rotating convection in Breeze support a persistent, concentrated anticyclonic large-scale vortex, and what changes its preference relative to a cyclone?

The immediate objective is an affordable and interpretable local model of that transition. Exact replication of a published calculation is not required. However, adding stratification, lowering a nominal Rossby number, or inserting a vortex does not guarantee that the mechanism in a published anticyclonic case has been retained.

The experiments remain on a **constant-f plane**. No beta, imposed shear, moist physics, or artificial sign-selective vorticity forcing is proposed in this sequence. The beta-plane Ferrel-cell work is separate. A successful anticyclone would be relevant to one ingredient of Jovian vortex formation, not a reproduction of the Great Red Spot.

## What went wrong in the decision process

The first anelastic pilot was checked for implementation correctness before submission, but its parameters were not sufficiently justified as a scientifically suitable convective state. These are different requirements.

The immediate response—reduce the initial temperature contrast by ten—would change a startup transient without changing the bottom heat supply. That might avoid an initial overshoot, but it might also postpone the same outcome or spend the allocation waiting for a new thermal boundary layer. There was not enough evidence to choose that factor.

A further reporting error was to treat the recorded gradient-based Rossby number as directly comparable with published local Rossby thresholds. Definitions of length, velocity, averaging and rotation must be matched first. In particular, a gradient length and a wavelength differ even for a single sinusoidal mode.

## What the completed pilot actually establishes

Unity job **2029809**, on an A100 80 GB, used **4 minutes 32 seconds** of its one-hour allocation. GPU rest, restart and initial timing checks passed. The run then stopped at its deliberately imposed Mach-number guard at nondimensional time **t=74**, with peak Mach number approximately **0.326**. It had accumulated only approximately **0.36 depth-based turnover times**.

This is neither a successful convective parent nor evidence against anticyclonic LSV formation. It is an early transient that exceeded a chosen approximation-validity screen. It did not establish a settled convective state, an inverse transfer, or a sign preference. A Mach threshold of 0.3 is an operational stop criterion, not a universal boundary of validity of the anelastic equations.

Two configuration facts make a blind initial-condition adjustment particularly weak:

- With bottom reference density 1, the imposed flux and diffusivity require a diffusive bottom gradient of magnitude `F/(rho*kappa) = 0.0375`. The original linear lower-layer profile has `Delta_theta/d = 0.03/0.8 = 0.0375`. Reducing only its contrast to 0.003 creates a tenfold mismatch with that boundary requirement. The original profile is not a globally steady conductive solution either, because density varies with height and the upper layer relaxes thermally.
- The depth diffusion time `d^2/kappa = 25,600` is much longer than the pilot target `t=2000`. Convection can accelerate thermal adjustment, so this is not a prediction of equilibration time; it is a warning that an arbitrary initial thermal profile may retain a long memory.

The last saved diagnostics show convective transport at mid-depth substantially exceeding the imposed bottom input during startup. Local flux need not equal boundary input in a transient. This is consistent with release and redistribution of stored thermal buoyancy, but does not identify the cause of the peak Mach number by itself.

## 1. Audit the existing run before another simulation

Use existing data and source on Unity. Do not advance the model or allocate a new GPU for this audit.

1. **Verify the reported Mach number.** Check the dry thermodynamic constants, the temperature used to calculate sound speed, velocity staggering, and interpolation. Locate high-speed regions in saved three-dimensional checkpoints; distinguish bulk convection, the low-density upper layer and localized plumes. The final state at t=74 was not saved as a full checkpoint, so do not claim its exact peak location from earlier snapshots.
2. **Separate startup from sustained forcing.** Inspect horizontally averaged potential temperature, convective and diffusive transport, cap relaxation and thermal storage. Quantify how far the starting profile is from a state compatible with the imposed transport. A closed integral thermal budget alone does not establish local balance or equilibration.
3. **Measure rotation and scale consistently.** Compute a vertical-velocity spectrum from a saved full state, specify cyclic versus angular wavenumber, and derive a spectral horizontal scale. Report both depth-based and convection-scale Rossby numbers. For literature comparison, use that paper's precise definition; report stratified volume-weighted and mass-weighted variants explicitly rather than mixing them.
4. **Check resolved scales.** Inspect spectral falloff and vertical boundary-layer structure. The existing gradient-derived wavelength is only about nine horizontal cells near the stop. That is a reason to inspect resolution, not proof of either convergence or failure. Stable integration is not a resolution test.

**Deliverable:** a compact audit stating what is measured, what is inferred, and what cannot be recovered. If a diagnostic is wrong, correct it without changing the physical equations. If the saved data are insufficient, name the missing quantity before proposing a short diagnostic run.

## 2. Choose an admissible regime, not a convenient temperature perturbation

Before selecting new numbers, document a consistent parameter set in terms of:

- density stratification and stable-layer strength;
- Ekman and Prandtl numbers, with explicit conventions;
- thermal driving and a justified estimate of distance from convective onset;
- predicted velocity, Mach number and convective horizontal scale, including uncertainty;
- anticipated Reynolds number and convection-scale Rossby number;
- cells across convection and thermal boundary layers, plus horizontal scale separation;
- predicted cost per convective turnover and thermal-adjustment time.

Boussinesq onset formulae may provide an order-of-magnitude screen, but cannot certify onset in this stratified, flux-driven, capped model. Likewise, a literature threshold measured in a different geometry is guidance, not a pass certificate.

The tradeoffs must be stated explicitly. Increasing rotation may suppress convection or push its horizontal scale below the mesh. Reducing thermal forcing may lower Mach and Rossby numbers but also remove the turbulence needed to feed a large-scale mode. Reducing viscosity may restore Reynolds number while increasing resolution requirements. There must be a plausible overlap of these requirements at the affordable grid.

**Decision:** if such an overlap cannot be defended, do not launch another 128²×64 run simply because it is cheap. Choose between a justified grid change, a physically different parameter set, or a return to a compressible formulation. Any changes to bottom forcing, cap, viscosity, rotation or geometry must be listed before submission.

Anelastic validity also requires more than small flow Mach number. Rapid rotation, stratification, superadiabaticity and thermal diffusion matter; the cited sound-proof studies below should constrain the choice rather than being treated as blanket approval.

## 3. One calibration experiment, only after the audit

If the design survives the preceding checks, run one bounded calibration with the exact proposed equations and thermal setup. Choose its initial mean thermal profile using the audit, rather than another arbitrary factor. Do not import beta or sign-selective forcing.

**Proposed cost limit:** at most 10 GPU-minutes including startup, with no automatic retry or continuation. This is a proposed ceiling, not a statement that useful thermal equilibration can always be achieved in that time. If the measured cost makes even calibration infeasible, stop before submission and revise the plan.

Before launching, define the observation window and tolerances for:

- low Mach number and small reference-state departures;
- sustained convective transport rather than initial potential-energy release;
- sufficiently stationary mean temperature and velocity statistics;
- a measured, resolved convective scale and sufficiently strong rotational influence;
- a stable upper layer;
- cost per turnover at the developed velocity, not at near-zero startup velocity.

These are scientific acceptance conditions. They are additional to rest, restart, finite-state and conservation checks. The numerical values of the acceptance tolerances must be justified and recorded before the run; this README intentionally does not pretend they are already established.

A smaller horizontal box could make calibration cheaper, but would only be a local-convection test: suppressing large scales changes feedback on convection. It could not certify the existence or absence of an LSV in the larger box.

## 4. A matched vortex-response experiment

Only once a usable parent exists, copy the same convecting state into three branches:

| Branch | Initial large-scale perturbation | Question |
|---|---|---|
| Control | None | What does convection generate without an imposed vortex? |
| Cyclonic | Weak compensated cyclone | Does convection amplify or erode this sign? |
| Anticyclonic | Corresponding compensated anticyclone | Does convection amplify or erode the opposite sign? |

Match seed radius, depth structure and kinetic energy. Respect periodic circulation constraints. Construct an appropriate pressure/buoyancy balance for each sign; a naive velocity addition may generate adjustment transients that masquerade as sign selection. There is no sustained seed forcing and no artificial anticyclone damping.

Measure the large-scale kinetic-energy tendency and the resolved convective momentum transfer into the seeded mode. Distinguish maintenance by convection from simple persistence of an initially imposed vortex. Compare matched turnover windows and quantify differences relative to background variability; one realization is exploratory evidence, not a statistical demonstration.

A 10–20-turnover window is an initial response test, not a guarantee of LSV formation. A slowly growing instability may require longer; continuation must be justified by a measured growth rate, uncertainty and estimated cost, not by waiting for a desired-looking movie. Branch budgets are to be specified after calibration, not reserved automatically now.

## 5. What would count as progress toward the sign transition?

- **Anticyclone survives but receives no demonstrated convective support:** insufficient evidence; possibly just slow viscous decay.
- **Both signs gain comparable support:** evidence for large-scale amplification, not sign selection.
- **Anticyclone receives systematically stronger support:** a candidate sign-selection mechanism. Test seed-amplitude dependence and independent parent states before concluding robustness.
- **An unseeded run produces a persistent anticyclone in the same regime:** stronger evidence of spontaneous selection. Check that its sign is relative to f, and follow with a properly reflected f-reversal control.
- **Neither sign grows:** report the negative response at those parameters. Do not keep changing unrelated terms to manufacture a vortex.

A seeded response and spontaneous formation are distinct questions. The seed can shorten a response experiment; it cannot eliminate the need to demonstrate spontaneous selection if that is the final claim.

## Literature anchors and their limits

- [Guervilly, Hughes & Jones (2014), Large-scale vortices in rapidly rotating Rayleigh–Bénard convection](https://arxiv.org/abs/1403.7442): motivates checking both vigorous convection and strong rotational influence at the convective scale. Their local Rossby definition uses a spectral scale (equations 15–16). Their cyclonic Boussinesq result does not predict anticyclones in this model.
- Chan & Mayr (2013), the supplied `chan2013.pdf`: motivates a stratified, heated-from-below model with upper cooling and a rotation-dependent transition. Their background Rossby number uses the nonrotating reference velocity, not the instantaneous rms velocity of a vortex-bearing run. Their total-energy equation, diffusion and density structure differ from this pilot; retaining only the geometry is not reproducing their mechanism.
- [Verhoeven, Wiesehöfer & Stellmach (2015), Anelastic versus fully compressible turbulent Rayleigh–Bénard convection](https://arxiv.org/abs/1501.01237): supports testing weak-superadiabatic anelastic models, while warning against assuming that they are always more computationally efficient.
- [Verhoeven & Glatzmaier (2017), Validity of sound-proof approaches in rapidly rotating compressible convection](https://arxiv.org/abs/1701.04575): distinguishes marginal-stability difficulties from supercritical turbulence and identifies additional rotational validity considerations. It does not validate every low-Mach rotating calculation.

## Execution and reporting discipline

All scientific computation, diagnostics and rendering remain on Unity. Preserve prior source and outputs. Do not relax a stopping threshold to obtain a desired result. Report queue time separately from allocated compute time. A technical preflight pass is never described as physical validation.

Publish code, parameters, provenance and conclusions to the existing GitHub repository. Any new simulation-results presentation must include synchronized three-depth vorticity and passive-tracer videos with fixed, documented color scales. Status and design-review reports need not include movies.

**Next action: complete the existing-data audit and a defensible parameter/cost assessment. No new simulation is released by this document.**
