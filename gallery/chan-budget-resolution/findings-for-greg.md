# Breeze compressible energy-budget deficit: short resolution and timestep checks

**24 November 2026 — for discussion with Greg Wagner**

**Finding:** a reconstructed total-energy budget in our dry compressible convection experiment has a negative residual. A short, controlled test found essentially no improvement from a smaller timestep, no improvement from a 1.5× refinement in all three directions, and negligible error from the heat-input sampling cadence. This narrows the diagnosis, but does not establish whether the residual is a true numerical energy sink, an inappropriate discrete energy diagnostic, or an omitted contribution.

- [Results and movies on GitHub Pages](https://npizzo.github.io/jupiter-ffr-breeze/gallery/chan-budget-resolution/)
- [Exact tested source](https://github.com/npizzo/jupiter-ffr-breeze/tree/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_v3)
- [Parent clean compressible/anelastic configuration](https://github.com/npizzo/jupiter-ffr-breeze/tree/d84454f7caddd9645cce7721b91234f14d514c9e/scripts/chan_clean_pair_v1)

## Scientific purpose

We are testing whether otherwise comparable compressible and anelastic rotating convection can favor different signs of large-scale vortices. The present pair removes beta and imposed jets to simplify that comparison. The energy audit is intended to establish how much confidence to put in that interpretation; it is not itself a test of vortex-sign selection.

## Configuration actually used

- Breeze 0.15.0, pinned repository revision `0479d9c8ce006adcb1e3314c324b84980ad4871b` (the launch receipt also fingerprints the installed source files).
- Three-dimensional dry ideal-gas atmosphere, potential-temperature formulation (`ρθ` prognostic); not a prognostic conservative total-energy formulation.
- `CompressibleDynamics(SplitExplicitTimeDiscretization(damping=NoDivergenceDamping()))`; no upper sponge.
- WENO9 momentum and potential-temperature advection; `closure=nothing`. Thus there is no explicit SGS closure, but advection still introduces numerical dissipation.
- Box `4 × 4 × 1`; periodic x/y, impermeable stress-free horizontal boundaries.
- Rotation rate Ω=3 at 67.5°N, including both vertical and horizontal Coriolis components; f_z≈5.5433. Beta=0, no jet, shear forcing, or momentum relaxation.
- R=1, c_p=2.5, c_v=1.5, γ=5/3, g=12; energy reference temperature=0.
- Initial reference profile: lower polytropic index 1.5 and upper index 24, interface z=0.78. The native discrete hydrostatic reference is constructed on each grid. The parent starts from rest with weak thermal noise; these tests branch from its developed state at t=200.
- Uniform bottom energy flux F_b=0.03125. Horizontal area A=16, so bottom power is AF_b=0.5.
- Weak conductive coefficient κ=0.0004 F_b/(g/c_p), with a conductive top boundary referenced to T=1.
- Local cap source is `Q_cap = -ρ (T-1) S(z)/20`, where S is cubic smoothstep from zero at z=0.78 to one at z=1. No c_v multiplier is applied. This term cools T>1 and heats T<1; it is not a high-temperature clipping rule.

All times and energies below are in the simulation's nondimensional units. These are evolving convective states, not demonstrated thermal equilibria. All computation and rendering were performed on Unity A100 GPUs / CPU nodes; active long runs were not modified.

## What the residual measures

We reconstruct

\[
E(t)=\int_V\rho(c_vT+gz)\,dV+K(t).
\]

K uses face velocities with arithmetic density interpolation to their staggered locations. Horizontal periodic faces are counted once; vertical boundary velocities are zero. This is a physically motivated discrete energy diagnostic, **not yet verified to be the exact invariant of the discrete Breeze equations**.

The imposed net power is evaluated from the same local cooling law and conductive boundary flux used in the driver:

\[
P(t)=16\,[F_b-F_t(t)]-\int_V\rho(T-1)S(z)/20\,dV,
\]

with

\[
F_t=2\kappa\,[\langle T_{\rm top\ cell}\rangle-1]/\Delta z.
\]

For each short branch, the energy ledger starts at its own actual initial state:

\[
R(t)=E(t)-E(t_0)-\int_{t_0}^{t}P(s)\,ds.
\]

Negative R means less reconstructed energy accumulated than the prescribed net heat input predicts. Reported percentages divide R by **gross bottom heat**, `0.5 × (t-t0)`, not by total atmospheric energy. Prescribed cap cooling and top exchange have already been included.

Power is integrated by a trapezoidal rule at **every outer timestep**. An independent ledger samples it every one time unit, reproducing the parent's budget cadence. Energies are recorded every 0.1 time unit. This tests sampling frequency; it is not an exact Runge–Kutta-stage-integrated source ledger.

## Controlled tests and results

All three ran from t=200 to t=210, with unchanged physics and the same donor state.

| Test | Grid | Maximum outer Δt | Incremental R at t=210 | R / gross bottom heat |
|---|---|---:|---:|---:|
| Control | 256 × 256 × 80 | 1/120 | −1.2667467 | **−25.3349%** |
| Refined | 384 × 384 × 120 | 1/180 | −1.2854029 | **−25.7081%** |
| Smaller timestep | 256 × 256 × 80 | 1/180 | −1.2667810 | **−25.3356%** |

Gross bottom heat over the interval is 5. The control's net prescribed heat is 4.9631513, while its reconstructed energy increases by about 3.6964045. Their difference is −1.2667467.

The finer run occasionally used slightly smaller timesteps than its cap (reported minimum about 0.0054753 versus the cap 0.0055556), so the two small-timestep trajectories are matched in maximum timestep, not guaranteed identical in every step.

**Heat sampling:** the final once-per-unit versus every-step heat discrepancy is −4.35×10⁻⁶ energy units for the control and +5.75×10⁻⁵ for the refined run. Both are negligible compared with a residual of order 1.27. Sampling frequency does not explain the deficit.

**Transient sensitivity:** fitting the residual slope only over t=202–210 gives deficits equivalent to 25.2780%, 27.5981%, and 25.2788% of bottom power for control, refined, and smaller-timestep cases respectively. Excluding the initial refinement adjustment therefore does not reveal an improvement.

The earlier quoted 16% was the parent's cumulative residual at t≈134. The table above measures newly accumulated error over t=200–210. These are different averaging intervals, so they should not be presented as a direct change caused by refinement.

## Transfer, validation and computational cost

- A private native t=200 checkpoint was preserved; the long run continued untouched.
- Same-grid reinstallation reproduces prognostic fields, temperature and velocity exactly. Its initial reconstructed energy agrees with the parent's recorded t=200 energy to about 9×10⁻¹³.
- Fine-grid transfer conserves integrated density, potential-temperature density, dye density and horizontal momentum dual-volume averages to roundoff, using limited piecewise-linear reconstruction. Vertical face momentum uses linear interpolation with exactly zero boundary values.
- Transfer changed total energy by **−0.0956731**, and kinetic energy by **−3.8663%**. The refined ledger is reset after that transfer: the quoted deficit does not include the instantaneous energy jump. Subsequent adjustment remains a limitation, addressed partly by the late-window comparison.
- All three completed with positive density/temperature, peak sampled Mach below 0.20, and relative mass errors ≤3.4×10⁻¹⁶. These checks do not prove energy fidelity.
- CPU and GPU output/restart checks passed; restarted fields and accumulated heat matched uninterrupted tests exactly.
- Completed allocations took 13m22s (control), 28m41s (refined), and 13m52s (smaller timestep), including startup and validation. Integration/output portions were about 5.6, 20.9 and 6.1 minutes.
- The previous submission failed in the small GPU validation because `Statistics.mean` on a nested GPU slice attempted scalar access. That diagnostic implementation was repaired using a preallocated GPU buffer and sum. Those failed allocations consumed about 14 GPU-minutes total and produced no full-size comparison data. The successful rerun first passed one GPU validation before releasing the remaining two cases.

## What this establishes—and what remains open

The tested timestep reduction has essentially no effect. A 1.5× spatial refinement does not reduce the residual. Therefore these results do **not support simply increasing resolution as a remedy**. They do not constitute an asymptotic convergence study or rule out all resolution-dependent contributions.

The residual's cause remains unidentified. In particular, we have not demonstrated that it is WENO kinetic-energy dissipation, heat being silently clipped, or an incorrect heating source. Possible explanations include the thermodynamic formulation's discrete energy properties, mismatched discrete pressure/gravity work, the reconstructed energy definition, or an omitted contribution. Numerical stability and mass conservation do not distinguish these.

The analogous compressible energy expression is not an established conserved invariant for the anelastic model. Its residual must not be used as a directly comparable numerical energy-loss measure.

## Specific questions for Greg

1. For this dry `ρθ` formulation with split-explicit compressible dynamics, which discrete total-energy diagnostic should close against a prescribed `ρE` forcing? Is the face-density-weighted K plus cell-centered internal/gravitational energy above appropriate?
2. The native source conversion in our preflight is `G_ρθ = Q × θ/(c_p T)` and agrees with the supplied forcing. Is that the intended volumetric-energy conversion for this setup, including our nondimensional constants? Is there a stage-integrated heating diagnostic we should use instead of reconstructing source work afterward?
3. Which native tendencies best expose energy exchange in pressure work, gravitational work, acoustic substeps, and thermal/momentum advection? Our next useful check would be a short term-by-term accounting from an existing developed state.
4. We retain the original `MomentumCoriolis` wrapper: it passes `ρu,ρv,ρw` under the `u,v,w` aliases to `ConstantCartesianCoriolis`. Does that preserve the intended discrete kinetic-energy property with spatially varying face densities, or require a different operator/diagnostic pairing?

The immediate request is help identifying the appropriate discrete energy balance and any missing source/work terms. We have not changed the thermodynamic formulation, cooling, closure, or long-run configuration in response to this result.

## Reproducibility links

- [Budget calculation](https://github.com/npizzo/jupiter-ffr-breeze/blob/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_v3/diagnostics.jl)
- [Ledger and output callbacks](https://github.com/npizzo/jupiter-ffr-breeze/blob/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_v3/run.jl)
- [Model and explicit heat/cooling source](https://github.com/npizzo/jupiter-ffr-breeze/blob/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_v3/model.jl)
- [State transfer](https://github.com/npizzo/jupiter-ffr-breeze/blob/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_v3/transfer.jl)
- [Analysis script](https://github.com/npizzo/jupiter-ffr-breeze/blob/b78fcc68289c2b9edb40d6f9773fe092292281ac/scripts/chan_budget_resolution_analysis_v2/analyze.py)
- [Machine-readable results](https://npizzo.github.io/jupiter-ffr-breeze/gallery/chan-budget-resolution/summary.json)

Successful short jobs: control `2221768_0`, refined `2222246_1`, smaller timestep `2222246_2`. The original t=200 donor checkpoint SHA-256 is `4031aed6757ce94eb271af4890452ffe5040e6af4e60ded28783d70c3a0c54d5`.
