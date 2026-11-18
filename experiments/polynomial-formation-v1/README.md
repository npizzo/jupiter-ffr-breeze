# What organizes the polynomial model's large vortices?

Authorized November 18, 2026. This campaign studies the original successful horizontal 2D polynomial model, without changing its forcing, saturation, viscosity, geometry, or initial state. All computation, analysis, and rendering stay on Unity. No new 3D convection run is part of this campaign.

## Matched formation control

The native Oceananigans centered momentum-advection term is disabled. Everything else is retained: periodic 2π square, 512² grid, f=+1, beta=0, viscosity 10⁻⁴, linear instability a=0.0002 in 12≤|k|≤16 excluding kx=0, and mean-subtracted quadratic/cubic feedback 0.002ζ²−0.0002ζ³. The original completed biased-plus trajectory is the advection-on control; it is not rerun. Both use `data/viscous2d-refine-v1-1848240/initial.jld2`.

Only one new trajectory is advanced from t=0 to 400 in five 80-unit segments. Each uses one CPU, 4 GB, a 36-minute scheduler limit and a 2050-second simulation limit, preserving a checkpoint on interruption. Total reserved production time is at most three CPU-hours. A failed segment holds its dependents; no automatic resubmission or extension. The original timestep ceiling and conservative CFL/reaction controls are retained. Vorticity frames are saved every0.5; native source work and the mathematical potential every0.25.

The existing cubic saturation is part of the expressly authorized polynomial hypothesis. No velocity drag, selective viscosity, hyperviscosity, beta, or new forcing is added. This is a diagnostic advection-off experiment, not a proposed physical fluid without advection.

Release requires native diffusion-only tendency verification, zero advection, source provenance, restart and timestep checks, potential descent, speed/allocation limits, and an actual production-driver output check. The frozen original code is not edited.

## Existing-output formation budget

Every saved half-unit field through t=400 is processed for all three original polynomial cases. Native-grid inverse curl and centered momentum flux reconstruct advection, which is checked against the actual Oceananigans operator at eight saved states spanning formation and controls. Source work is checked against saved scalar diagnostics. Absolute low-mode energy change is compared with separate advective, linear, quadratic, cubic, and viscous work. Half-unit and unit-interval quadrature are compared; the remaining residual is reported, never assigned to a physical source.

Advective transfer into several wavenumber ranges distinguishes nonlinear redistribution from direct low-mode polynomial input. Morphology uses the previously declared fixed threshold and periodic connectivity. Neither skewness nor an attractive movie establishes an inverse cascade or a unique LSV.

## Constrained potential

For the discrete symmetric Laplacian and band multiplier, the no-advection vorticity equation is projected gradient descent of

F = ν⟨|∇ζ|²⟩/2 − ⟨ζ Lζ⟩/2 − B⟨ζ³⟩/3 + C⟨ζ⁴⟩/4,

with the gradient term interpreted through the native discrete Laplacian. Mean subtraction is the projection enforcing zero circulation. The potential is not kinetic energy. With advection restored it need not decrease. Descent does not itself prove coarsening or a single localized state.

## Galerkin projection of the same equation

`galerkin.py` derives mode couplings directly from the original 512-grid centered momentum-flux symbols and pointwise quadratic/cubic polynomial. It preserves conjugate pairs and omits the zero mode. The linear coefficients, viscosity, nonlinear coefficients, and phase information are unchanged. Products are evaluated by exact finite convolutions; cubic convolution retains all intermediate pair modes, even when those intermediates are absent from the state. Candidate |k|≤64 support avoids aliasing at the original 512 grid.

Selection uses only t≤200: retain all 0<|k|≤2 modes and three bounded sets of 16,32,64 conjugate pairs, divided between the forcing band and other modes ranked by their absolute contribution to quadratic low-mode work. t=210–400 is held out. Low-mode total, quadratic, and cubic tendency errors must each be ≤20% in an energy-weighted RMS norm. The threshold is a declared operational screen, not a theorem.

Only the smallest passing projection is integrated from the projected original initial condition to t=400, with an independent low-energy-history prediction screen. No parameters are fitted. Failure of all three sizes stops this line; it does not automatically launch more modes or a new PDE. The entire ROM job is limited to15 minutes on one CPU. Time holdout is not an independent-seed validation.

## Convection connection: a separate gate

The original capped f-plane sign audit closes its chosen native numerical budget, but its physical source partition failed: unassigned corrections were about0.72 and5.91 times the convective central-third-moment contribution for +f and −f. This campaign audits that evidence and does not fit a signed polynomial coefficient from it. A physical closure requires a sufficiently resolved source partition and a prediction of signed, scale-dependent convective work in a held-out interval or case. An f flip imposed in a phenomenological source is not independent physical evidence. No additional GPU simulation is released by this campaign.

## Output

The existing-output analysis, ROM prediction tests, provenance, and matched advection-on/off vorticity movie belong on the GitHub Pages gallery. This is a genuinely horizontal 2D model: one plane and no tracer or fabricated depth panels. Movies use a single symmetric asinh range fixed across time and cases. Pending and failed physical-prediction gates must remain visible.
