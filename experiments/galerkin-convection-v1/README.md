# Coupled nonlinear Galerkin convection

This bounded experiment tests whether rotating convection can generate a preferred sign of depth-averaged vorticity when the barotropic flow, three-dimensional fluctuations and mean thermal profile all evolve together. The forcing of the barotropic mode is the projected Reynolds stress from convection. No sign-dependent vorticity reaction is inserted.

The [derivation](DERIVATION.md) specifies the equations, exact transfer budgets, symmetry, truncation caveats and how a later adaptation to the original capped Breeze setup would differ. This first experiment is a **standard fixed-temperature Boussinesq reference**, with free-slip impermeable walls and no cap. It is not a quantitative reproduction of the earlier Breeze runs.

## Assumptions fixed before execution

- Three-dimensional incompressible rotating Boussinesq equations; constant \(f=+1\) or \(-1\), with \(H=1\). Positive vorticity is cyclonic at positive \(f\).
- Periodic horizontal square; fixed-temperature, free-slip top and bottom. Heating enters through an unstable conductive profile and the evolving temperature equation.
- Ordinary Laplacian viscosity and thermal diffusion, initially equal. No beta, drag, hyperviscosity, moisture, stochastic stirring or polynomial feedback.
- All nonlinear couplings among retained velocity and temperature modes, including reciprocal transfer to the barotropic flow and evolution of the mean temperature.
- Initial mode cutoff \((K_h,N_v)=(4,2)\); nested extension \((6,3)\). Horizontal-mean velocity and temperature extend to \(2N_v\) to retain required balance terms. All horizontal phases and conjugates are kept.
- Two independent positive-\(f\) seeds, one exactly reflected negative-\(f\) counterpart, one matched zero-imposed-gradient control, and one extension using the same physical parameters and embedded initial state.
- \(k_0=20\), \(L/H=2\pi/20\), \(\nu=\kappa=10^{-5}\). Heating is \(G=2\min G_c\) over the retained base convective modes, evaluated on Unity and fixed for the extension.
- Velocity RMS \(10^{-4}\), buoyancy RMS \(10^{-5}\); seed identifiers 20261118 and 20261119. Initially zero barotropic flow and horizontal-mean fields. Newly added refinement modes start at zero.
- Duration \(|f_0|t=200\), timestep 0.02, with timestep validation and a projected total integration budget of 1800 CPU seconds. Actual coefficients and measured cost are recorded before release. A negative result does not trigger automatic longer runs.
- All model integration, diagnostics and rendering run on Unity. Local work is limited to source and publication editing.

## Required release checks

1. Reality, boundary parity, incompressibility and zero net vertical circulation are satisfied by the basis and projected tendencies.
2. Instantaneous kinetic energy, thermal variance and physical potential moment close independently. The potential moment is integrated on \(0\leq z\leq1\), not incorrectly over the doubled Fourier domain.
3. Barotropic-to-fluctuation energy transfer is reciprocal. The resolved depth-averaged vorticity equation closes using the actual retained-mode projection.
4. The independently calculated linear eigenvalues agree with the implemented operator, including the selected convective onset.
5. Reflection with \(f\to-f\) commutes with the full nonlinear tendency; short matched integrations preserve that symmetry and converge when the timestep is reduced.
6. A measured short Unity timing establishes a finite wall-clock budget before the five nonlinear cases are released. Failed balance, finite-state or cost checks block that release.

## Comparison design

| Case | Rotation | Background gradient | Initial condition | Cutoff |
|---|---|---|---|---|
| Convective seed A | \(f=+1\) | \(2\min G_c\) | Seed 20261118 | \((4,2)\) |
| Convective seed B | \(f=+1\) | Same | Seed 20261119 | \((4,2)\) |
| Reflected sign control | \(f=-1\) | Same | Exact spatial reflection of seed A | \((4,2)\) |
| Zero-gradient control | \(f=+1\) | \(G=0\) | Exact seed A state | \((4,2)\) |
| Nested extension | \(f=+1\) | Same as seed A | Exact embedded seed A state | \((6,3)\) |

This table specifies the planned cases; execution and completion require separate recorded status. The zero-gradient control keeps the initial buoyancy anomaly and its force, so it has a finite initial thermal reservoir. It is not a pure velocity-only decay experiment.

## How to interpret the result

The useful quantity is the signed stress forcing of the depth-averaged vorticity and its associated power, not a visually impressive vortex alone. Report absolute barotropic energy, convective conversion, reciprocal exchange, skewness and the finite-truncation contribution to the skewness budget. A high barotropic energy fraction by itself does not demonstrate an inverse transfer.

Finite Fourier-Galerkin advection preserves relevant quadratic balances but generally does not preserve the cubic vorticity moment. The projected-advection contribution to skewness must therefore be shown separately from convective production. A preference that disappears when additional modes are retained is not a reliable minimal mechanism.

Two seeds and one extension can support a preliminary mechanism test, not a statistical claim about Jupiter. The narrow box, \(L/H\simeq0.314\), contains only a few critical convective wavelengths. It is a finite-mode hypothesis test, not LSV DNS or a resolved Jovian layer; the limited scale separation prevents a broad inverse-cascade claim. If the dynamics reduce to steady rolls, lack sustained transfer, or fail the extension, the appropriate outcome is a documented negative test and a revision of the reduction.

When outputs are available, the gallery should show actual vertical vorticity at three heights in one synchronized frame, with fixed physical color scales and the recorded dimensionless time. This model has no passive tracer; any displayed thermal field must be labeled buoyancy perturbation rather than tracer concentration. A planned experiment or a queued job is not a published result.

## Primary reference

[Rademacher and Welter, *Energetic consistency and heat transport in Fourier-Galerkin truncations of free slip 3D rotating convection*](https://arxiv.org/abs/2505.07678) supplies the consistency framework for selecting vertical mean modes. The present anisotropic full-phase implementation uses a separately derived projection of the PDEs and must pass its own tests. It does not assume that a truncation validated for heat transport automatically captures large-vortex sign selection.
