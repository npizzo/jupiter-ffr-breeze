# Bounded f-plane vorticity-sign experiment

## Question and unchanged physics

Identify which terms select the sign of the density-weighted depth-mean vertical vorticity Z in the existing heated, capped Breeze anelastic channel. This is not a new convection model, beta experiment, or literature-reproduction campaign.

Both cases restart the actual day-90 heated state without remapping: 256×256×64, 6000×6000×120 km, original 30 s timestep ceiling, 162 W/m² bottom heating, the existing reference density, stable cap, heat controller, viscosity/diffusivity and vertical sponge. Both run exactly six additional Earth hours. `plus` retains f; `minus` changes only f to its negative. All initial prognostic fields and the initial pressure guess are identical. The negative-f trajectory includes adjustment of a mature positive-f state; it is not an unbiased long-time sign-selection ensemble. No beta or imposed jet is allowed by the driver.

## Independent implementation checks

A small CPU comparison requires baseline and instrumented dynamics to give bitwise-identical fields. Analytic stress/reflection tests check derivative conventions. A native-grid GPU preflight repeats baseline versus instrumented dynamics for ten minutes, plus a separate ten-minute reflected-state experiment. The reflection is x→−x, u→−u, v,w and scalars reflected with unchanged sign, and f→−f. Native staggered face locations are respected. This tests symmetry, not an independent physical explanation.

Full-channel native Z increments are compared against the time integral of actual Breeze momentum tendencies at every SSPRK3 stage. Explicit diffusion and horizontal forcing are included. Pressure-gradient curl, column mass divergence and complete-channel closure are checked, including boundary strips. The native closure tolerance is 1e−8 of the largest increment/term scale plus a specified Float64 roundoff allowance; no source is fitted to force this identity. An independent saved-field reader checks Z. Reflection tolerances are 1e−6 for normalized fields and 1e−5 for accumulated source maps. CPU/GPU checks must pass before production can execute; the queued pair has an after-success dependency and checks the frozen validation record.

## Physical attribution

For U=〈u_h〉ρ and R_ij=〈u'_i u'_j〉ρ, the continuum convective source is C=−curl_h(div_h R). Reference density depends only on z. Impermeable boundaries imply div_h U=0, so the explicit constant-f contribution vanishes in the continuum depth-mean vorticity equation. f can still influence the 3D correlations R.

Native advection is partitioned into four explicitly named pieces:

1. Centered mean-flow transport −U·grad_h Z.
2. Centered mean momentum-flux curl minus that transport; a discrete product/centering correction, not a new physical source.
3. Centered convective stress divergence curl C, evaluated from density-weighted 3D products reduced on the GPU at every integration stage.
4. Actual native advection minus the centered total momentum-flux curl. This is a measured reconstruction/discretization difference, not identified automatically as physical viscosity or numerical dissipation.

Native Coriolis, explicit diffusion (including boundary stresses) and horizontal forcing complete the accounting. Normal momentum boundary updates are respected; pressure curl is tested without zeroing the normal pressure-gradient faces. The full-domain budget is required, not just an interior-only gate. Full channel, central 80% and wall strips retain separate contributions; transport across the cropped interior is not discarded. The full stress term is evaluated before making any statement about a local stretching contribution.

For each actual step the raw third-moment contributions use the exact finite-increment identity Δ(Z³)=(Z_end²+Z_end Z_start+Z_start²)ΔZ. Central third moments use centered end states and subtract each region's mean tendency. This avoids a changing regional mean masquerading as skewness production. These diagnostics do not assert a conserved third invariant.

Conditional mean tendencies E[S_i|Z] are accumulated at every RK stage in fixed symmetric bins of Z/|f|, for full/interior/wall regions. Positive/negative matched bins give the even component of C. Raw and central third-moment contributions distinguish convective production from explicit diffusion, transport and unassigned corrections.

## Bounds and stop rules

- Small CPU validation: 30-minute allocation, four CPU cores, no GPU.
- Native-grid GPU preflight: one GPU, 45-minute hard limit.
- Two production cases: one GPU each, 90-minute hard limit each, serialized.
- Total requested GPU ceiling for this initial chain: 3.75 GPU-hours. No requeue or automatic extension. Warm-run extrapolation must predict each six-hour case in under one wall hour before submission.
- Postprocessing/movie review: CPU only, 60-minute hard limit.
- Physical attribution is screened separately: unassigned mean-flux and advection corrections must contribute less than 10% of the claimed convective central-third-moment production, whose sign must persist across six one-hour blocks. This is an operational screen, not a statistical confidence interval. Otherwise report unresolved attribution and stop; do not automatically refine or extend.
- All computation and media validation stay on Unity. Raw fields remain there. Temporary start/end exports are reduced and removed only after successful validation; compact diagnostics, movie fields and final checkpoint remain.

## Output

Integrated Z and raw/central-third-moment budgets, separate wall/interior contributions, conditional drift and its even part, one-hour source blocks, and the independent symmetry-check record. The paired movie is four rows by three depths: +f vorticity/tracer, then −f vorticity/tracer. Vorticity is normalized by |f| so a physical sign reversal remains visible; positive means counterclockwise in x-east, y-north, z-up coordinates. Tracer uses the existing fixed absolute-concentration contrast scale. No per-frame renormalization.

The earlier 512²/1024² budget campaign remains held and is not released by this experiment.
