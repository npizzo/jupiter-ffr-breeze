# Testing a polynomial closure against existing heated Breeze flow

Analysis only, authorized November 17. No simulation is advanced, no checkpoint is modified, and the existing cap, heating, rotation and numerical setup remain unchanged. The target is the mature dry anelastic f-plane run at 55°N: 256×256×64 cells in a 6000×6000×120 km channel, periodic zonally, with 162 W/m² bottom heating. This is separate from beta-plane Ferrel-cell work.

## Equation and known attribution limitation

With reference density depending on height only, U=<u_h>_rho and Z=curl_h U. The centered resolved convective source is C=-curl_h(div_h R), where R_ij=<u'_i u'_j>_rho. Saved velocities and stresses are independently checked against two saved full 3D states. Heating influences C through the 3D correlations. The vertically integrated constant-f stretching contribution vanishes for impermeable boundaries; it is not replaced by a fitted source.

The six-hour +f/−f diagnostic restarts retain native every-stage budgets. Their exact accounting separately retains mean transport, mean-flux/transport correction, centered convective stress, native-minus-centered advection, Coriolis, diffusion and horizontal forcing. Earlier work found a substantial native/centered difference. It is not relabeled as convection or discarded. The negative-f run inherits the same positive-f mature state and includes adjustment; it is not a mirrored independent realization.

## Predeclared prediction test

- Primary input: all 81 three-hourly states on days 90–100 of the heated run.
- Fit only days 90 <= day < 95; evaluate on days 95 <= day <= 100. No random spatial split and no independent-gridpoint significance claim.
- Fit on the central 80% of the meridional channel. Report held-out interior, full-channel and wall errors separately.
- Normalize z=Z/|f| and q=C/|f|². Coefficients are learned from Breeze, not taken from the unrelated nondimensional 2D coefficient values.
- Compare a linear baseline (z and h² Laplacian z), unconstrained quadratic/cubic feedback, the same feedback constrained to cyclonic bias and cubic saturation, and a fixed extended model including strain, vorticity gradients and convective intensity. This flexible linear baseline does not derive the 2D band-limited forcing.
- Polynomial features are sign(f) P0(z²) and -P0(z³). P0 removes the full-channel mean of each feature. The measured target mean is retained: channel circulation is not presumed zero, and a source-mean mismatch remains an error.
- Primary analysis uses unfiltered centered fields. One fixed sensitivity uses a Gaussian sigma of four grid cells (93.75 km), with periodic x and parity-aware channel y extension. Filtered convective stress stays distinct from any barotropic subfilter stress. This is not a resolution test.
- Metrics: held-out source R²/normalized RMSE, conditional drift, central-third-moment source and work in large channel modes. The latter use orthogonal channel streamfunctions with matched centered derivative symbols and horizontal wavenumber magnitude <=2(2π/Lx). They are source contributions, not a complete ten-day energy budget.
- Audit ten-minute trapezoidal source samples against saved every-stage integrals in both six-hour sign restarts. Keep temporal sampling error and native-advection reconstruction error separate. Compare transferred fits against centered C and the explicitly labeled C-plus-numerical-remainder effective source.

The fixed feature library is documented in `scripts/fplane_closure_v1/analyze.py`. Fits use scaled small linear systems; positive quadratic and damping-cubic coefficients are imposed only in the explicitly constrained model. Negative unconstrained coefficients are evidence against the proposed signs, not silently corrected.

## Validation, cost and publication

One Unity CPU, 8 GB memory, 30-minute hard ceiling, no requeue. Checks cover analytic centered curl, reflection of all features under f reversal, zero-mean projection, synthetic coefficient recovery, the quadratic third-moment identity, low-mode orthogonality and momentum/curl work agreement, independent 3D density weighting, and the existing native full-domain budget gate. Computation and figure rendering stay on Unity. Raw fields are not downloaded.

No successful fit alone establishes a physical closure: prediction must survive held-out times, capture the source contributions of interest, and be assessed against the numerical reconstruction error. A failed fit is a valid scientific result and does not trigger a new simulation or parameter search. Published results will embed the existing same-run three-depth vorticity/tracer movie and link the existing sign-test movie; no redundant movie rendering is required.
