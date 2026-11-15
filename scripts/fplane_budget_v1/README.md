# Dense native-vorticity budget and horizontal refinement

This isolated campaign responds to the large residual and resolution sensitivity in the initial f-plane sign-selection analysis. It does not change any existing simulation source or checkpoint.

Two 18-hour heated restarts, 512²×64 and 1024²×64, begin from the same actual day-90 f-plane state. Domain, bottom heat flux, cap, thermal controller, physical mixing, reference density and 15-second timestep ceiling are held fixed. Conservative commuting horizontal refinement is applied once or twice to the original 256² staggered state. The restriction error, native mass divergence, scalar conservation and initial energy change are checked. The last six hours, after twelve hours of adjustment, form the initial comparison window; this is not a claim of equilibrated resolution convergence.

## What changes in the diagnostic

The new observer evaluates the actual pinned Breeze horizontal momentum operators: explicit advection, Coriolis, diffusion and external forcing. It reduces native staggered momentum and its tendencies vertically on the GPU, divides by the mean reference density, and applies the same linear cell-centered curl used for the existing depth-mean Z diagnostic. A separate pressure-curl check tests the projection cancellation in the interior.

Every SSPRK3 stage contributes with weights 1/6, 1/6 and 2/3. Each independently observed state increment must match the integrated native tendency, with a scale-aware roundoff allowance. The exact finite-step change in the raw third moment uses

    Δ(Z³) = (Z_end² + Z_end Z_start + Z_start²) ΔZ.

This avoids replacing a finite change with infrequently sampled 3 Z² ∂tZ. It is an accounting identity for the discrete model, not evidence that Z³ is conserved. The actual WENO advection term must not be relabeled as solely physical transport or numerical dissipation.

The V/H continuum reconstruction is retained as a separate diagnostic, sampled every ten minutes rather than daily. Comparing 10-, 20- and 60-minute subsampling against the every-stage native integral measures sampling sensitivity without confusing it with native-budget closure. A successful native balance does not guarantee that the centered V/H decomposition is accurate.

## Validation and execution

1. A 64²×64 CPU observer test compares all five final native fields bit-for-bit with an otherwise identical run. Its seed wavelength is 2400 km so it satisfies the existing 16-cell initializer guard; this is only a diagnostic test, not the production initial condition. This passed in job 1767120. Independent saved-field reader/reduction validation also passed in job 1767470: maximum disagreement in Z was 1.22×10⁻¹⁹ s⁻¹, and the small-test native balance residual was 3.31×10⁻¹² relative to the measured increment scale. These are diagnostic validation numbers, not full-resolution physical results.
2. A GPU test repeats the unchanged-state comparison. Full 512² and 1024² mature-state pilots then test memory, remapping, every-stage closure, temporary-field reduction and runtime. Production requires an estimated runtime below 14 hours within a 16-hour allocation.
3. An independent Python gate checks cumulative native balance, the third-moment identity and agreement of Z with the existing independently reduced three-dimensional calculation.
4. Only successful gates release the two production cases. Their dependent review compares exact budgets, output cadence and common-grid sensitivity, and renders a synchronized four-row, three-column vorticity/tracer movie.

The current chain is preflight 1767196, production array 1767197 (one GPU at a time), review 1767284. All computation and rendering stay on Unity. GPU allocation requests at least 80 GB VRAM. The earlier CPU development checks failed before science production (a syntax error, then an under-resolved small-test seed); those were corrected. An obsolete pending GPU pilot was cancelled before it executed. Existing science jobs were not cancelled.

## Bounded storage

A temporary Float64 3D export is written every ten minutes, reduced on Unity into source maps and scalar diagnostics, checked, then retired. If reduction fails, the export remains and the job stops. Native cumulative maps, reduced diagnostics, three-depth movie fields and restart checkpoints are retained. No raw field downloads are part of publication.

## Decisions after review

The user authorized monitoring and longer runs if the checks justify them. First verify numerical closure and measured cost/storage. If the centered reconstruction changes materially with sampling, add a short finer-cadence probe. If the fine state is still adjusting or statistics differ materially across grids, longer matched continuations may be useful, but should not be called convergence. Vertical resolution, cap/wall sensitivity and independent-seed robustness remain separate physical checks; horizontal refinement alone does not resolve them. Central and nonzonal skewness and coherent-structure evolution must accompany raw third-moment statistics before claiming reversal or destruction of a vortex.
