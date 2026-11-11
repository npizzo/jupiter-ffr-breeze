# Campaign review and next experiments — November 11, 2026

All simulation, diagnostic and movie computation remains on Unity. Finished publication assets are uploaded to this repository; the playable gallery is https://npizzo.github.io/jupiter-ffr-breeze/.

## Ferrel mechanism screening

The four six-hour 512² × 96 pilots completed and passed native numerical checks (job1409341). [Watch the three-depth vorticity/tracer comparisons](https://npizzo.github.io/jupiter-ffr-breeze/#ferrel-pilot). The [scientific review](../gallery/ferrel2026-pilot/review.md) distinguishes numerical validation from fidelity: weak lower-plane vertical motion peaks near two horizontal cells. The heated case has98.3% of its zonal-mean-removed vertical-velocity variance below four-cell radial wavelengths at6h. Upper-level heat cospectra are precision-indeterminate in the old Float32 absolute-temperature outputs. These pilots do not establish resolved convection, persistent Ferrel cells or an inverse cascade.

The48h unheated and heated cases completed (1409343 and1409344); boundary sensitivity1409345 is running at this record. Their reports are waiting for that comparison to finish. The existing1024² × 96 six-hour refinement1409342 remains pending for H100 resources; it was not cancelled or duplicated. Its queued start estimate can change.

Selected next batch: four fresh six-hour thermal-seed tests, preserving existing velocity initialization, jets, rotation, forcing, thermal balancing and numerical operators. Two512² zero-thermal-noise cases compare0/162W/m², and heated512²/1024² cases share the same analytic physical horizontal modes and the exact measured original512 per-layer variance. Added Float64 temperature-anomaly and vertical-velocity planes address diagnostic precision and near-bottom coverage. This is a horizontal-seed test; smoothing initial noise cannot prove intrinsic unstable scales are resolved. There is no vertical refinement or complete heating-by-seed factorial in this batch.

The new allocated initialization gate1422465 passed2613 checks on Unity. Production has not yet been submitted: first review existing fine-pilot output, then require its unchanged numerical gate and the new source/reference-table gate. [Definition, limits and launch instructions](../scripts/ferrel2026_seedtest_v1/README.md). Independent static review caught and corrected the launch-environment routing before deployment. Existing frozen drivers and jobs were untouched.

## Gamma-plane reference

All six1024² cases (polar/flat controls, seeds11/22/33) completed75Earth days. Reporting uses actual2D vorticity fields; the reference has no convection, vertical depths or passive tracer. The first-pair report1412700 completed. The new three-seed report1422447 is processing all completed pairs.

The original report stopped because the flat seed11 day1 total-PV snapshot/solver comparison differed by5.896×10⁻⁶, exceeding the original5×10⁻⁶ tolerance. Energy and relative-vorticity enstrophy checks pass at every saved time. Reporting v2 keeps the original thresholds and exposes the PV failure rather than labeling all checks passed. An independent native checkpoint audit1412654 reproduced final exports exactly and found native Fourier PV agreement near roundoff atday75; the exactday1 Float64 state was not retained, so that intermediate mismatch remains unexplained. Both solver and reconstructed PV series are retained. Verified vorticity movies may be used descriptively; tight PV-budget claims are not supported by a failed consistency audit.

The2048² L40S gate1402735 passed its numerical tests but estimated7.46768h for75days, exceeding the frozen throughput budget. Consequently dependent simulations1402736/1402737 were automatically cancelled without starting. This was a cost gate, not a user cancellation or demonstrated numerical instability.

At the user's request the unchanged2048² tests were resubmitted on H100s:

| Job | Purpose | Dependency |
|---|---|---|
|1422498|H100 numerical/throughput gate|none|
|1422516|75-day polar2048² seed11|successful1422498|
|1422518|75-day flat2048² seed11|successful1422516|

All are queued at this record. The same source, configurations, numerical tolerances and throughput criterion are retained; only the requested GPU type and job names change. If the H100 benchmark still exceeds budget, it must fail rather than weaken a gate. No1024 state is interpolated onto2048.

Longer gamma continuation code is prepared but inactive. A scientific extension target has not been submitted; actual paired movies, budgets, seed dependence and refinement evidence must inform that choice. Seventy-five days is early organization, not a mature vortex-crystal or inverse-cascade demonstration.
