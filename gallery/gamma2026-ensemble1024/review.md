# Review of the completed 1024² gamma-plane ensemble

Recommendation: finish the already queued, unchanged 2048² matched pair before selecting a long extension. The 1024² results contain a repeatable descriptive effect of polar PV geometry, but the trap boundary and numerical scale dependence remain important alternative explanations for some of the visual organization. More elapsed time cannot resolve those uncertainties.

This review reads the finished Unity publication package at `work/finished-publication/gamma2026-ensemble1024`, including all three contact sheets, the diagnostic/radial-profile figures, JSON summaries and README. No raw fields were transferred for this review; no scientific computation or rendering was performed locally. The scheduler job identifiers below were supplied by the parent task and were not independently queried here.

## What the completed runs establish

All six 1024² cases have reached 75 Earth days, with 76 daily saved frames. The publication manifest records full Unity decoding of the three synchronized, two-dimensional polar/control movies. These are unforced barotropic QG calculations: they have neither convection nor physical depth levels, and they do not yet test Breeze's convective effect on polar organization.

The matched-pair changes in the sign-weighted radial statistic are positive for all three seeds:

| Seed | Polar change / Lγ | Flat change / Lγ | Polar minus flat change / Lγ |
|---|---:|---:|---:|
| 11 | −0.62043 | −1.13847 | +0.51804 |
| 22 | −0.34294 | −0.95470 | +0.61176 |
| 33 | −0.51050 | −0.96559 | +0.45509 |

The already-computed ensemble mean is +0.52830 Lγ, with sample standard deviation 0.07884 Lγ and range 0.45509–0.61176 Lγ. This is a consistent effect across three realizations, not a significance test or an established population distribution.

The sign is crucial. The statistic is the negative-vorticity mean radius minus the positive-vorticity mean radius inside r < 5 Lγ. Both individual changes are negative. The supported statement is that the polar cases undergo a **less negative change than their controls**. It would be incorrect to describe this number as positive absolute radial segregation or proof that the bulk positive circulation has migrated inward. Positive-vorticity mean radii at day 75 remain about 3.6–3.8 Lγ and are not systematically smaller than the controls.

There is complementary descriptive evidence: positive circulation inside r < Lγ comprises about 10.1–14.2% of in-trap positive circulation in the polar cases, versus 3.0–7.1% in the controls. The contact sheets show several compact positive-vorticity structures near the center of each polar realization and negative structures farther out. These are field-based observations, not object tracking, a demonstrated mature vortex crystal, or an inverse-cascade measurement.

## Main adversarial concerns

1. **The common trap creates a strong ring.** A conspicuous positive ring just inside r = 5 Lγ and an opposing exterior band occur in both geometries. The frozen background PV has a jump at this radius; the trap is not an impermeable wall. The whole-trap sign-weighted statistic includes this diffuse annular circulation as well as compact vortices. Therefore its sign and changes cannot be read as the motion of individual cyclones or anticyclones.

2. **The integration region is not closed.** At day 75 the reported fraction of relative enstrophy outside the trap is 10.5–22.6% for the polar cases and 40.1–47.3% for the flat controls. The matched difference may combine interior rearrangement, boundary production/deformation and export. Outside-trap enstrophy is a descriptive vorticity measure, not an energy-leakage budget. A later extension should also watch approach to the periodic domain boundaries.

3. **The initial relaxation is substantial.** Energy retention is about 93.64–93.73% at day 75, with most loss early in the displayed histories. Relative enstrophy retention is about 5.97–6.96%. Relative enstrophy is not conserved in this spatially varying background PV, so its entire decline cannot be assigned to dissipation. Nevertheless, the filtered small-scale evolution and sharply varying trap make refinement a high-value check. The nominal initial inverse wavenumber is 200 km: the corresponding full wavelength is about 1,257 km, roughly 10.7 cells on the present grid and twice that on the fine grid. The taper broadens that initial band.

4. **The PV comparison is not completely validated.** Saved-field energy and relative-enstrophy reconstructions pass the unchanged 5×10⁻⁶ criterion at every saved time. Total-PV reconstruction fails at day 1 for flat seeds 11 and 33, at 5.896×10⁻⁶ and 6.626×10⁻⁶. These are unresolved diagnostic discrepancies, distinct from the much larger evolution of total-PV enstrophy. They must remain visible in publication, without asserting that all checks passed or assuming Float32 rounding explains them.

5. **Atmospheric interpretation remains separate.** This is a controlled mathematical benchmark with different code/filter/seeds from the historical calculation. Strong localized vorticity and velocity maxima in the summaries warrant a local-Rossby/scale check before identifying them with quantitative Jovian winds. Three 75-day realizations do not establish mature long-lived polygonal order, global Jupiter fidelity or any convective mechanism.

## Next experiment and decision criteria

The queued H100 gate 1422498 followed by polar 1422516 and flat 1422518 is the appropriate next step. Keep the frozen core, initialization, configuration and 75-day duration. The older 2048² L40S gate passed its numerical checks but projected 7.46768 hours, exceeding its allowed budget; this was a throughput failure, not a reason to weaken numerical gates. Require the H100 gate's actual runtime/headroom result before its dependent production jobs can proceed.

Compare 1024² and 2048² seed 11 at common physical times using identical physical masks and normalization. Prioritize the early energy/PV histories, whole-trap and core statistics, exterior enstrophy fraction, radial profiles and compact-vortex scale distribution. Exact pointwise agreement at day 75 is not an appropriate expectation for chaotic evolution; differences in statistical organization and their size relative to the three-seed spread are the useful evidence. One fine pair tests sensitivity, not convergence.

Before declaring the segregation interpretation robust, add allocated-Unity diagnostics that separate r < Lγ, the interior, the trap-transition annulus and the exterior. Report the existing r < 5 Lγ statistic alongside prespecified alternative interior radii, not whichever radius gives the preferred sign. Coherent-vortex counts, size distributions and centroids would be useful if segmentation thresholds are fixed across cases and their sensitivity is reported. Keep the existing Eulerian summaries as independently defined quantities.

If the fine pair retains the qualitative paired effect and reasonable budget behavior, the next bounded temporal test should extend a matched polar/control pair from 75 to 125 days, with the existing separate extension plan and same-GPU uninterrupted/restarted validation. Reassess before activating 175 or 200 days. Do not select only the most visually attractive seed. The extension must retain original provenance, exact output/timestep barrier cadence, late-flow runtime estimates and daily native Float64 state for diagnostic audits.

If refinement changes the sign of the paired effect or makes the boundary/circulation diagnostics strongly resolution-dependent, investigate that sensitivity before committing to a long extension. If the H100 gate is too slow, design and validate scientifically equivalent segmented execution; do not loosen the wall-time or numerical criteria silently.

For an inverse-transfer claim, add an allocated-Unity spectral-transfer observer. Separate transfer from nonlinear relative-vorticity advection from redistribution associated with the imposed background PV, and account for filter/removal effects. Larger-looking structures or declining enstrophy alone do not establish an inverse cascade. A convection experiment should follow only after this dry two-dimensional reference and an independently resolved three-dimensional convective response are established.

No simulations were launched, modified or cancelled during this review.
