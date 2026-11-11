# Polar gamma-plane and flat-trap pilots: completed pairs

These are **two-dimensional, unforced barotropic QG** simulations, with no convection, vertical levels or passive tracer. They independently implement the reference supplementary-information initialization: a ±5% spectral annulus, quartic radial taper and common initial RMS speed. This is a 75-Earth-day early-organization test, not a mature vortex-crystal reproduction. The current solver/filter, adaptive timestep and new seeds differ from the historical calculation.

## Movies

- [gamma_flat_n1024_seed11_0_75days.mp4](gamma_flat_n1024_seed11_0_75days.mp4): synchronized polar/control comparison, all 76 saved times, two actual 2D fields.
- [gamma_flat_n1024_seed22_0_75days.mp4](gamma_flat_n1024_seed22_0_75days.mp4): synchronized polar/control comparison, all 76 saved times, two actual 2D fields.
- [gamma_flat_n1024_seed33_0_75days.mp4](gamma_flat_n1024_seed33_0_75days.mp4): synchronized polar/control comparison, all 76 saved times, two actual 2D fields.

Movies use a symmetric asinh color map with physical vorticity labels, fixed across all times and selected cases. The range includes all selected saved extrema. No field clipping, duplicated depths or temporal interpolation is used. Each MP4 was fully decoded on Unity and its frame count and duration verified.

**Total-PV diagnostic audit has unresolved discrepancies.** The original relative tolerance of 5×10⁻⁶ was retained and failed at the times listed below. This publication does not certify all budget checks or make a tight total-PV conservation claim. Independently reconstructed kinetic energy and relative-vorticity enstrophy pass that same threshold at every saved time, supporting publication of the actual vorticity fields and descriptive vorticity statistics. The total-PV issue is not assumed to be Float32 roundoff.

| Case | Maximum relative PV mismatch | Failed saved times / total | Failed Earth days |
|---|---:|---:|---|
| polar_n1024_seed11 | 2.6839888e-06 | 0 / 76 | none |
| flat_n1024_seed11 | 5.8960232e-06 | 1 / 76 | 1 |
| polar_n1024_seed22 | 2.4364616e-06 | 0 / 76 | none |
| flat_n1024_seed22 | 2.7894501e-06 | 0 / 76 | none |
| polar_n1024_seed33 | 1.2707765e-06 | 0 / 76 | none |
| flat_n1024_seed33 | 6.6257974e-06 | 1 / 76 | 1 |

## Numerical and descriptive diagnostics

| Case | Grid | Days | E/E₀ | Solver PV enstrophy / initial | Relative enstrophy / initial | Mean radius (negative − positive) / Lγ |
|---|---:|---:|---:|---:|---:|---:|
| polar_n1024_seed11 | 1024² | 75 | 0.936903 | 0.071875 | 0.064930 | -0.62051 |
| flat_n1024_seed11 | 1024² | 75 | 0.937343 | 0.097034 | 0.069219 | -1.13856 |
| polar_n1024_seed22 | 1024² | 75 | 0.937154 | 0.075089 | 0.067055 | -0.34293 |
| flat_n1024_seed22 | 1024² | 75 | 0.937259 | 0.097582 | 0.069566 | -0.95469 |
| polar_n1024_seed33 | 1024² | 75 | 0.936431 | 0.067285 | 0.059740 | -0.51045 |
| flat_n1024_seed33 | 1024² | 75 | 0.936984 | 0.096449 | 0.068752 | -0.96554 |

| Grid / seed | Polar change from day 0 | Flat-control change | Polar minus flat change |
|---|---:|---:|---:|
| 1024² / 11 | -0.62043 | -1.13847 | 0.51804 |
| 1024² / 22 | -0.34294 | -0.95470 | 0.61176 |
| 1024² / 33 | -0.51050 | -0.96559 | 0.45509 |

All separation changes are in Lγ. Subtracting day 0 removes the finite-realization initial sign imbalance; subtracting the paired flat-control change helps isolate the effect of the polar PV geometry. It is a descriptive matched-pair effect, not a significance test.


Every native Float32 snapshot was independently compared against the solver's Float64 kinetic-energy, relative-enstrophy and total-PV-enstrophy diagnostics (unchanged relative tolerance 5×10⁻⁶, absolute tolerance 10⁻¹⁸). The per-field results above distinguish passed checks from failed PV checks. Kinetic energy was reconstructed by periodic spectral inversion of saved relative vorticity at infinite deformation radius. Relative enstrophy is not an invariant when the background PV varies in space; its decline is not, by itself, numerical dissipation. The recorded total-PV diagnostic is ½⟨(ζ+η)²⟩. The PV panel shows both solver and saved-field reconstructions against the same solver initial denominator; CSVs preserve both series and their discrepancy at each time. Finite-grid filtering can remove energy and PV enstrophy, but this observation does not explain the solver/snapshot mismatch.

The radial-separation statistic is the mean radius of negative vorticity minus that of positive vorticity, each weighted by its sign's vorticity magnitude **inside r < 5 Lγ**. Positive values indicate that positive vorticity is more centrally weighted. CSVs also contain sign-weighted Cartesian centroids, RMS radii, core fractions and the initial separation. These quantities describe an Eulerian field; they do not identify or track individual vortices, and their values should be compared against the matching flat-trap control and seed variation. No inverse-cascade claim is made: this publication does not calculate nonlinear spectral transfer.

Resolution differences between the 1024² and 2048² matched-seed runs are sensitivity evidence, not a convergence proof. Day-75 morphology alone cannot establish a mature polar crystal. Comparisons with convection require a separate, dynamically matched 3D experiment.

## Campaign coverage at rendering

- polar_n2048_seed11: NOT_FOUND (not included in this publication).
- flat_n2048_seed11: NOT_FOUND (not included in this publication).

## Provenance

All analysis, rendering and media validation ran in an allocated Unity CPU job. Raw model fields remain on Unity. Per-case metadata, native diagnostic CSVs, frame hashes and video normalization records accompany this report. The JSON manifest records the exact renderer and source revisions.
