# 512² unforced resolution check

Authorized November 16, 2026, after the 256² unforced test passed the inverse-transfer checks but exceeded the preset high-wavenumber enstrophy limit (1.85898% versus 1%). The three short forcing probes had also completed, but the dependent long runs and movie were cancelled before starting.

Assumptions disclosed before this submission:

- Only the unforced case; no production restart or forcing sweep.
- Same periodic 2π square, f=1, beta=0, ordinary Laplacian ν=10⁻⁴, and exactly zero drag and forcing.
- Grid 512²; the timestep remains 0.0025. Model time 0–8, matching the existing coarse check.
- Preserve the actual coarse initial corner-vorticity Fourier coefficients at k=12–16 by zero-padding. Reconstruct the fine-grid velocity using the compatible native curl inversion. No new random realization or amplitude renormalization. Record the small discrete kinetic-energy/RMS-velocity change.
- Retain all existing inverse-transfer, budget, CFL and high-wavenumber thresholds. A pass alone does not authorize long production or establish convergence.

The initializer checks restriction of the refined vorticity to coarse corner coordinates, the recovered curl, divergence, zero means, spectral support, and that the discrete initial-energy difference is below 1%. Initialization failure prevents integration.

The analysis applies the same criteria: negative mean native Π(K) at K=4,6,8 for t≥4; at least 75% of those samples negative; low-k energy gain at least 1% of initial total energy; under 1% enstrophy at k≥N/3; zero forcing/drag work; integrated relative energy closure below 10⁻⁸; CFL below 0.5. The first and last spectral sample times are recorded because they do not coincide exactly with t=0 and t=8. This matches the coarse check's existing definition.

Because refining the grid moves N/3 to a higher physical wavenumber, the report also compares enstrophy above the common fixed cutoff k=256/3, the low-wavenumber fluxes and absolute large-scale energy. Review these quantities before inferring convergence or choosing production resolution. No cutoff or threshold is adjusted after inspecting the fine result.

The job renders a synchronized 256²/512² vorticity comparison on Unity: 17 saved states at t=0,0.5,…,8, two frames/s, with a shared fixed asinh scale spanning full extrema. This is one actual horizontal 2D plane per grid, not a three-depth simulation. Diagnostics, finished movie, provenance and a PASS/FAIL report remain on Unity for review. Analysis writes the report and movie even when the scientific screen fails, then exits nonzero.

Resource limit: one CPU, 16 GB, two-hour job ceiling and 90-minute internal integration stop. No automatic extensions, additional simulations or production release. No raw model data is downloaded.

## Metadata correction during startup

The frozen initializer has two incorrect printed RMS expressions: Julia parses `2E0` and `2E1` as exponent literals. Its velocity fields, energy calculation, Fourier-mode checks and simulation are unaffected. The separate Unity metadata job runs `repair_metadata.py`, which preserves `initialization.original.toml`, corrects the two RMS entries from the recorded energies, and records `metadata_correction.json`. The simulation driver independently records the actual RMS velocity in its configuration. No simulation was repeated for this reporting correction.
