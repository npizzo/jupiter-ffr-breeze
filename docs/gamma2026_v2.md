# SI-informed polar gamma-plane reproduction pilot

This isolated v2 campaign follows the corrected main paper and recovered
supplementary information, section 1 equations 1–3, of
[Siegelman, Young and Ingersoll (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170063/),
DOI 10.1073/pnas.2120486119. It is an independent reproduction attempt, not a
bitwise historical rerun. No explicit convection or sustained forcing is present.

The official GeophysicalFlows SingleLayerQG solver is pinned to commit
`51ca98a0d471ce4ef8bbba8b73daf2bad54d53fb`, in the separate gamma2026 Julia
environment and depot. Existing Breeze and v1 source files are unchanged.
The evolution is Dq/Dt=0, with q=zeta+eta and infinite deformation radius.
Inside a circular trap, polar eta=-gamma*r^2/2; the flat control has constant
eta=-gamma*r_trap^2/2. Both vanish outside the trap on the doubly periodic square.

## Verified reference setup and declared departures

- gamma=7.8690e-20 m^-2 s^-1; domain RMS velocity U=80m/s;
  Lgamma=(U/gamma)^(1/3); square width12Lgamma; trap radius5Lgamma.
- Linit=200km is an inverse wavenumber, with wavelength about1257km. The
  November2022 correction removing the spurious2pi convention is respected.
- Independent complex-normal Fourier amplitudes with variance1 occupy
  0.95/Linit <= |k| <= 1.05/Linit. Hermitian symmetry gives a real field.
- SI quartic taper (r_trap-r)^4 inside the trap, zero outside; subtract the
  periodic inversion's zero mode, then normalize domain KE to U^2/2.
- The integer-mode random draw order is independent of grid size. Each polar/
  flat pair shares its initial realization, as does the1024/2048 seed11 pair.
- Float64, FilteredRK4, no explicit viscosity or drag; strict two-thirds
  dealiasing. Historical filter coefficients, software version and seeds are
  unavailable, so the current pinned implementation is documented instead.
- The paper uses fixed62.5s. This pilot uses adaptive advective CFL0.25 with
  a62.5s ceiling for stability. This is an explicit numerical departure.

## Finite first batch

Six75-Earth-day reference-grid cases at1024²: polar and flat traps for seeds
11,22,33. A matched seed11 pair at2048² follows a separate2048² GPU validation
and throughput gate. SI reports1024² as default and2048² early tests as nearly
identical. Eight-year integrations are not silently included in this pilot.

The first75days test early mixing and cyclone/anticyclone segregation. The SI
states mature crystals take roughly2–3Earthyears; absence at75days cannot
reject the mechanism. Longer runs and any4096² stage follow a review of
throughput, energy loss, physical diagnostics and the paired grid comparison.

## Validation, accounting and resource bounds

CPU and GPU tests check an analytic stationary Fourier mode, invariant budgets,
energy normalization, identical modal realizations across grids, matched
initialization within roundoff, finite polar-trap evolution, restart equivalence
and CPU/GPU agreement. Each grid's GPU gate benchmarks100 steps after warmup.
A run requires matching core-source hash and grid; estimated75day cost must
be below80% of its4.5hour soft wall budget. Each production job requests one
L40S,4CPUs,32GB and6hours. Gates are numerical sanity checks, not proof of
statistical convergence or reproduction.

Runs stop if energy exceeds102% or falls below70% of initial, becomes nonfinite,
or exceeds4.5wallhours; a wall-limit stop leaves a resumable checkpoint.
The energy bounds detect gross problems, not acceptable paper agreement.
The SI reference loses6.05% KE and92.18% relative-vorticity enstrophy by75days;
those are comparison targets across realizations, not exact-seed assertions.
Initial and daily full-grid energy, total-PV enstrophy and relative-vorticity
enstrophy are recorded separately. Do not compare total-PV enstrophy to the
SI relative-vorticity loss table.

Daily native1024² or2048² vorticity/PV fields and latest full spectral restart
stay on Unity. No stride downsampling or local field download is required.
Expected output is approximately9GB and below12GB for the complete batch.
Rendering and all diagnostics remain on Unity. Only finished publication assets
may be transferred for the gallery upload if needed.

The sequential dependency chain stops after a failed gate or run. It neither
changes nor cancels any old job. Submit once to a new output root:

```sh
bash scripts/gamma2026_v2/submit.sh data/gamma2026-v2a [SETUP_JOB_ID]
```

When the setup job has aged out of the Slurm controller, verify that the shared
environment is installed, then omit the optional SETUP_JOB_ID. Otherwise the
CPU unit gate has an afterok dependency on that setup job. All subsequent gates
and runs always retain their successful-predecessor dependencies.

Checkpoints permit restart only with the same core source and configuration.
The exclusive submission record prevents accidental duplicate chains.

## Interpretation and presentation

Required assessment includes energy/enstrophy loss, radial PV profiles,
cyclone/anticyclone segregation, crystal radius, vortex counts/spacing, merger
rates and seed/grid sensitivity. Morphology alone does not demonstrate an
inverse energy transfer. The initial phase is a baseline for subsequent
convection experiments, not a model of heat-driven convection itself.

This is a genuinely2D barotropic model with no depth coordinate. Its movies
must be labeled2D and cannot honestly supply three physical depth panels.
The parallel Breeze experiments retain real three-depth synchronized
vorticity/tracer movies. Neither output represents cloud radiance directly.
