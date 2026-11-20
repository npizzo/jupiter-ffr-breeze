# Bounded physical feasibility pilots

On 2026-11-20 UTC the user requested preliminary evidence before committing to the five long compressible-convection runs. Those five production jobs are held: 2020719, 2020720, 2023193, 2023194, 2023196. Successful validation must **not** be interpreted as authorization to release these holds.

One Unity CPU job runs two pilots sequentially, each with a hard 1800-second process limit. Allocation maximum 62 minutes, 16 CPU cores, 24 GB memory, short QoS. The driver stops gracefully after 1200 seconds of integration wall time to leave room for compilation and checkpointing. No automatic extension or production submission. All computation and output remain on Unity.

- Käpylä A1-target and Chan–Mayr Case 4: same frozen model, box, vertical grid (128/80), thermodynamics, rotation, boundary conditions, heating, dissipation parameters and tiny random seed as the queued production configurations.
- Horizontal grid reduced to 64×64. Target at most t=500 in each paper's nondimensional time; a wall-limited shorter interval is explicitly incomplete evidence.
- Chan SGS viscosity is grid dependent: coarsening changes its effective viscosity. This is a feasibility pilot, not a quantitative reproduction or a vortex-sign prediction. The coarse grid may underresolve convection in either setup.
- Three-depth vorticity/tracer planes every 5 units; diagnostics every unit; restart checkpoints every 25 units and at graceful end. No full 3D dump campaign.
- End-of-pilot profiles: vertical RMS velocity, density, temperature, mass flux and enthalpy flux relative to horizontal mean temperature. These distinguish convective transport from velocity amplitudes alone.

## Review before any release

1. Check numerical stability, mass, boundary conditions and temperature/density positivity. Inspect the saved three-depth vorticity/tracer evolution for grid artifacts.
2. Determine whether actual convection and upward convective heat transport have developed. If the pilot has not reached this state, report it as inconclusive; do not call a quiescent startup a success or a lack of LSVs a negative result.
3. Inspect thermal evolution and total-energy change against integrated imposed flux, top flux and cooling, using matching intervals. Diagnostic endpoint quadrature and the model's nonconservative energy discretization must be distinguished. Do not claim closed energy/vorticity budgets from these checks.
4. Measure achieved timestep and time per simulated unit after startup. Extrapolating this CPU coarse-grid speed directly to full-grid A100 cost is invalid: retain the independent full-grid GPU timing checks. Cost estimates must also account for later turbulent CFL and the paper's long thermal relaxation.
5. Review pilot evidence and GPU validation before proposing one longer case. The current five production holds remain in force until that review is discussed with the user. A favourable early pilot cannot guarantee an LSV or the expected sign reversal.

Existing bounded GPU checks (2020718, 2023191) remain queued in normal QoS. The short-job service permits only one submitted job per user; use it for this CPU pilot bundle. No duplicate GPU jobs were submitted.
