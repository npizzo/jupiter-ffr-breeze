# Convection acting on initially seeded turbulence

The aim is to measure how imposed bottom heating modifies an existing random velocity field. Three cases share exactly the same initial velocity, tiny thermal perturbation, neutral lower reference layer, Jovian rotation and numerical settings. They differ only in bottom and upper thermal sources: heated/cooled (+5.4/−5.4 W/m²), cooled-only (0/−5.4), and unforced (0/0). This isolates the incremental effect of heating from freely decaying turbulence and rotational adjustment. Rotation is present in every case; this campaign does not isolate its effect against a nonrotating control.

## Initial field

Independent random Fourier phases and Gaussian amplitudes populate horizontal wavenumbers within 25% of 2π/(300 km), with three smooth vertical envelopes below 80 km. The 1,200-km domain and 128² horizontal grid resolve the seeded wavelengths with roughly 26–43 cells. Scalar toroidal and poloidal potentials generate a discretely divergence-free **mass flux**, then division by the anelastic reference density produces velocity. Vertical boundary flow vanishes exactly; the potentials taper smoothly before the stable cap. The field includes horizontal and vertical motion, rather than a purely two-dimensional barotropic pattern.

Define the reference eddy length as L = λ/(2π) = 47.75 km. The requested Rossby number is U/(|f|L) = 0.1, where U is lower-layer horizontal RMS speed measured on velocity faces. At f = 2.88096 × 10⁻⁴ s⁻¹, U is approximately 1.38 m/s. Lower-layer vertical RMS speed at cell centers is set to 0.2U, approximately 0.275 m/s. This convention is explicit: using the full wavelength as the length in the Rossby number would give a different value by 2π. The RMS vorticity/f diagnostic need not equal the normalization Rossby number.

The velocity seed is 31415, separate from the thermal seed 42. The thermal perturbation parameter remains 10⁻⁵ K. A random field with this prescribed spectrum is an initial condition, not already equilibrated turbulence. There is no continuous mechanical forcing. Early evolution includes adjustment; comparisons use identical initial fields so that adjustment is shared. Matching the field to Rossby 0.1 does not make the entire atmosphere geostrophically balanced.

## Runs and checks

Each production case uses 128 × 128 × 64 cells and lasts 48 Earth hours. Physical and mixing parameters remain those of the neutral-layer tests, including the 5.4 W/m² flux; beta and imposed shear remain zero. The initial eddies are well resolved relative to the mesh, but the small convection-driven structures identified in the neutral tests remain a separate horizontal-resolution concern.

Before the production cases, test exact seed reproducibility, the mass-flux divergence, boundary velocity, and RMS/Rossby normalization. A short GPU run then verifies that the actual solver's initial output retains the prescribed amplitudes, followed by the existing finite-output and source-budget checks for every case. The job stops if any check fails.

Compare kinetic-energy decay, horizontal/vertical energy profiles, buoyancy work, resolved heat transport, vorticity and horizontal spectra. Snapshot maps move from the weakly moving cap base to 20 km, inside the initially turbulent layer. Sections and snapshots at 6, 12, 24, and 48 hours accompany the comparisons. Small differences between chaotic realizations require care; a single matched seed and a 48-hour window cannot establish a statistically robust convective modification of turbulence.

Submit `sbatch slurm/turbulence48.sbatch`. One L40S, four CPUs, 8 GB host memory, and a one-hour limit cover validation and three sequential cases. Output is isolated under `data/turbulence48-v1`. The previous neutral tests and snapshots remain available.
