# Matched shear tests with stronger random flow at fixed Rossby number

On November 18, 2026 the user authorized adding two independent, concurrently eligible shear cases to the existing no-shear job 1983355. This explicitly expands the previous single-case scope to three total cases, not a wider sweep.

- Dry anelastic Breeze on the 55 N f-plane, beta zero; 6000 x 6000 x 120 km channel; 128 x 128 x 64 cells (46.875-km horizontal spacing).
- Identical random field across the three cases: seed 31415, current odd-zonal-mode selection, 1200-km wavelength, RMS approximately 5.50 m/s, initial random-eddy Ro=0.1 with length wavelength/(2 pi).
- The only physical difference from the current no-shear control is adding the initial, full-depth background U(y)=-S Ly/pi sin(pi y/Ly), with S/f=+0.05 (cyclonic) or -0.05 (anticyclonic). The pressure is initialized in discrete geostrophic balance. There is no maintained shear forcing.
- Keeping S/f fixed while widening the box increases the background peak speed from approximately 6.9 to 27.5 m/s. This is additional to the random RMS; Ro=0.1 describes the random eddies, not the total sheared flow. Total vorticity movies include the background. The existing diagnostics also retain vorticity with the instantaneous zonal mean subtracted.
- Same heating (162 W/m2), stable cap, sponge, heat controller, mechanical boundaries, viscosity, diffusivities, initial temperature perturbations and passive tracer as the no-shear control. No extra drag.
- Fifty Earth days per case; each runs a one-day stability and cost gate, then continues its checkpoint. Three GPU-hours maximum per case, nine across the three simulations. One L40S per job; scheduler permitting, all can run concurrently. No automatic extensions or retries.
- Compared with the original cyclone-forming run, the horizontal grid is coarser by two and the newer random-mode selection remains. A missing LSV at this resolution is not a converged negative result.
- A dependent CPU job joins days 0-1 and 1-50, checks timestamps and duplicate restart frames, renders the three individual synchronized 2x3 movies (vorticity above tracer at three heights) on shared fixed scales, and verifies full video decoding. It has a 30 CPU-minute limit. Publication to the existing GitHub gallery follows completed-media verification.

Existing model and no-shear source files are frozen. Each shear config is checked against the no-shear config and must differ only in name and shear_rossby. Inherited model metadata uses the older pilot's campaign label; the new config, experiment.toml and this document identify these runs.
