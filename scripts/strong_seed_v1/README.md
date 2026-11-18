# Single no-shear, fixed-Rossby-number test

User-authorized revision, November 18, 2026. Exactly one case. The initially submitted amplitude-only Ro=0.3 job was held before allocation and revised in place; it never ran.

- Dry three-dimensional anelastic Breeze, f-plane at 55 N, no beta or initial/maintained shear.
- Relative to the recent smaller-box stable no-shear run: horizontal domain increases from 1500 to 6000 km in both directions; random-eddy wavelength increases from 400 to 1200 km. At unchanged initial Ro=0.1 and f, random-flow RMS increases from approximately 1.83 to 5.50 m/s.
- Retain 128 x 128 x 64 cells for speed. Horizontal spacing is 46.875 km: twice the original LSV run's spacing and four times the recent smaller-box spacing. Initial wavelength spans 25.6 cells; shortest seeded mode spans at least 20.48 cells. This is a coarse screening test, not resolution convergence.
- Same 120-km height, 162 W/m2 heating, cap stratification and relaxation, vertical sponge, heat controller, free-slip impermeable walls/top/bottom, viscosity and diffusivities. Same current random seed and odd-zonal-mode selection. No added drag.
- Compared with the original 6000-km LSV case, restore the domain, eddy scale and initial RMS, but use half the horizontal resolution and retain the newer sign-symmetric random-mode selection. It is therefore not an exact reproduction of its initial field.
- Fifty Earth days. Run one day first; require finite integration, intact cap, correct initialized RMS, frozen-source checks and measured speed consistent with the compute limit before continuing that checkpoint to day 50. Existing driver enforces CFL and cap checks throughout.
- One L40S GPU, four CPUs, 24 GB. Hard limit three GPU-hours including the gate; no automatic extension, retries or shear-case expansion.
- Unchanged output cadence and three-depth vorticity/tracer fields support synchronized movies. All computation and rendering stay on Unity.

The frozen existing driver is reused without edits. Its inherited campaign-description strings refer to the older pilot; `experiment.toml`, this README and the new config identify the present experiment. The segments join at day 1; remove the duplicate timestamp when rendering.
