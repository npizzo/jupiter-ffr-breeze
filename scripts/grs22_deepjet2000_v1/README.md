# One GRS-oriented compressible beta-channel experiment

User authorization: run a compressible case to t=2000 closer to GRS conditions.

This is a fresh, reduced, dry-convection experiment at 22 degrees south, not a dimensional Jupiter reproduction and not a continuation with a changed rotation rate. Omega=3 stays fixed. A deep opposing jet replaces the previous weak cap-confined target. The target is eastward south of the center and westward north of it: dU/dy<0, positive relative vorticity, anticyclonic for negative f. The same mixed-sign eddies are rebalanced at the new latitude; they do not predetermine a dominant vortex sign.

## Explicit configuration

- Compressible Breeze; WENO9; no closure, added drag, or sponge.
- 128 x 256 x 80 in 4H x 8H x H; periodic east-west, impermeable stress-free north/south and top/bottom.
- Latitude -22 degrees; Omega=3; full tilted rotation; beta=6*cos(22 degrees)*250/71492. The 250-km length mapping sets beta/f only; it does not calibrate the temperature, sound speed, gravity, heat flux or all dimensional units.
- Retain R=1, cp=2.5, g=12, initial stratification/cap, bottom heat flux0.03125, conduction and horizontal-mean cap thermal relaxation.
- Initial mixed-sign eddy RMS0.136, same phases and scales; jet peak0.125, widthH, depth-independent, smoothly tapered near walls. Full-column zonal-mean relaxation time200. Jet mean flow is maintained, not every instantaneous eddy velocity.
- Fresh t0 to2000; output interval2, three vorticity/tracer heights; full fields every100; rotating checkpoint every20.
- Stop for nonfinite/nonpositive fields, Mach>=0.8, mass error>1e-7, temperature-profile departure>50%, or theta departure>factor4. Record energy residual without artificial heat compensation. These are exploration guards, not accuracy guarantees.
- One case only. Existing high-latitude runs provide context, not a matched control isolating latitude versus shear.

## Why this is closer, and where it is not

The GRS is near22S, and its surrounding jets penetrate deeper than the observed vortex. A full-depth target is a better first representation of that relationship than a weak target confined above z/H0.82. At one jet width, imposed speed/(beta*width^2) is about6.4; central shear/(beta*width) is about10.6. These are scaling estimates, not evidence of trapping. A uniform current cannot prevent beta drift; the shear, PV gradient, and vortex interaction matter. Vortex position and persistence must be measured relative to the moving zonal flow.

The domain aspect ratio and thermodynamics are still the prior idealized convection model. In particular, mapping both Omega and the 250-km length to Jupiter would not also match Jupiter's gravity or sound speed. This run cannot be called a quantitative GRS analogue. There is no inserted mature GRS, moisture, or PV staircase. Positive zeta is anticyclonic in this southern-hemisphere run; plot zeta/f so negative remains anticyclonic.

## Launch and review

CPU preflight parses all Julia source, runs the actual driver and its output/diagnostic paths, verifies restart equivalence, verifies heat-source accounting and resting balance, and checks jet orientation and balance. An exact-source/environment receipt permits one GPU submission. On the GPU, repeat restart checks then integrate the production grid to t10. Require positive mean cap N2, Mach<0.8 and mass error<1e-7; only then restart the unchanged state to2000 in the same allocation. No automatic retuning or parameter changes.

All computation and rendering stay on Unity. The inherited unresolved compressible energy residual remains a limitation. Future interpretation must distinguish vortex formation from initial large-eddy survival, and zonal-flow advection from meridional confinement. The t1000 experiments are preserved unchanged.

References: Shetty, Asay-Davis & Marcus (2007), https://arxiv.org/abs/0707.4193 ; Parisi et al. (2021), https://doi.org/10.1126/science.abf1396 ; Cai, Chan & Chow, https://arxiv.org/abs/2111.11003 .
