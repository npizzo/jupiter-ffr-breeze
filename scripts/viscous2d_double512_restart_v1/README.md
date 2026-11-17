# Doubled-forcing restart and t=0–160 movies

The first doubled-forcing run reached its CFL guard near t=161 (CFL0.50007249; limit0.5). Preserve its full output. Publish only the frozen 321-frame prefix through t=160; label this partial result clearly. Reuse the unforced and original-forcing controls in a synchronized three-column vorticity movie. No depths or tracer exist in this horizontal 2D model. Shared fixed asinh colors of width0.1 cover all extrema across the three cases and published times.

Restart from the exact saved t=160 velocity checkpoint with dt=0.00125, half the previous0.0025. Retain512², periodic2π×2π, f=1,beta=0, viscosity0.0001, no drag, coefficient0.0004, forcing band12–16 excluding kx=0, RK3 and centered advection. The same CFL<0.5 guard remains. Short validation verifies checkpoint/frame agreement, finite fields, serialized restart, a local timestep comparison, efficient execution, and evolution past the old stopping point to t=162. Production is gated on PASS and split into t=160–240,240–320,320–400. Checkpoint every10 model units; movies every0.5. No broader convergence or physical-attribution claim is made from the short check.

All computation and rendering stay on Unity. Failed checks prevent production. No automatic timestep/forcing change. Preserve existing simulations and published media.

