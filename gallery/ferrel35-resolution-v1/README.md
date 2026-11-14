# 35°N horizontal-resolution audit, days 300–301

This directory holds the synchronized coarse-over-fine three-depth movie and native momentum/thermal budget comparison for the authorized Unity resolution test. Both cases use a 15-second timestep ceiling, the same mature checkpoint, physical domain, forcing and 64 vertical levels. The existing coarse 30-second restart is a timestep control.

The 512² state uses conservative divergence-preserving refinement of the 256² native mass fluxes; its initial kinetic energy is 1.7075% lower. This is a one-day sensitivity test, not a statistically equilibrated resolution-convergence claim.

The movie, numerical audit and publication validation are added after the queued Unity review completes. Source and methods: [resolution experiment](../../scripts/ferrel35_resolution_v1/README.md).
