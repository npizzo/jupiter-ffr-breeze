# Publication and next-experiment review

All scientific analysis, rendering and media validation run in allocated Unity jobs. Model sources, configurations, checkpoints and existing experiments remain unchanged. Finished assets and compact diagnostics are published to this repository and its GitHub Pages gallery.

| Report | Reporting job | Simulation prerequisite |
|---|---:|---|
| Gamma first matched 75-day pair | 1410628 | 1402730 |
| Ferrel four six-hour pilot cases | 1410629 | 1409341 |
| Ferrel three 48-hour cases | 1410630 | 1409345, with successful accounting checks for the complete chain |
| Gamma three-seed 1024² ensemble | 1410651 | 1402734 |

The independent Ferrel1024² and gamma2048² refinements remain separate validation stages. Their results will be included when complete; they do not prevent publication of explicitly labeled coarse milestones.

The gamma publisher independently checks native saved-field kinetic energy and both relative and total-PV enstrophy against the solver diagnostics. Its allocated first-case probe1410627 completed successfully. Gamma movies show two actual two-dimensional fields, without depth or tracer panels. Radial sign separation and centroids are descriptive diagnostics, not individual vortex tracking or proof of an inverse cascade.

Ferrel reports require complete outputs and successful simulation accounting. They retain exact sampled momentum/mass checks, direct mean/eddy tracer fluxes and explicitly approximate centered eddy-to-mean energy transfer. Each movie has three depth columns with vorticity above tracer, fixed scales and actual saved times. Startup behavior, symmetric controls, boundary sensitivity and refinement are distinguished. Encoded movies require full decode, frame-count, duration and size checks on Unity, followed by visual review and public playback verification.

Next simulations will be selected from the completed diagnostics and their numerical limitations. The decision record will state the evidence, alternatives and bounded compute plan before submission on Unity. No additional scientific conclusion follows merely from a successful reporting job or visual resemblance.
