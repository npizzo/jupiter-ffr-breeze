# Matched horizontal-resolution test at day 300

Restart the 35°N heated beta-plane channel at Earth day 300 and compare one day on 256×256×64 and 512×512×64 grids. Horizontal spacing is 23.4375 versus 11.71875 km. Both retain the 6000×6000×120 km domain, original 64 vertical levels, bottom forcing, stable cap, thermal feedback, explicit viscosity/diffusivity, and native-face diagnostics. Both use a 15-second maximum timestep. The already completed coarse 30-second result supplies a timestep-sensitivity control.

This targets the predominantly horizontal numerical reconstruction correction found in the day-300–301 audit. It is a short resolution-sensitivity test, not proof of convergence of a statistical Ferrel circulation; it does not test vertical resolution or channel width.

## Restart transfer

Only Unity reads and remaps model fields. The original day-300 checkpoint is exported to a frozen bundle of native prognostic mass densities, pressure, clock, and thermal-controller state. No fine-scale random perturbations are added.

Momentum uses a conservative staggered prolongation satisfying D_f P_F = P_C D_c. Center quantities are prolonged with linear conservative slopes; normal-face midpoint values receive the matching second-difference correction. Applying the same transverse center operator in all components preserves discrete anelastic mass continuity without a new pressure projection. Wall-normal boundary fluxes remain zero. Restricting the fine fluxes to coarse faces recovers the original fields. Scalar mass densities use conservative minmod-limited slopes to preserve extrema and total heat/tracer content. The vertical grid and reference density remain identical.

Validation checks the commuting-divergence identity on manufactured arrays; restriction/conservation and scalar bounds on the actual checkpoint; a same-grid bundle restart against a direct checkpoint restart; native-stage momentum and potential-temperature closure; and a full fine-grid one-hour pilot. Production requires a successful gate, source/bundle hashes, a bounded runtime estimate, and an initial kinetic-energy change below 10%, with the actual change reported. Passing this generous remap guard is not evidence of dynamical convergence.

## Comparison

Compare both grids on their native grid and after averaging fine diagnostic tendencies over the same coarse meridional cells. Use the fixed interior |y|≤2400 km and height 10–75 km. Report absolute mean, eddy and reconstruction accelerations, not only their ratio; compare the 15- and 30-second coarse cases. Separate the first six simulated hours from the final eighteen to expose remapping adjustment. Preserve native heat-budget closure and report enthalpy storage and imposed-source contributions.

Movies use actual saved fields at the same three vertical levels and times. Stack coarse vorticity/tracer above fine vorticity/tracer, with three synchronized columns, shared fixed vorticity ranges per depth, and the existing fixed tracer contrast. All simulation, analysis, and rendering stay on Unity. The wider-channel campaign is unchanged.
