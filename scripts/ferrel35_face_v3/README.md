# Day-300 native-face diagnostic restart

This isolated observer experiment restarts the completed 35°N heated beta-plane simulation at Earth day 300 and advances to day 301. It retains the original 256×256×64 grid, 6000×6000×120 km domain, bottom forcing, cap, thermal feedback, diffusion, initial state and time-step controls. It does not modify the wider-channel campaign.

Before production, `execute.jl` compares three hours of evolution against the frozen original `cascade35` driver. A failed physical-field comparison or budget check prevents production. The original checkpoint and CUDA environment are checked. All execution, analysis, and rendering occur on Unity.

## What is measured

At every SSPRK3 stage, the observer evaluates the pinned model's actual meridional and vertical momentum fluxes at their native faces. At those same faces it separately averages:

- the advecting mass transport;
- linearly interpolated zonal velocity;
- their instantaneous product;
- the actual WENO flux.

The product of the first two means defines mean transport. Product mean minus product of means defines eddy transport. WENO flux minus the linear product defines an explicit reconstruction remainder. All are differenced across the native momentum control volumes, then integrated with SSPRK3 weights. The periodic zonal flux divergence averages to zero. Their sum is checked against the independent native advection operator. This decomposition states its interpolation convention; the reconstruction remainder is not automatically numerical dissipation or a physical Reynolds stress.

The potential-temperature budget independently resolves native advection, explicit diffusion, stable-cap relaxation, thermal feedback, and bottom heat flux. The callback precedes application of flux boundary tendencies, so the prescribed bottom rho-theta flux is included explicitly. Both the RHS partition and the actual prognostic rho-theta increments are checked. Heat transport receives the same native-face decomposition. The column enthalpy diagnostic uses cp times reference Exner times rho-theta; it is not a complete total-energy conservation test.

Every step is audited. Three-hour frames retain synchronized three-depth vorticity and tracer output. A one-day movie uses one frame per second to show all nine available states, including the original day-300 frame. Analysis uses a fixed interior region |y|≤2400 km and height 10–75 km for future width comparisons.

## Interpretation limits

This is a diagnostic validation, not a one-day statistical test for Ferrel cells. Establishing eddy maintenance still requires sustained averages, the wider-channel comparison, and eventually cap and resolution sensitivity. The mean and eddy decomposition must not be replaced by separate nonlinear WENO evaluations of mean and fluctuating states.

## Execution

From the Unity project root:

```
sbatch slurm/ferrel35-face-v3.sbatch data/ferrel35-face-v3
sbatch --dependency=afterok:SIMULATION_JOB slurm/ferrel35-face-review-v2.sbatch data/ferrel35-face-v3 rendered/ferrel35-face-v3
```

Output directories must be new. Do not rerun these commands if a job is already active. The development checks corrected mass-momentum indexing and a snapshot-time comparison against a stale thermal tendency. Native-stage audits remain strict; snapshot writing does not recalculate a thermal RHS against a different feedback value. Production remains gated behind the complete mature-state comparison and numerical audits.
