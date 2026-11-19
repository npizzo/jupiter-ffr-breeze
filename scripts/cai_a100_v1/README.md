# A100 relaunch of the three Cai benchmark cases

Authorized November 19 after the prior batch stopped at a cost forecast. This supersedes the earlier timing-job dependency, not the physical experiment.

## What happened and what changes

GPU validation 2005414 passed conduction, initialization, Coriolis and serialized-restart checks; the production-grid smoke integration completed through t=2. It then failed a four-hour cost threshold using **L40S** timing. Production had requested A100/H100, so that was not a hardware-matched runtime measurement. All three dependent production jobs and the renderer were cancelled with zero runtime. The validation used 211 GPU-seconds.

The relaunch uses one **A100** per independent production case. There is no separate timing allocation and no early runtime-forecast stop. Timing is reported inside each production allocation after t>=100 using several output intervals. The hard four-hour allocation and the model's 13,800-second stop remain; unfinished runs retain checkpoints and are explicitly partial. This does not promise t=5000 within four hours.

The passed implementation tests are reused through a receipt linked to their actual log and unchanged model SHA. The prior cost failure remains recorded as a failure; it is not relabeled as an all-tests pass. Source/configuration hashes and GPU type are checked before execution.

## Unchanged assumptions

The exact frozen `scripts/cai_benchmark_v1/model.jl` and three original configurations are used: dry compressible ideal gas, 4H×4H×H periodic box, 128×128×64, constant bottom flux and fixed top temperature, local Smagorinsky viscosity, zero initial velocity with seed 31415 thermal perturbations, and Omega=0 / 0.64 / 2.56 at 90 degrees. No cap, shear, beta, moisture, numerical-method change or extra case. Target t=5000 in benchmark units. See the [benchmark specification](../cai_benchmark_v1/README.md) for all parameter values and differences from the paper.

Maximum production allocation: 3×4 = **12 GPU-hours**. No automatic extensions or retry. Movies and diagnostic review use the existing renderer on Unity after production, preserving synchronized vorticity/tracer at three heights. A wall-limited result is labeled partial, not an absence of LSV formation.

## Rules for subsequent work

- Separate equation/restart correctness from cost estimates and from physical evidence.
- Reuse valid tests when model/configuration hashes are unchanged; test changed code where relevant.
- Obtain runtime estimates from the actual grid, output cadence, GPU class, and settled stepping. State their uncertainty.
- Prefer one checkpointed production allocation to repeated validation/continuation queues. Retain hard resource limits and normal stability checks.
- Do not call jobs running until the scheduler confirms it; distinguish queued, running, stopped and cancelled.
- Before another campaign, disclose every physical/numerical change and its resource bound. Never interpret an under-resolved or unfinished case as a physical null result.
