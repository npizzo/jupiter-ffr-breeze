# 55°N unforced beta-plane control, 100 Earth days

This experiment starts from the same initial random velocity, potential temperature, reference density, and passive-tracer fields as `cascadebeta100_heated`. It is a fresh day-0 run, not a restart from the heated day-100 state. Grid: 256 × 256 × 64; domain: 6000 × 6000 × 120 km; seeded eddy wavelength: 1200 km; Rossby number: 0.1. Rotation and beta are identical to the heated 55°N case.

**Fully unforced** means zero bottom heat flux, no cap-temperature relaxation, no global thermal-budget correction, no vertical-velocity sponge, and no imposed momentum forcing. Molecular/eddy viscosity and scalar diffusivity remain as in the heated run. Gravity, Coriolis acceleration, and buoyancy remain active. The initial stable cap above 80 km remains, but is not restored. Walls remain impermeable and free-slip; x is periodic.

This is a freely evolving/decaying control. Its comparison with the heated run changes the full external forcing package, not just bottom heating; it cannot alone isolate the effect of bottom convection from cap relaxation and sponge effects.

A GPU preflight checks that all unchanged configuration fields and initialization code match the heated case, compares actual saved initial velocity, temperature, vorticity and tracer fields, verifies uninterrupted-versus-checkpoint-restarted evolution, and runs one day at production resolution. It checks finite fields, tracer bounds and conservation, potential-temperature integral conservation, zero imposed sources, cap stability, wall conditions and checkpoint output. Production is ten 10-day segments, gated by successful validation and each preceding segment. Source/config/environment hashes are frozen by the preflight. All computation and rendering run on Unity; no CPU fallback.

The dependent renderer produces a synchronized 4 × 3 comparison: unforced vorticity and tracer above heated vorticity and tracer, at three depths, over days 0–100. Fixed physical asinh color scales are shared between cases. Movies and metadata can be published after validation; no raw fields belong in GitHub.

Submit from the Unity repository root:

```bash
bash scripts/beta55_unforced100_v2/submit.sh data/beta55-unforced100-v2 data/beta55-unforced100-preflight-v2
```

## Relaunch correction

The first validation (1526400) stopped before model initialization because the shared CUDA launcher did not list the new campaign entry points. Version 2 uses a dedicated launcher that preserves the verified CUDA 12.2 environment and GPU checks, and permits only this campaign’s two drivers. The shared launcher and active 35°N simulations are unchanged. The renderer now includes its missing collection helper and permits the intended cooling-model difference between the heated and unforced cases. No physical parameters or initial conditions changed.
