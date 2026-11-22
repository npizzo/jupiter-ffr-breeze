# Ω=3 anelastic continuation

Continue the validated `omega3_anelastic_v4` full-resolution t=5 checkpoint to t=1000, with the same equations, grid, rotation, IC history, heating, cooling and transport. The gate verifies that `model.jl` is identical to the validated pilot apart from its experiment-scope description. The exact checkpoint is SHA-256 bound to the release record.

Diagnostics every 0.5, three-depth vorticity/tracer/temperature frames every 1, full fields every 20, checkpoint every 10. Model wall limit 7.5 hours inside an 8-hour allocation; the pilot suggests about 5 hours, with rate varying as the flow develops. Physical parameters are unchanged.

See the pilot README for initial-state mapping and limits: this compares anelastic and fully compressible equations while preserving entropy and velocity; it does not remove background stratification or isolate only acoustic waves. No automatic continuation beyond t=1000.
