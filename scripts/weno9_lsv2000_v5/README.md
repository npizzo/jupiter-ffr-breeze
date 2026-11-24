# Unchanged four-case continuation to t2000

User authorization: “launch the longer runs as well”. Continue all four completed cases from their own t1000 checkpoints: compressible nojet/jet, anelastic nojet/jet. Fixed Omega3, latitude67.5N, original beta, original heating/cap and full-column mean-flow relaxation in jet cases, WENO9/no closure, grid128x256x80, identical guards and outputs. No reset, parameter switch, ramp, or new seeds. Preserve all t0-1000 output.

These are persistence extensions. The compressible cap has weakened relative to the anelastic cap and the compressible energy residual remains unresolved. The anelastic compressible-energy proxy is not an invariant. Longer integration does not fix these limitations. Positive zeta/f remains cyclonic. The separate grs22_deepjet2000_v1 is a fresh experiment and is not part of this matched continuation.

Exact-driver CPU preflight (including all enabled diagnostics/outputs and restart) precedes an idempotent GPU-array submission. Each GPU task re-verifies receipt and restart paths, then restores its own t1000 state and original budget and health references. No further automatic extension after2000. Scientific computation stays on Unity.
