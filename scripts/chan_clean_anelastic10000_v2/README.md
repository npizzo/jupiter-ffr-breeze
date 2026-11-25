# Anelastic continuation t6000 to10000

Authorized November25: continue completed job2249776 from its native t6000 checkpoint to10000. All physics, grid256x256x80, Omega3, beta0, nojet, WENO9/no closure, bottom flux, local cap cooling, environment and timestep safeguards remain unchanged. No reseeding. Time is nondimensional. Maximum24GPUhours; no automatic extension beyond10000.

A new driver version refreshes derived state after checkpoint load and performs the native endpoint comparison BEFORE any timestep, replacing the older exact-time callback check that could be skipped. CPU/GPU restart tests explicitly exercise this branch and verify the receipt. This also makes the diagnostic profile-drift reference the actual restarted state; the scientific equations and checkpoint are unchanged. The compressible energy proxy remains nonconserved for anelastic dynamics and is not a physical energy-closure diagnostic.

Depth-mean KE was still growing through6000. At6000 the thermal profile has drifted and minimum cap-mean N2 is negative; continuation is exploratory, not a thermal-equilibrium claim. Rendered three-depth movies through6000 are published separately. All scientific computation remains on Unity.
