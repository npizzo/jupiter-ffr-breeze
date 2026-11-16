# Weaker forcing: one controlled change

Authorized November 16, 2026. This is a reduced horizontal 2D turbulence experiment, not a rapidly rotating convection simulation with temperature replaced by noise.

Three new one-CPU runs continue the identical mature σ_ref t=120 checkpoint to σ_ref t=1200 (reference σ_ref=0.6). Grid 128², periodic 2π square, f=1, beta=0, ν=0.001, cubic velocity drag coefficient 1, forcing band 6–8 excluding k_x=0, and random streams are unchanged. All science and rendering run on Unity. The existing unforced continuation is reused, saving another simulation; its timestep is 0.01 versus 0.005 for the new runs. The prior stochastic validation included a timestep-refinement check.

| Quantity | Original | New |
|---|---:|---:|
| Instability coefficient a | 0.009375 | 0.001875 |
| Stochastic expected energy input | 0.0318727811741 | 0.00637455623482 |
| Pulse interval | 0.2 | 0.2 |
| Stochastic velocity-kick amplitude relative to old | 1 | sqrt(0.2) |
| Viscosity | 0.001 | 0.001 |
| Cubic drag coefficient | 1 | 1 |

The instability term is a P_band(-∇²)ζ using physical Fourier k²; viscous damping uses the compatible native-grid Laplacian. Net linear growth is a k²−νλ_native. The fivefold coefficient reduction retains positive net growth, about tenfold slower than before; a tenfold coefficient reduction would be linearly damped on this grid. Actual energy input in the instability case depends on the evolving state and is not constrained to one fifth of its original long-time mean. Stochastic input is matched in expectation between sign laws; realized cross work is measured. Small-scale polarity is prescribed in the fixed-sign case, not derived from f-sensitive convection.

## Gate and resource bound

The check job tests scaling, net growth and pulse normalization, then runs the actual new driver for 40.2 time units in all three cases. It verifies matched initial fields and random locations, all 40 pulse phases, finite diagnostics, exact integrated energy accounting, and the existing native vorticity/third-moment/CFL guards. Production requires the PASS.json file and an afterok dependency. Three runs use one CPU each, at most two concurrently, four-hour Slurm ceilings and internal 3.5-hour stops. Partial results are retained; no automatic extension or parameter sweep. Existing simulations and raw outputs remain untouched.

## Diagnostics and movies

Diagnostic spacing is 201 steps with dt=0.005 and pulses every 40 steps. Successive diagnostic samples traverse all 40 discrete pulse phases; window means use complete 40-sample cycles. The previous 200-step spacing always sampled immediately after a pulse. Record full and low-pass vorticity skewness, native spectra and nonlinear flux, absolute large-scale energy, integrated energy work, RMS native/forcing tendencies, energy-input time, turnover time and pulse-to-flow velocity ratio. Budget rates use cumulative all-step increments and report the actual saved endpoints used for each window. Numerical closure alone is not physical attribution or resolution convergence.

The dependent Unity renderer combines new instability/random-sign/cyclonic-core cases with the reused unforced control in a synchronized 2×2 vorticity movie. All 901 actual fields use one fixed asinh scale, full extrema, no clipping, no interpolation. The movie samples post-pulse states; phase-balanced skewness statistics come from the separate 201-step diagnostics, not movie frames. There are no vertical layers or tracer in this explicitly 2D model.

## Connection to convection and limits

In the rapidly rotating reduced convection equations, depth averaging yields

∂t〈ζ〉 + J(〈ψ〉,〈ζ〉) = −〈J(ψ′,ζ′)〉 + ν ∇²〈ζ〉.

The forcing of the depth-independent flow is a correlation of depth-dependent motions; temperature and vertical velocity influence those motions. Replacing that term with random forcing is a statistical closure assumption, not the substitution θ→noise in Rayleigh–Bénard equations. It loses correlations, memory, and feedback of the large-scale flow on convection. A uniform-f incompressible horizontal 2D model does not itself distinguish the sign of f in its vorticity dynamics. A sign-biased stochastic source tests transfer of an imposed asymmetry; it cannot establish why real convection generates that asymmetry.

Van Kan et al. show that instability forcing and nonlinear damping can favor forcing-scale shielded vortices and suppress the inverse cascade. The current cheap pilot is inspired by that model, not a quantitative reproduction of its scale separation, dissipation or parameters. Weaker forcing does not guarantee an LSV: lower Reynolds number, cubic damping and limited scale separation remain possible restrictions. This bounded test isolates injection strength before considering any separate damping change.

Sources:
- [Maffei et al., inverse cascade in rapidly rotating convection, equation 2.9](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/on-the-inverse-cascade-and-flow-speed-scaling-behaviour-in-rapidly-rotating-rayleighbenard-convection/65EA7F62ED610264B8588035016CFBF3)
- [Guervilly, Hughes & Jones, large-scale vortices and cyclone/anticyclone asymmetry](https://arxiv.org/abs/1403.7442)
- [van Kan et al., spontaneous suppression of inverse energy cascade](https://arxiv.org/abs/2206.12355)
