# The same polynomial equation: potential and finite-mode projection

These are mathematical reductions of the imposed 2D equation, not a derivation from convection. The forcing band, coefficients and native-grid symbols remain fixed.

## A precise phase-ordering limit

With advection suppressed, define P₀g=g−⟨g⟩ and the discrete self-adjoint operator A=−νΔ−L. The equation is ζ̇=−P₀g, where

g=Aζ−Bζ²+Cζ³=δF/δζ,

F=⟨ζAζ⟩/2−B⟨ζ³⟩/3+C⟨ζ⁴⟩/4.

Consequently dF/dt=−⟨(P₀g)²⟩. This is a semi-discrete identity; a finite timestep may introduce an error. Production records F and the gradient norm, while the release checks independently compare the potential gradient to the equation's actual source. F is not kinetic energy. Restoring self-advection generally removes monotonicity.

One analytically soluble limiting comparison omits L, considers well-separated interfaces, and retains the zero-mean constraint. It is not a change to the experiment. Set m=B/(3C), ζ=m+φ. Up to a constant and a term fixed by the mean constraint, the local potential becomes

V(φ)=Cφ⁴/4−B²φ²/(6C).

Its common-tangent coexistence values are ζ±=m(1±√3). For B>0 and ⟨ζ⟩=0, the positive phase occupies fraction (√3−1)/(2√3), approximately0.211; the ideal two-level field has skewness √2. Finite interfaces, finite domain, the actual band forcing, and advection invalidate treating these as predictions of the full experiment. They are useful falsifiable benchmarks: substantially different core fractions or skewness require more than this simplest equilibrium picture. The full band-forced potential can instead support finite-wavelength or localized states.

## Exact finite-mode equation for the native discretization

Represent ζ(x)=Σₖ zₖ exp(ik·x), with z₋ₖ=zₖ* and z₀=0. Use the original spacing d=2π/512, even when retaining few modes. Define

Dⱼ⁺(k)=(exp(ikⱼd)−1)/d, Dⱼ⁻(k)=(1−exp(−ikⱼd))/d,

λ(k)=4[sin²(kₓd/2)+sin²(kᵧd/2)]/d².

Native velocity symbols are Uₖ=Dᵧ⁺zₖ/λ and Vₖ=−Dₓ⁺zₖ/λ. Apply the same face/center averages and conservative momentum-flux differences as the original centered Oceananigans operator. Their curl gives a quadratic interaction tensor A(k;p,q), evaluated explicitly in `galerkin.py`. The reduced equation for k in the retained set S is

żₖ=[a|k|²M(k)−νλ(k)]zₖ
    + Σ_{p,q∈S; p+q=k}[A(k;p,q)+B]zₚz_q
    − C Σ_{p,q,r∈S; p+q+r=k}zₚz_qz_r.

M is the unchanged band12–16 mask excluding kx=0. There is no constant-f term in the solenoidal horizontal vorticity equation. The phenomenological B=bf remains imposed. Setting A=0 is precisely the advection-off intervention.

The cubic term must keep all intermediate p+q values, not project ζ² back onto S before multiplying by ζ. Premature projection would change saturation. The implementation builds that complete pair convolution, then contracts it with the retained z. Retaining conjugate pairs preserves real fields; excluding k=0 enforces circulation. With |k|≤64, all cubic sums fit without wraparound on the original512 grid.

An independent full-grid evaluation of each nonlinear term verifies the coupling coefficients. The advection tensor is checked for zero kinetic-energy work, using E=Σ|zₖ|²/(2λₖ). Selection and forecast errors are separate from these algebraic checks: exact projection of the equations does not guarantee an accurate low-dimensional trajectory.

## Physical interpretation remains conditional

The polynomial source acts directly in vorticity. A convective Reynolds-stress closure acts through curl of stress divergence and generally carries spatial structure, memory and scale-dependent coefficients. Global mean removal alone does not establish such a closure. The present tests identify how the polynomial organizes and sustains structures; a separate, trustworthy convective-stress prediction is required to establish why the same mechanism operates in the original heated flow.
