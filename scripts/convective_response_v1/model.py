"""Finite-f Boussinesq vertical Galerkin operators, in H and |f| units.

Scientific execution belongs on Unity.  The frozen measured N2 profile is an
input, not a tunable instability parameter.  Thermal perturbations obey
homogeneous Neumann conditions, appropriate to spatially uniform imposed flux.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as la

HEIGHT_M = 120000.0
F_ABS = 0.0002880957739764384
NU_SI = 1000.0


def source_hashes():
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parent.glob("*.py"))}


def read_reference(path):
    path = Path(path)
    d = json.loads(path.read_text())
    if "N2_s2" in d:
        n2 = np.asarray(d["N2_s2"], dtype=float)
        f_input = d.get("f_s_inv", d.get("f_abs_s_inv", F_ABS))
        assert d["source_day"] == 100.0
    else:
        records = d["extremes_by_day"]
        row = next(x for x in records if x["day"] == 100.0)
        n2 = np.asarray(row["mean_N2_profile"], dtype=float)
        f_input = d["f_s_inv"]
    assert len(n2) == 64 and np.isfinite(n2).all()
    assert abs(f_input - F_ABS) < 1e-15
    return n2, {"source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "source_day": 100.0,
                "source_public_path": "gallery/fplane-centrifugal-v1/analysis.json",
                "height_m": HEIGHT_M, "f_abs_s_inv": F_ABS}


def quadrature(ncell=64, order=12):
    # Split at all profile interpolation knots and all taper/relaxation knots.
    knots = np.unique(np.r_[0., (np.arange(ncell)+.5)/ncell, .5, 2/3, .75, 1.])
    xx, ww = np.polynomial.legendre.leggauss(order)
    z = np.concatenate([.5*(b-a)*xx+.5*(a+b) for a, b in zip(knots[:-1], knots[1:])])
    weights = np.concatenate([.5*(b-a)*ww for a, b in zip(knots[:-1], knots[1:])])
    return z, weights


def vertical_tables(nv, nb, z):
    q = np.arange(1, nv+1)*np.pi
    r = np.arange(nb+1)*np.pi
    sine = np.sqrt(2.)*np.sin(z[:, None]*q)
    cosine = np.sqrt(2.)*np.cos(z[:, None]*q)
    beta = np.sqrt(2.)*np.cos(z[:, None]*r)
    beta[:, 0] = 1.
    dbeta = -np.sqrt(2.)*np.sin(z[:, None]*r)*r
    return q, sine, cosine, beta, dbeta


def reference_arrays(n2_si, z, control=False):
    cells = (np.arange(len(n2_si))+.5)/len(n2_si)
    n2 = np.interp(z, cells, n2_si)/(F_ABS**2)
    if control:
        n2 = np.maximum(n2, 0.)
    s = np.clip((z-.5)/(1/6), 0., 1.)
    kappa = (1000. + (1.-1000.)*s*s*(3.-2.*s))/(F_ABS*HEIGHT_M**2)
    s = np.clip((z-2/3)/(1/12), 0., 1.)
    relaxation = s*s*(3.-2.*s)/(F_ABS*86400.)
    return n2, kappa, relaxation


@dataclass
class VerticalSystem:
    """g sine coefficients and h=g'/q; mass matrices are not assumed diagonal."""
    g_coeff: np.ndarray
    nb: int
    n2_si: np.ndarray
    control: bool = False
    thermal_diffusion: bool = True
    thermal_relaxation: bool = True
    nu: float = NU_SI/(F_ABS*HEIGHT_M**2)

    def __post_init__(self):
        self.g_coeff = np.asarray(self.g_coeff, dtype=float)
        self.nv, self.m = self.g_coeff.shape
        self.z, self.weight = quadrature(len(self.n2_si))
        n, sine, cosine, self.beta, dbeta = vertical_tables(self.nv, self.nb, self.z)
        dg_coeff = n[:, None]*self.g_coeff
        self.q = np.sqrt(np.sum(dg_coeff**2, axis=0))
        assert np.all(self.q > 0)
        self.h_coeff = dg_coeff/self.q
        self.g = sine@self.g_coeff
        self.h = cosine@self.h_coeff
        dg = cosine@dg_coeff
        dh = -sine@(n[:, None]*self.h_coeff)
        gram = lambda a, b, v=1.: a.T@(self.weight[:, None]*np.asarray(v).reshape(-1, 1)*b) if np.ndim(v) else a.T@(self.weight[:, None]*v*b)
        self.H = gram(self.h, self.h)
        self.W = gram(self.g, self.g)
        self.Hinv = la.inv(self.H)
        self.Winv = la.inv(self.W)
        self.C = gram(self.g, self.beta)
        n2, kappa, relaxation = reference_arrays(self.n2_si, self.z, self.control)
        self.n2 = n2
        self.B = gram(self.beta, self.g, n2)
        self.Vh = gram(dh, dh)
        self.Vw = gram(dg, dg)
        self.Kx = gram(self.beta, self.beta, kappa) if self.thermal_diffusion else np.zeros((self.nb+1, self.nb+1))
        self.Kz = gram(dbeta, dbeta, kappa) if self.thermal_diffusion else np.zeros_like(self.Kx)
        self.Gamma = gram(self.beta, self.beta, relaxation) if self.thermal_relaxation else np.zeros_like(self.Kx)
        self.M = la.block_diag(self.H, self.H, self.W)
        self.Minv = la.block_diag(self.Hinv, self.Hinv, self.Winv)
        self.beta_mass_error = float(la.norm(gram(self.beta, self.beta)-np.eye(self.nb+1), 2))

    def divergence(self, kx, ky=0.):
        eye = np.eye(self.m)
        return np.hstack((1j*kx*eye, 1j*ky*eye, np.diag(self.q)))

    def projector(self, kx, ky=0.):
        D = self.divergence(kx, ky)
        return np.eye(3*self.m)-self.Minv@D.conj().T@la.solve(D@self.Minv@D.conj().T, D, assume_a="her")

    def full_operator(self, kx, ky=0., f=1.):
        """Projected momentum + Neumann buoyancy coefficient generator."""
        m, nb = self.m, self.nb+1
        k2 = kx*kx+ky*ky
        eye = np.eye(m)
        A = np.zeros((3*m+nb, 3*m+nb), dtype=complex)
        Hdiff = -self.nu*(k2*eye+self.Hinv@self.Vh)
        Wdiff = -self.nu*(k2*eye+self.Winv@self.Vw)
        A[:m, :m] = Hdiff
        A[m:2*m, m:2*m] = Hdiff
        A[2*m:3*m, 2*m:3*m] = Wdiff
        A[:m, m:2*m] = f*eye
        A[m:2*m, :m] = -f*eye
        A[2*m:3*m, 3*m:] = self.Winv@self.C
        A[3*m:, 2*m:3*m] = -self.B
        A[3*m:, 3*m:] = -(k2*self.Kx+self.Kz+self.Gamma)
        A[:3*m] = self.projector(kx, ky)@A[:3*m]
        return A

    def reduced_operator(self, kh, f=1.):
        """Constraint-free state (w, v, b), wavevector (kh, 0)."""
        assert kh > 0
        m, nb = self.m, self.nb+1
        E = np.zeros((3*m+nb, 2*m+nb), dtype=complex)
        E[:m, :m] = 1j*np.diag(self.q)/kh
        E[m:2*m, m:2*m] = np.eye(m)
        E[2*m:3*m, :m] = np.eye(m)
        E[3*m:, 2*m:] = np.eye(nb)
        rows = np.r_[np.arange(2*m, 3*m), np.arange(m, 2*m), np.arange(3*m, 3*m+nb)]
        A = (self.full_operator(kh)@E)[rows]
        if f != 1.:
            A = (self.full_operator(kh, f=f)@E)[rows]
        return A, E

    def spectrum(self, kh, vectors=False, f=1.):
        A, E = self.reduced_operator(kh, f=f)
        if vectors:
            values, states = la.eig(A, check_finite=False)
            order = np.argsort(values.real)[::-1]
            return values[order], states[:, order], E
        return la.eigvals(A, check_finite=False)


def full_system(nv, nb, n2_si, **kwargs):
    return VerticalSystem(np.eye(nv), nb, n2_si, **kwargs)


def real_velocity_basis(system, kh, count=2):
    values, states, _ = system.spectrum(kh, vectors=True)
    basis, info = [], []
    for value, state in zip(values, states.T):
        w = system.g_coeff@state[:system.m]
        if la.norm(w) < 1e-9*la.norm(state):
            continue
        w *= np.exp(-1j*np.angle(w[np.argmax(abs(w))]))
        complex_ratio = float(la.norm(w.imag)/la.norm(w))
        candidates = [w.real] if complex_ratio < 1e-7 else [w.real, w.imag]
        for candidate in candidates:
            candidate = candidate.copy()
            for old in basis:
                candidate -= old@candidate*old
            norm = la.norm(candidate)
            if norm < 1e-7*la.norm(w):
                continue
            basis.append(candidate/norm)
            info.append({"eigenvalue_real": float(value.real), "eigenvalue_imag": float(value.imag),
                         "w_imaginary_fraction_after_phase": complex_ratio,
                         "requires_real_imag_pair": complex_ratio >= 1e-7})
            if len(basis) == count:
                return np.column_stack(basis), info
    raise RuntimeError("Insufficient independent velocity eigenfunctions")


def algebra_checks():
    """Independent operator identities; no fit to production output."""
    rng = np.random.default_rng(20261118)
    n2_si = np.full(64, .7*F_ABS**2)
    coefficients = rng.normal(size=(8, 2))
    coefficients, _ = la.qr(coefficients, mode="economic")
    s = VerticalSystem(coefficients, 18, n2_si, thermal_diffusion=False,
                       thermal_relaxation=False, nu=0.)
    kx, ky = 2.4, -.7
    P = s.projector(kx, ky)
    D = s.divergence(kx, ky)
    p_error = float(la.norm(D@P)/la.norm(D))
    orth_error = float(la.norm(s.M@P-P.conj().T@s.M)/la.norm(s.M))
    idempotence = float(la.norm(P@P-P))
    state = rng.normal(size=3*s.m+s.nb+1)+1j*rng.normal(size=3*s.m+s.nb+1)
    state[:3*s.m] = P@state[:3*s.m]
    derivative = s.full_operator(kx, ky)@state
    energy_power = float(np.real(np.vdot(state[:3*s.m], s.M@derivative[:3*s.m])
                                 + np.vdot(state[3*s.m:], derivative[3*s.m:])/.7))
    scale = float(la.norm(state)**2)
    # Separate exact uniform-N2 Fourier-symbol check, not a thermal-BC check.
    k, q, n2, f = 1.8, 2.7, .7, 1.
    a = np.array([[0, -1j*k*q*f/(k*k+q*q), k*k/(k*k+q*q)],
                  [-1j*f*q/k, 0, 0], [-n2, 0, 0]], dtype=complex)
    eig = la.eigvals(a)
    omega = np.sqrt((n2*k*k+f*f*q*q)/(k*k+q*q))
    dispersion = max(min(abs(eig-x)) for x in (0, 1j*omega, -1j*omega))
    checks = {"continuity_relative": p_error, "pressure_metric_orthogonality": orth_error,
              "pressure_idempotence": idempotence, "constant_N2_energy_relative": abs(energy_power)/scale,
              "uniform_fourier_dispersion_absolute": float(dispersion),
              "thermal_cosine_mass_error": s.beta_mass_error}
    assert max(checks.values()) < 2e-10, checks
    return checks
