"""Complete vertical eigenmode subspaces, retaining toroidal structure.

This supplements the frozen full-reference screen. One selected eigenmode may
require two real functions when its vertical eigenfunction is genuinely complex.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as la

from model import (F_ABS, HEIGHT_M, NU_SI, quadrature, reference_arrays,
                   vertical_tables)


def append_orthogonal(columns, vector, relative=1e-9):
    vector = np.asarray(vector, dtype=float).copy()
    original = la.norm(vector)
    if original < 1e-14:
        return False
    # Reorthogonalization avoids a near-dependent derivative/toroidal direction.
    for _ in range(2):
        for old in columns:
            vector -= np.dot(old, vector)*old
    norm = la.norm(vector)
    if norm < relative*original:
        return False
    columns.append(vector/norm)
    return True


def complete_eigenmode_bases(reference, kh, counts=(1, 2)):
    """Nested spaces containing selected full-reference velocity eigenvectors."""
    values, states, _ = reference.spectrum(kh, vectors=True)
    g_columns, h_columns, accepted, snapshots = [], [], [], {}
    q = np.arange(1, reference.nv+1)*np.pi
    selected_values = []
    for value, state in zip(values, states.T):
        # Skip the conjugate representation of the same complex normal mode.
        if any(abs(value-v.conjugate()) < 1e-7*max(1., abs(value)) for v in selected_values):
            continue
        w = reference.g_coeff@state[:reference.m]
        v = reference.h_coeff@state[reference.m:2*reference.m]
        if la.norm(w) < 1e-9*la.norm(state):
            continue
        phase = np.exp(-1j*np.angle(w[np.argmax(abs(w))]))
        w, v = w*phase, v*phase
        added = append_orthogonal(g_columns, w.real)
        added |= append_orthogonal(g_columns, w.imag)
        if not added:
            continue
        # Exact continuity requires all derivatives; toroidal profiles retain
        # the full Coriolis-coupled eigenvector rather than guessing v~g'.
        for g in g_columns:
            append_orthogonal(h_columns, q*g)
        append_orthogonal(h_columns, v.real)
        append_orthogonal(h_columns, v.imag)
        selected_values.append(value)
        accepted.append({"eigenvalue_real": float(value.real),
                         "eigenvalue_imag": float(value.imag),
                         "w_imaginary_fraction_after_phase": float(la.norm(w.imag)/la.norm(w)),
                         "cumulative_vertical_velocity_functions": len(g_columns),
                         "cumulative_horizontal_velocity_functions": len(h_columns)})
        n = len(accepted)
        if n in counts:
            snapshots[n] = {"g_coeff": np.column_stack(g_columns),
                            "h_coeff": np.column_stack(h_columns),
                            "selected_modes": list(accepted)}
        if n == max(counts):
            return snapshots
    raise RuntimeError("Could not select independent complete vertical modes")


@dataclass
class AugmentedSystem:
    g_coeff: np.ndarray
    h_coeff: np.ndarray
    nb: int
    n2_si: np.ndarray
    control: bool = False
    thermal_diffusion: bool = True
    thermal_relaxation: bool = True
    nu: float = NU_SI/(F_ABS*HEIGHT_M**2)

    def __post_init__(self):
        self.g_coeff = np.asarray(self.g_coeff, dtype=float)
        self.h_coeff = np.asarray(self.h_coeff, dtype=float)
        self.nv, self.nw = self.g_coeff.shape
        assert self.h_coeff.shape[0] == self.nv
        self.nh = self.h_coeff.shape[1]
        self.z, self.weight = quadrature(len(self.n2_si))
        n, sine, cosine, self.beta, dbeta = vertical_tables(self.nv, self.nb, self.z)
        self.g = sine@self.g_coeff
        self.h = cosine@self.h_coeff
        dg = cosine@(n[:, None]*self.g_coeff)
        dh = -sine@(n[:, None]*self.h_coeff)

        def gram(a, b, factor=1.):
            return a.T@((self.weight*factor)[:, None]*b)

        self.H = gram(self.h, self.h)
        self.W = gram(self.g, self.g)
        self.Hinv, self.Winv = la.inv(self.H), la.inv(self.W)
        self.J = self.Hinv@gram(self.h, dg)
        self.derivative_span_error = float(la.norm(dg-self.h@self.J)/la.norm(dg))
        assert self.derivative_span_error < 1e-8
        self.C = gram(self.g, self.beta)
        n2, kappa, relaxation = reference_arrays(self.n2_si, self.z, self.control)
        self.B = gram(self.beta, self.g, n2)
        self.Vh, self.Vw = gram(dh, dh), gram(dg, dg)
        zero = np.zeros((self.nb+1, self.nb+1))
        self.Kx = gram(self.beta, self.beta, kappa) if self.thermal_diffusion else zero.copy()
        self.Kz = gram(dbeta, dbeta, kappa) if self.thermal_diffusion else zero.copy()
        self.Gamma = gram(self.beta, self.beta, relaxation) if self.thermal_relaxation else zero.copy()
        self.M = la.block_diag(self.H, self.H, self.W)
        self.Minv = la.block_diag(self.Hinv, self.Hinv, self.Winv)
        self.nvelocity = 2*self.nh+self.nw

    def divergence(self, kx, ky=0.):
        eye = np.eye(self.nh)
        return np.hstack((1j*kx*eye, 1j*ky*eye, self.J))

    def projector(self, kx, ky=0.):
        if kx == 0 and ky == 0:
            return la.block_diag(np.eye(self.nh), np.eye(self.nh), np.zeros((self.nw, self.nw)))
        D = self.divergence(kx, ky)
        return np.eye(self.nvelocity)-self.Minv@D.conj().T@la.solve(D@self.Minv@D.conj().T, D, assume_a="her")

    def full_operator(self, kx, ky=0., f=1.):
        nh, nw, nb = self.nh, self.nw, self.nb+1
        nv = self.nvelocity
        k2 = kx*kx+ky*ky
        A = np.zeros((nv+nb, nv+nb), dtype=complex)
        Hdiff = -self.nu*(k2*np.eye(nh)+self.Hinv@self.Vh)
        Wdiff = -self.nu*(k2*np.eye(nw)+self.Winv@self.Vw)
        A[:nh, :nh] = Hdiff
        A[nh:2*nh, nh:2*nh] = Hdiff
        A[2*nh:nv, 2*nh:nv] = Wdiff
        A[:nh, nh:2*nh] = f*np.eye(nh)
        A[nh:2*nh, :nh] = -f*np.eye(nh)
        A[2*nh:nv, nv:] = self.Winv@self.C
        A[nv:, 2*nh:nv] = -self.B
        A[nv:, nv:] = -(k2*self.Kx+self.Kz+self.Gamma)
        A[:nv] = self.projector(kx, ky)@A[:nv]
        return A

    def reduced_operator(self, kh, f=1.):
        assert kh > 0
        nh, nw, nb = self.nh, self.nw, self.nb+1
        nv = self.nvelocity
        E = np.zeros((nv+nb, nw+nh+nb), dtype=complex)
        E[:nh, :nw] = 1j*self.J/kh
        E[nh:2*nh, nw:nw+nh] = np.eye(nh)
        E[2*nh:nv, :nw] = np.eye(nw)
        E[nv:, nw+nh:] = np.eye(nb)
        rows = np.r_[np.arange(2*nh, nv), np.arange(nh, 2*nh), np.arange(nv, nv+nb)]
        return (self.full_operator(kh, f=f)@E)[rows], E

    def spectrum(self, kh, vectors=False, f=1.):
        A, E = self.reduced_operator(kh, f=f)
        if vectors:
            values, states = la.eig(A, check_finite=False)
            order = np.argsort(values.real)[::-1]
            return values[order], states[:, order], E
        return la.eigvals(A, check_finite=False)
