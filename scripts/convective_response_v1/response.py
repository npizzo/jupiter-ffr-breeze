"""Dealiased horizontal response to a fixed vortex; scientific runs on Unity."""
from __future__ import annotations

import numpy as np
import scipy.fft as fft
import scipy.linalg as la

from model import HEIGHT_M
from optimized import AugmentedSystem


def load_basis(path):
    with np.load(path) as data:
        nb = int(data["thermal_nb"])
        systems = {m: AugmentedSystem(data[f"g{m}"], data[f"h{m}"], nb, data["n2_si"])
                   for m in (1, 2)}
        kh = float(data["kh_H"])
    return systems, {"kh_seed": kh, "length": 16*2*np.pi/kh,
                     "radius": 2*2*np.pi/kh}


def with_control(system):
    return AugmentedSystem(system.g_coeff, system.h_coeff, system.nb,
                           system.n2_si, control=True)


class Response:
    def __init__(self, system, n=128, length=None, radius=None, vortex_sign=1,
                 f=1., core_rossby=.25, kh_seed=None):
        assert n % 2 == 0 and vortex_sign in (-1, 1) and abs(f) == 1.
        self.system, self.n = system, n
        self.length, self.radius = float(length), float(radius)
        self.f, self.vortex_sign = float(f), vortex_sign
        self.kh_seed = float(kh_seed)
        self.nh, self.nw = system.nh, system.nw
        self.nv = system.nvelocity
        self.nvar = self.nv+system.nb+1
        self.dx = self.length/n
        self.x = self.y = (np.arange(n)-n/2)*self.dx
        self.X, self.Ygrid = np.meshgrid(self.x, self.y, indexing="ij")
        self.kx = (2*np.pi*fft.fftfreq(n, d=self.dx))[:, None]
        self.ky = (2*np.pi*fft.rfftfreq(n, d=self.dx))[None, :]
        self.k2 = self.kx**2+self.ky**2
        ix = fft.fftfreq(n)*n
        iy = fft.rfftfreq(n)*n
        self.mask = (abs(ix[:, None]) < n/3)&(abs(iy[None, :]) < n/3)
        self.weights = np.ones((1, n//2+1))
        self.weights[:, 1:-1] = 2.
        self.normalization = float(n**4)
        # A periodic Gaussian sum, then the same spectral truncation as dynamics.
        gauss = sum(np.exp(-(self.x+a*self.length)**2/(2*self.radius**2)) for a in (-1, 0, 1))
        psihat = self.spectral(-gauss[:, None]*gauss[None, :])*self.mask
        z0 = self.physical(-self.k2*psihat)
        self.gaussian_amplitude = core_rossby/float(z0.max())*vortex_sign
        psihat *= self.gaussian_amplitude
        self.Uhat, self.Vhat = -1j*self.ky*psihat, 1j*self.kx*psihat
        self.U, self.V = self.physical(self.Uhat), self.physical(self.Vhat)
        self.Z = self.physical(-self.k2*psihat)
        self.Ux, self.Uy = self.physical(1j*self.kx*self.Uhat), self.physical(1j*self.ky*self.Uhat)
        self.Vx, self.Vy = self.physical(1j*self.kx*self.Vhat), self.physical(1j*self.ky*self.Vhat)
        self.base_rms = float(np.sqrt(np.mean(self.U**2+self.V**2)))
        self.velocity_mass = np.asarray(system.M, dtype=float)
        self.H = np.asarray(system.H, dtype=float)
        self.W = np.asarray(system.W, dtype=float)
        self.Hinv = np.asarray(system.Hinv, dtype=float)
        self.Winv = np.asarray(system.Winv, dtype=float)
        self.J = np.asarray(system.J, dtype=float)
        self.C = np.asarray(system.C, dtype=float)
        self.B = np.asarray(system.B, dtype=float)
        self.Kx = np.asarray(system.Kx, dtype=float)
        self.Kz = np.asarray(system.Kz, dtype=float)
        self.Gamma = np.asarray(system.Gamma, dtype=float)
        self.Hlap = np.asarray(system.Hinv@system.Vh, dtype=float)
        self.Wlap = np.asarray(system.Winv@system.Vw, dtype=float)
        self.b_to_w = np.asarray(system.Winv@system.C, dtype=float)
        # Solve k² H^-1 + J W^-1 J^T through a fixed small generalized eigensystem.
        chol = la.cholesky(system.Hinv, lower=True)
        ci = la.inv(chol)
        matrix = ci@(system.J@system.Winv@system.J.T)@ci.T
        mu, vec = la.eigh(matrix)
        mu[np.abs(mu) < 1e-11*max(1., max(abs(mu)))] = 0.
        assert min(mu) > -1e-9
        mu = np.maximum(mu, 0.)
        self.pressure_T = np.asarray(ci.T@vec, dtype=float)
        denom = mu[:, None, None]+self.k2
        self.pressure_inverse = np.divide(1., denom, out=np.zeros_like(denom), where=denom > 0)

    def physical(self, field):
        return fft.irfft2(field, s=(self.n, self.n), axes=(-2, -1), workers=1)

    def spectral(self, field):
        return fft.rfft2(field, axes=(-2, -1), workers=1)

    @staticmethod
    def multiply(matrix, field):
        flat = np.ascontiguousarray(field.reshape(field.shape[0], -1))
        if field.dtype == np.complex128 and not np.iscomplexobj(matrix):
            # Exact real-coefficient complex multiplication via one DGEMM on
            # interleaved real/imaginary columns; avoid redundant complex BLAS.
            product = np.ascontiguousarray(matrix@flat.view(np.float64)).view(np.complex128)
        else:
            product = matrix@flat
        return product.reshape((matrix.shape[0],)+field.shape[1:])

    def inner(self, a, b):
        return float(np.real(np.sum(np.conj(a)*b*self.weights))/self.normalization)

    def project_velocity(self, velocity):
        nh, nv = self.nh, self.nv
        result = velocity.copy()
        div = 1j*self.kx*result[:nh]+1j*self.ky*result[nh:2*nh]+self.multiply(self.J, result[2*nh:nv])
        pressure = self.multiply(self.pressure_T,
                                self.pressure_inverse*self.multiply(self.pressure_T.T, div))
        hp = self.multiply(self.Hinv, pressure)
        result[:nh] += 1j*self.kx*hp
        result[nh:2*nh] += 1j*self.ky*hp
        result[2*nh:nv] -= self.multiply(self.Winv@self.J.T, pressure)
        result[2*nh:nv, 0, 0] = 0.
        return result*self.mask

    def project(self, state):
        result = state.copy()*self.mask
        result[:self.nv] = self.project_velocity(result[:self.nv])
        return result

    def energy(self, state):
        return {"kinetic": .5*self.inner(state[:self.nv], self.multiply(self.velocity_mass, state[:self.nv])),
                "thermal_variance": .5*self.inner(state[self.nv:], state[self.nv:])}

    def seed(self, seed, amplitude_ratio=1e-6):
        rng = np.random.default_rng(seed)
        shell = (self.k2 >= (.8*self.kh_seed)**2)&(self.k2 <= (1.2*self.kh_seed)**2)&self.mask
        psi = self.spectral(rng.standard_normal((self.n, self.n)))*shell/self.kh_seed
        w = self.spectral(rng.standard_normal((self.n, self.n)))*shell
        state = np.zeros((self.nvar, self.n, self.n//2+1), dtype=complex)
        state[2*self.nh] = w
        jw = self.multiply(self.J, state[2*self.nh:self.nv])
        invk2 = np.divide(1., self.k2, out=np.zeros_like(self.k2), where=self.k2 > 0)
        state[:self.nh] = 1j*self.kx*invk2*jw
        state[self.nh:2*self.nh] = 1j*self.ky*invk2*jw
        state[0] += -1j*self.ky*psi
        state[self.nh] += 1j*self.kx*psi
        state = self.project(state)
        state *= amplitude_ratio*self.base_rms/np.sqrt(2*self.energy(state)["kinetic"])
        return state

    def reflect_state(self, state):
        reflected = state[:, (-np.arange(self.n)) % self.n, :].copy()
        reflected[:self.nh] *= -1
        return reflected

    def rhs(self, state, with_work=False):
        nh, nv = self.nh, self.nv
        u, v, w, b = state[:nh], state[nh:2*nh], state[2*nh:nv], state[nv:]
        # Pseudospectral Galerkin advection; all factors and output use 2/3 mask.
        adv = -self.U*self.physical(1j*self.kx*state)-self.V*self.physical(1j*self.ky*state)
        hp = self.physical(state[:2*nh])
        up, vp = hp[:nh], hp[nh:]
        adv[:nh] -= self.Ux*up+self.Uy*vp
        adv[nh:2*nh] -= self.Vx*up+self.Vy*vp
        tendency = self.spectral(adv)*self.mask
        du = -self.system.nu*(self.k2*u+self.multiply(self.Hlap, u))
        dv = -self.system.nu*(self.k2*v+self.multiply(self.Hlap, v))
        dw = -self.system.nu*(self.k2*w+self.multiply(self.Wlap, w))
        buoy = self.multiply(self.b_to_w, b)
        thermal_gradient = -self.multiply(self.B, w)
        thermal_diffusion = -(self.k2*self.multiply(self.Kx, b)+self.multiply(self.Kz, b))
        thermal_relaxation = -self.multiply(self.Gamma, b)
        tendency[:nh] += du+self.f*v
        tendency[nh:2*nh] += dv-self.f*u
        tendency[2*nh:nv] += dw+buoy
        tendency[nv:] += thermal_gradient+thermal_diffusion+thermal_relaxation
        tendency[:nv] = self.project_velocity(tendency[:nv])
        tendency *= self.mask
        if not with_work:
            return tendency
        Hu = self.multiply(self.system.H, up)
        Hv = self.multiply(self.system.H, vp)
        rxx = np.sum(up*Hu, axis=0)
        rxy = np.sum(up*Hv, axis=0)
        ryy = np.sum(vp*Hv, axis=0)
        powers = {"bt_power": float(np.mean(rxx*self.Ux+rxy*(self.Uy+self.Vx)+ryy*self.Vy)),
                  "buoyancy_power": self.inner(w, self.multiply(self.C, b)),
                  "viscous_power": self.inner(u, self.multiply(self.H, du))+self.inner(v, self.multiply(self.H, dv))+self.inner(w, self.multiply(self.W, dw)),
                  "thermal_gradient_power": self.inner(b, thermal_gradient),
                  "thermal_diffusion_power": self.inner(b, thermal_diffusion),
                  "thermal_relaxation_power": self.inner(b, thermal_relaxation)}
        return tendency, powers

    def rk4(self, state, dt):
        k1, p1 = self.rhs(state, True)
        k2, p2 = self.rhs(state+.5*dt*k1, True)
        k3, p3 = self.rhs(state+.5*dt*k2, True)
        k4, p4 = self.rhs(state+dt*k3, True)
        result = state+dt/6*(k1+2*k2+2*k3+k4)
        result *= self.mask
        work = {k: dt/6*(p1[k]+2*p2[k]+2*p3[k]+p4[k]) for k in p1}
        return result, work

    def diagnostics(self, state):
        rhs, powers = self.rhs(state, True)
        energy = self.energy(state)
        kinetic_rate = self.inner(state[:self.nv], self.multiply(self.velocity_mass, rhs[:self.nv]))
        thermal_rate = self.inner(state[self.nv:], rhs[self.nv:])
        div = 1j*self.kx*state[:self.nh]+1j*self.ky*state[self.nh:2*self.nh]+self.multiply(self.J, state[2*self.nh:self.nv])
        velocity_scale = np.sqrt(max(self.inner(state[:self.nv], state[:self.nv]), np.finfo(float).tiny))
        return {**energy, **powers,
                "kinetic_rhs_rate": kinetic_rate,
                "kinetic_identity_residual": kinetic_rate-(-powers["bt_power"]+powers["buoyancy_power"]+powers["viscous_power"]),
                "thermal_rhs_rate": thermal_rate,
                "thermal_identity_residual": thermal_rate-sum(powers[k] for k in ("thermal_gradient_power", "thermal_diffusion_power", "thermal_relaxation_power")),
                "divergence_relative": np.sqrt(max(self.inner(div, div), 0.))/velocity_scale,
                "perturbation_rms_over_base": np.sqrt(2*energy["kinetic"])/self.base_rms}

    def plane_fields(self, state, heights_m=(20000., 60000., 100000.)):
        z = np.asarray(heights_m)/HEIGHT_M
        n = np.arange(1, self.system.nv+1)*np.pi
        h = (np.sqrt(2.)*np.cos(z[:, None]*n))@self.system.h_coeff
        r = np.arange(self.system.nb+1)*np.pi
        beta = np.sqrt(2.)*np.cos(z[:, None]*r)
        beta[:, 0] = 1.
        qcoeff = self.physical(1j*self.kx*state[self.nh:2*self.nh]-1j*self.ky*state[:self.nh])
        zeta = self.multiply(h, qcoeff)
        buoyancy = self.multiply(beta, self.physical(state[self.nv:]))
        return {"heights_m": np.asarray(heights_m), "perturbation_zeta_over_fabs": zeta.astype(np.float32),
                "total_zeta_over_fabs": (zeta+self.Z).astype(np.float32),
                "buoyancy_over_fabs2_H": buoyancy.astype(np.float32)}

    def stable_step(self, ceiling=.1):
        kmax = float(np.sqrt(np.max(self.k2[self.mask])))
        samples = np.geomspace(2*np.pi/self.length, kmax, 10)
        values = np.concatenate([self.system.spectrum(float(k), f=self.f) for k in samples])
        cfl_rate = np.max(abs(self.U)+abs(self.V))/self.dx
        dt = min(ceiling, .4/max(cfl_rate, 1e-15))
        while True:
            z = dt*values
            amplification = 1+z+z*z/2+z*z*z/6+z*z*z*z/24
            if np.all(abs(amplification) <= np.exp(np.maximum(z.real, 0))*(1+1e-9)):
                break
            dt *= .8
            if dt < 1e-5:
                raise RuntimeError("Explicit response timestep impractical")
        # A divisor of .1 aligns both .2 scalar and .5 frame events. This avoids
        # extra short steps while retaining the independently tested stability
        # limit; the response validation still gates temporal accuracy.
        dt = .1/max(1, int(np.ceil(.1/dt-1e-12)))
        return float(dt)
