"""Nonlinear, energy-consistent Fourier--Galerkin rotating convection.

Independent implementation of the Boussinesq equations, using the mode-selection
principles of Rademacher & Welter, arXiv:2505.07678.  This is their standard
fixed-temperature, stress-free reference layer, not the capped anelastic Breeze
atmosphere.  Scientific execution and rendering belong on Unity.

Units: depth H=1; time is an externally chosen reference inverse Coriolis rate.
The conductive buoyancy is B0=-G*z.  Its perturbation theta vanishes at z=0,1.
The state contains normalized Fourier coefficients of (u,v,w,theta) on the
even/odd extension to 0<=z<2.  Horizontal coordinates have period 2*pi/kfund.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


def _even_above(value):
    """Small even integer strictly greater than value."""
    return 2 * (int(value) // 2 + 1)


@dataclass
class GalerkinModel:
    K: int = 3
    Nv: int = 2
    kfund: float = 1.0
    nu: float = 0.005
    kappa: float = 0.005
    G: float = 1.0
    f: float = 1.0
    grid_shape: tuple | None = None

    def __post_init__(self):
        assert self.K >= 1 and self.Nv >= 1
        assert self.kfund > 0 and self.nu >= 0 and self.kappa >= 0
        assert self.G >= 0 and np.isfinite(self.f)
        minimum = (_even_above(3*self.K), _even_above(3*self.K),
                   _even_above(6*self.Nv))
        self.shape = tuple(self.grid_shape or minimum)
        assert len(self.shape) == 3
        assert self.shape[0] > 3*self.K and self.shape[1] > 3*self.K
        assert self.shape[2] > 6*self.Nv
        self.ix, self.iy, self.iz = [np.rint(np.fft.fftfreq(n)*n).astype(int)
                                   for n in self.shape]
        self.neg = tuple((-np.arange(n)) % n for n in self.shape)
        ix, iy, iz = np.meshgrid(self.ix, self.iy, self.iz, indexing="ij")
        self.k = np.stack((self.kfund*ix, self.kfund*iy, np.pi*iz))
        self.k2 = np.sum(self.k*self.k, axis=0)
        self.inv_k2 = np.divide(1., self.k2, out=np.zeros(self.shape),
                                where=self.k2 != 0)
        hmask = (np.abs(ix) <= self.K) & (np.abs(iy) <= self.K)
        hzero = (ix == 0) & (iy == 0)
        # Mean shear and mean thermal modes are essential balance-law partners.
        vertical = np.where(hzero, np.abs(iz) <= 2*self.Nv,
                            np.abs(iz) <= self.Nv)
        self.mask_u = hmask & vertical
        self.mask_u[0, 0, 0] = False
        self.mask_theta = hmask & vertical & (iz != 0)
        self._parity = np.asarray([1., 1., -1., -1.])[:, None, None, None]
        self.x = np.arange(self.shape[0]) * (2*np.pi/self.kfund/self.shape[0])
        self.y = np.arange(self.shape[1]) * (2*np.pi/self.kfund/self.shape[1])
        self.z = np.arange(self.shape[2]) * (2/self.shape[2])

    @staticmethod
    def spectral(fields):
        return np.fft.fftn(fields, axes=(-3, -2, -1), norm="forward")

    @staticmethod
    def physical(state):
        return np.fft.ifftn(state, axes=(-3, -2, -1), norm="forward").real

    def _symmetries(self, state):
        """Orthogonal reality and vertical-parity projections."""
        partner = state[:, self.neg[0]][:, :, self.neg[1]][:, :, :, self.neg[2]]
        out = 0.5*(state + partner.conj())
        return 0.5*(out + self._parity*out[:, :, :, self.neg[2]])

    def velocity_project(self, velocity):
        out = np.asarray(velocity, dtype=complex).copy()
        out -= self.k*(np.sum(self.k*out, axis=0)*self.inv_k2)[None, ...]
        out *= self.mask_u[None, ...]
        return out

    def project(self, state):
        assert np.shape(state) == (4,)+self.shape
        out = self._symmetries(np.asarray(state, dtype=complex))
        out[:3] = self.velocity_project(out[:3])
        out[3] *= self.mask_theta
        return out

    def _advective(self, state, fields=None):
        fields = self.physical(state) if fields is None else fields
        # All products are evaluated on a quadrature grid large enough that no
        # discarded quadratic product can alias into a retained Fourier mode.
        gradients = self.physical(1j*self.k[:, None, ...]*state[None, ...])
        adv = -np.einsum("jxyz,jcxyz->cxyz", fields[:3], gradients,
                         optimize=False)
        return self.project(self.spectral(adv))

    def tendencies(self, state, nonlinear=True):
        """Separate actual RHS terms, all in normalized Fourier coefficients."""
        state = self.project(state)
        terms = {}
        terms["advective"] = self._advective(state) if nonlinear else np.zeros_like(state)
        rotation = np.zeros_like(state)
        rotation[0] = self.f*state[1]
        rotation[1] = -self.f*state[0]
        terms["coriolis"] = self.project(rotation)
        buoyancy = np.zeros_like(state)
        buoyancy[2] = state[3]
        terms["buoyancy"] = self.project(buoyancy)
        diffusion = np.empty_like(state)
        diffusion[:3] = -self.nu*self.k2[None, ...]*state[:3]
        diffusion[3] = -self.kappa*self.k2*state[3]
        terms["diffusion"] = diffusion
        gradient = np.zeros_like(state)
        gradient[3] = self.G*state[2]
        terms["background_gradient"] = self.project(gradient)
        return terms

    def rhs(self, state, nonlinear=True):
        """Same terms as tendencies(); input must be in the projected subspace."""
        # A projected initial state and every projected RK stage keep the input
        # in the selected real, divergence-free subspace.
        out = self._advective(state) if nonlinear else np.zeros_like(state)
        raw = np.zeros_like(state)
        raw[0] = self.f*state[1]
        raw[1] = -self.f*state[0]
        raw[2] = state[3]
        out[:3] += self.velocity_project(raw[:3])
        out[:3] -= self.nu*self.k2[None, ...]*state[:3]
        out[3] += self.G*state[2] - self.kappa*self.k2*state[3]
        return out

    def rk4_step(self, state, dt, nonlinear=True, return_stages=False):
        stages = []
        k1 = self.rhs(state, nonlinear)
        if return_stages:
            stages.append((state, k1, 1/6))
        s2 = self.project(state + .5*dt*k1)
        k2 = self.rhs(s2, nonlinear)
        if return_stages:
            stages.append((s2, k2, 1/3))
        s3 = self.project(state + .5*dt*k2)
        k3 = self.rhs(s3, nonlinear)
        if return_stages:
            stages.append((s3, k3, 1/3))
        s4 = self.project(state + dt*k3)
        k4 = self.rhs(s4, nonlinear)
        if return_stages:
            stages.append((s4, k4, 1/6))
        result = self.project(state + dt*(k1+2*k2+2*k3+k4)/6)
        return (result, stages) if return_stages else result

    def scalar_rates(self, state, total_rhs):
        """Budget rates using an already computed RHS, without additional FFTs.

        Barotropic work from buoyancy and projected Coriolis acceleration is
        zero. The nonlinear exchange is consequently its actual total work
        minus viscous work. Independently computed stress work is available in
        diagnostics(). The same decomposition is used for thermal variance.
        """
        kin_rate = float(np.vdot(state[:3], total_rhs[:3]).real)
        theta_rate = float(np.vdot(state[3], total_rhs[3]).real)
        bt_rate = float(np.vdot(state[:2, :, :, 0], total_rhs[:2, :, :, 0]).real)
        heat = float(np.vdot(state[2], state[3]).real)
        visc = float(-self.nu*np.sum(self.k2[None, ...]*abs(state[:3])**2))
        bt_visc = float(-self.nu*np.sum(self.k2[None, :, :, 0]
                                        *abs(state[:2, :, :, 0])**2))
        theta_diff = float(-self.kappa*np.sum(self.k2*abs(state[3])**2))
        moment_rate = self.potential_moment(total_rhs)
        moment_diff = self.potential_moment(-self.kappa*self.k2*state[3])
        return {
            "kinetic": kin_rate,
            "thermal_variance": theta_rate,
            "potential_moment": moment_rate,
            "barotropic_kinetic": bt_rate,
            "barotropic_advective": bt_rate-bt_visc,
            "barotropic_viscous": bt_visc,
            "buoyancy": heat,
            "viscous": visc,
            "kinetic_advective_and_coriolis": kin_rate-heat-visc,
            "thermal_gradient": self.G*heat,
            "thermal_diffusive": theta_diff,
            "thermal_advective": theta_rate-self.G*heat-theta_diff,
            "potential_advective": moment_rate-moment_diff,
            "potential_diffusive": moment_diff,
        }

    def rk4_step_with_work(self, state, dt, nonlinear=True):
        """Return next state and dt-integrated RK-stage budget terms.

        A quadratic-energy change differs from its RK quadrature by a finite
        time-integration defect; retain that defect rather than force closure.
        The linear potential moment closes to floating-point accuracy.
        """
        result, stages = self.rk4_step(state, dt, nonlinear, return_stages=True)
        work = {}
        for stage, tendency, weight in stages:
            for key, value in self.scalar_rates(stage, tendency).items():
                work[key] = work.get(key, 0.) + dt*weight*value
        return result, work

    @staticmethod
    def state_energy(state):
        """Physical-layer volume-mean kinetic energy (Parseval normalization)."""
        return float(.5*np.vdot(state[:3], state[:3]).real)

    @staticmethod
    def thermal_variance(state):
        return float(.5*np.vdot(state[3], state[3]).real)

    def potential_moment(self, state_or_theta):
        """Exact <(1-z)*theta> on 0<=z<=1, not on the doubled extension."""
        theta = state_or_theta[3] if np.ndim(state_or_theta) == 4 else state_or_theta
        positive = self.iz > 0
        return float(np.real(np.sum(2j*theta[0, 0, positive]/(np.pi*self.iz[positive]))))

    def reflect(self, state):
        """x -> -x with ux -> -ux. Pair with a model whose f is reversed."""
        out = state[:, self.neg[0], :, :].copy()
        out[0] *= -1
        return self.project(out)

    def initial_state(self, seed=20261118, velocity_rms=1e-3,
                      temperature_rms=1e-3, barotropic=False):
        """Matched broad modal disturbances; no continuing random forcing."""
        rng = np.random.default_rng(seed)
        state = self.project(self.spectral(rng.normal(size=(4,)+self.shape)))
        if not barotropic:
            state[:3, :, :, 0] = 0
        # Do not seed horizontally averaged thermal/shear changes. They evolve.
        state[:, 0, 0, :] = 0
        vrms = math.sqrt(2*self.state_energy(state))
        trms = math.sqrt(2*self.thermal_variance(state))
        state[:3] *= velocity_rms/vrms if vrms else 0
        state[3] *= temperature_rms/trms if trms else 0
        return self.project(state)

    def diagnostics(self, state):
        state = self.project(state)
        terms = self.tendencies(state)
        total_rhs = sum(terms.values())
        fields = self.physical(state)
        kin = self.state_energy(state)
        kt = float(.5*np.vdot(state[:2, :, :, 0], state[:2, :, :, 0]).real)
        variance = self.thermal_variance(state)
        heat = float(np.vdot(state[2], state[3]).real)
        powers = {name: float(np.vdot(state[:3], value[:3]).real)
                  for name, value in terms.items()}
        thermal_powers = {name: float(np.vdot(state[3], value[3]).real)
                          for name, value in terms.items()}
        pbt = float(np.vdot(state[:2, :, :, 0], terms["advective"][:2, :, :, 0]).real)
        pbt_diff = float(np.vdot(state[:2, :, :, 0], terms["diffusion"][:2, :, :, 0]).real)
        # Independently evaluate the stress-strain form of barotropic work.
        mean = fields[:2].mean(axis=-1)
        fluctuation = fields[:2] - mean[..., None]
        stress = np.einsum("ixyz,jxyz->ijxy", fluctuation, fluctuation)/self.shape[2]
        bt = state[:2, :, :, 0]
        gradient = np.fft.ifftn(1j*self.k[:2, None, :, :, 0]*bt[None, ...],
                                axes=(-2, -1), norm="forward").real
        pbt_stress = float(np.mean(np.einsum("ijxy,jixy->xy", stress, gradient)))
        zeta_hat = 1j*(self.k[0]*state[1]-self.k[1]*state[0])
        zeta = self.physical(zeta_hat)
        zeta_bt = np.fft.ifftn(zeta_hat[:, :, 0], axes=(-2, -1), norm="forward").real
        # The cubic vorticity moment is not a quadratic Galerkin invariant:
        # retain the actual projected barotropic self-advection contribution.
        bt_only = np.zeros_like(state)
        bt_only[:2, :, :, 0] = state[:2, :, :, 0]
        bt_self_advective = self._advective(bt_only)
        def third_moment_rate(tendency):
            zeta_rhs_hat = 1j*(self.k[0, :, :, 0]*tendency[1, :, :, 0]
                                - self.k[1, :, :, 0]*tendency[0, :, :, 0])
            zeta_rhs = np.fft.ifftn(zeta_rhs_hat, axes=(-2, -1), norm="forward").real
            return float(3*np.mean(zeta_bt*zeta_bt*zeta_rhs))
        third_full_advective = third_moment_rate(terms["advective"])
        third_bt_self = third_moment_rate(bt_self_advective)
        horizontal_mode2 = self.ix[:, None]**2+self.iy[None, :]**2
        lowest = (horizontal_mode2 > 0) & (horizontal_mode2 <= 2)
        outer = ((abs(self.ix[:, None]) == self.K)
                 | (abs(self.iy[None, :]) == self.K))
        mode_energy = .5*np.sum(abs(state[:3])**2, axis=0)
        lowest_energy = float(np.sum(mode_energy[lowest]))
        lowest_bt_energy = float(np.sum(mode_energy[:, :, 0][lowest]))
        outer_energy = float(np.sum(mode_energy[outer]))
        def skew(array):
            second = float(np.mean(array*array))
            return float(np.mean(array**3)/second**1.5) if second > 1e-30 else 0.
        div = np.sum(1j*self.k*state[:3], axis=0)
        grad_norm = float(np.sqrt(np.sum(self.k2[None, ...]*abs(state[:3])**2)))
        # Boundary conduction fluxes from the actual evolving mean temperature.
        derivative_bottom = float(np.real(np.sum(1j*np.pi*self.iz*state[3, 0, 0])))
        derivative_top = float(np.real(np.sum(1j*np.pi*self.iz*state[3, 0, 0]
                                              *np.exp(1j*np.pi*self.iz))))
        return {
            "kinetic": kin, "barotropic_kinetic": kt,
            "baroclinic_kinetic": kin-kt,
            "barotropic_energy_fraction": kt/kin if kin else 0.,
            "thermal_variance": variance,
            "potential_moment": self.potential_moment(state),
            "buoyancy_power": heat, "barotropic_power": pbt,
            "barotropic_stress_power": pbt_stress,
            "barotropic_viscous_power": pbt_diff,
            "kinetic_tendency": float(np.vdot(state[:3], total_rhs[:3]).real),
            "thermal_variance_tendency": float(np.vdot(state[3], total_rhs[3]).real),
            "potential_moment_tendency": self.potential_moment(total_rhs),
            "kinetic_advective_power": powers["advective"],
            "coriolis_power": powers["coriolis"],
            "viscous_power": powers["diffusion"],
            "thermal_advective_power": thermal_powers["advective"],
            "thermal_diffusive_power": thermal_powers["diffusion"],
            "thermal_gradient_power": thermal_powers["background_gradient"],
            "potential_advective_power": self.potential_moment(terms["advective"]),
            "potential_diffusive_power": self.potential_moment(terms["diffusion"]),
            "potential_buoyancy_exchange_residual": self.potential_moment(terms["advective"])+heat,
            "relative_divergence": float(np.linalg.norm(div)/(grad_norm or 1.)),
            "barotropic_vorticity_skewness": skew(zeta_bt),
            "vertical_vorticity_skewness": skew(zeta),
            "barotropic_vorticity_third_moment": float(np.mean(zeta_bt**3)),
            "barotropic_vorticity_third_moment_tendency": third_moment_rate(total_rhs),
            "barotropic_third_moment_full_advective": third_full_advective,
            "barotropic_third_moment_self_advective": third_bt_self,
            "barotropic_third_moment_convective": third_full_advective-third_bt_self,
            "barotropic_third_moment_viscous": third_moment_rate(terms["diffusion"]),
            "barotropic_third_moment_coriolis": third_moment_rate(terms["coriolis"]),
            "barotropic_third_moment_buoyancy": third_moment_rate(terms["buoyancy"]),
            "lowest_horizontal_shell_kinetic": lowest_energy,
            "lowest_barotropic_shell_kinetic": lowest_bt_energy,
            "outer_horizontal_band_kinetic": outer_energy,
            "outer_horizontal_band_energy_fraction": outer_energy/kin if kin else 0.,
            "barotropic_vorticity_rms": float(np.sqrt(np.mean(zeta_bt**2))),
            "velocity_rms": math.sqrt(2*kin),
            "maximum_velocity": float(np.max(np.sqrt(np.sum(fields[:3]**2, axis=0)))),
            "temperature_rms": math.sqrt(2*variance),
            "nusselt_volume": (1+heat/(self.kappa*self.G)) if self.kappa*self.G else None,
            "nusselt_bottom": 1-derivative_bottom/self.G if self.G else None,
            "nusselt_top": 1-derivative_top/self.G if self.G else None,
        }

    def planes(self, state, heights=(.2, .5, .8), nxy=64):
        """Exact modal reconstruction, not interpolation of a coarse picture."""
        heights = np.asarray(heights, dtype=float)
        assert np.all((0 <= heights) & (heights <= 1))
        assert nxy > 2*self.K
        phase = np.exp(1j*np.pi*self.iz[:, None]*heights[None, :])
        zeta_hat = 1j*(self.k[0]*state[1]-self.k[1]*state[0])
        coefficients = np.stack((zeta_hat, state[3]))
        horizontal = np.einsum("cxyz,zh->chxy", coefficients, phase)
        pad = np.zeros((2, len(heights), nxy, nxy), dtype=complex)
        for i, ix in enumerate(self.ix):
            if abs(ix) > self.K:
                continue
            for j, iy in enumerate(self.iy):
                if abs(iy) <= self.K:
                    pad[:, :, ix % nxy, iy % nxy] = horizontal[:, :, i, j]
        values = np.fft.ifftn(pad, axes=(-2, -1), norm="forward").real
        bt = np.zeros((nxy, nxy), dtype=complex)
        for i, ix in enumerate(self.ix):
            if abs(ix) > self.K:
                continue
            for j, iy in enumerate(self.iy):
                if abs(iy) <= self.K:
                    bt[ix % nxy, iy % nxy] = zeta_hat[i, j, 0]
        return {"heights": heights, "vorticity": values[0],
                "buoyancy_perturbation": values[1],
                "total_buoyancy": values[1]-self.G*heights[:, None, None],
                "barotropic_vorticity": np.fft.ifftn(bt, norm="forward").real,
                "x": np.arange(nxy)*(2*np.pi/self.kfund/nxy),
                "y": np.arange(nxy)*(2*np.pi/self.kfund/nxy)}

    def metadata(self):
        nh_nonzero = (2*self.K+1)**2-1
        return {
            "equations": "nonlinear Boussinesq-Coriolis with conductive buoyancy B0=-G*z",
            "K": self.K, "Nv": self.Nv, "kfund": self.kfund,
            "nu": self.nu, "kappa": self.kappa, "G": self.G, "f": self.f,
            "quadrature_grid": list(self.shape), "horizontal_period": 2*np.pi/self.kfund,
            "physical_depth": 1., "vertical_extension_period": 2.,
            "nonzero_horizontal_vertical_cutoff": self.Nv,
            "mean_velocity_vertical_cutoff": 2*self.Nv,
            "mean_temperature_vertical_cutoff": 2*self.Nv,
            "velocity_real_degrees_of_freedom": nh_nonzero+2*self.Nv*(nh_nonzero+2),
            "thermal_real_degrees_of_freedom": self.Nv*(nh_nonzero+2),
            "total_real_degrees_of_freedom": nh_nonzero+3*self.Nv*(nh_nonzero+2),
            "momentum_boundaries": "impermeable and stress-free",
            "buoyancy_boundaries": "B(0)=0, B(1)=-G; theta(0)=theta(1)=0",
            "nonlinear_projection": "orthogonal Galerkin with alias-free quadratic quadrature",
            "reference": "https://arxiv.org/abs/2505.07678",
            "cap": False, "anelastic": False, "beta": 0., "drag": 0.,
            "hyperviscosity": 0., "passive_tracer": False,
            "caveat": "reference-layer reduced model; not a quantitative reproduction of the Breeze atmosphere",
        }
