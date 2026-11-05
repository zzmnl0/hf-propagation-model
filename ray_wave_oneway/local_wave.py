"""Power-normalized local Maxwell ports and a uniform planar slab.

The state uses (s, u, v), exp(i k.r - i omega t), SI distances and SI
wavevectors.  Each exterior port carries 1 W through a 1 m^2 reference area.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import c, epsilon_0, mu_0
from scipy.linalg import expm

from models.dispersion_ah import appleton_hartree_roots
from models.ne_provider_geocentric import COORDINATE_SYSTEM
from models.ray_tracer_ah_3d_ecef import COORDINATE_MODEL, WGS84_COORDINATE_MODEL
from models.vector_pe_3d import (
    ETA0_OHM,
    berreman_matrix_uniform,
    cold_plasma_dielectric_tensor,
    dielectric_tensor_local,
)


@dataclass(frozen=True)
class MediumSnapshot:
    frequency_hz: float
    time_s: float
    position_ecef_km: np.ndarray
    Ne_m3: float
    B_ecef_T: np.ndarray
    collision_s: float
    plasma_coefficient: float
    earth_radius_km: float
    coordinate_model: str
    near_cutoff_n2: float

    def __post_init__(self):
        position = np.asarray(self.position_ecef_km, dtype=float)
        field = np.asarray(self.B_ecef_T, dtype=float)
        scalars = np.asarray(
            [self.frequency_hz, self.time_s, self.Ne_m3,
             self.collision_s, self.plasma_coefficient, self.earth_radius_km,
             self.near_cutoff_n2], dtype=float
        )
        if (position.shape != (3,) or field.shape != (3,)
                or not np.all(np.isfinite(position))
                or not np.all(np.isfinite(field))
                or not np.all(np.isfinite(scalars))
                or self.frequency_hz <= 0 or self.Ne_m3 < 0
                or self.collision_s < 0 or self.plasma_coefficient <= 0
                or self.earth_radius_km <= 0 or self.near_cutoff_n2 < 0
                or self.coordinate_model not in (COORDINATE_MODEL, WGS84_COORDINATE_MODEL)):
            raise ValueError("invalid local medium snapshot")
        object.__setattr__(self, "position_ecef_km", position.copy())
        object.__setattr__(self, "B_ecef_T", field.copy())

    @classmethod
    def from_providers(cls, frequency_hz, time_s, position_ecef_km,
                       ne_provider, geomagnetic_provider, *, collision_s=0.0,
                       plasma_coefficient, earth_radius_km, coordinate_model,
                       near_cutoff_n2):
        position = np.asarray(position_ecef_km, dtype=float)
        expected_coordinates = {
            COORDINATE_MODEL: COORDINATE_SYSTEM,
            WGS84_COORDINATE_MODEL: "wgs84_geodetic",
        }
        if coordinate_model not in expected_coordinates:
            raise ValueError("unsupported AH coordinate model")
        metadata = ne_provider.get_metadata()
        if metadata.get("coordinate_system") != expected_coordinates[coordinate_model]:
            raise ValueError("Ne Provider coordinate model differs from AH model")
        if (coordinate_model == COORDINATE_MODEL
                and metadata.get("earth_radius_km") != earth_radius_km):
            raise ValueError("Ne Provider Earth radius differs from AH model")
        ne = ne_provider.get_Ne_ecef(position, earth_radius_km, t=time_s)
        field = geomagnetic_provider.get_B_ecef_T(position, t=time_s)
        return cls(frequency_hz, time_s, position, ne, field, collision_s,
                   plasma_coefficient, earth_radius_km, coordinate_model, near_cutoff_n2)

    @classmethod
    def from_ah_model(cls, model, position_ecef_km, time_s):
        return cls.from_providers(
            model.freq_MHz * 1e6, time_s, position_ecef_km,
            model.provider, model.geomagnetic_provider,
            collision_s=model.collision_s, plasma_coefficient=model.plasma_coefficient,
            earth_radius_km=model.earth_radius_km, coordinate_model=model.coordinate_model,
            near_cutoff_n2=model.near_cutoff_n2,
        )

    @property
    def omega(self):
        return 2 * np.pi * self.frequency_hz

    def epsilon(self):
        return cold_plasma_dielectric_tensor(
            self.omega, self.Ne_m3, self.B_ecef_T, self.collision_s,
            plasma_coefficient=self.plasma_coefficient,
        )


def _basis(basis_suv_ecef):
    basis = np.asarray(basis_suv_ecef, dtype=float)
    if (basis.shape != (3, 3) or not np.all(np.isfinite(basis))
            or not np.allclose(basis.T @ basis, np.eye(3), rtol=0, atol=1e-10)
            or abs(np.linalg.det(basis) - 1) > 1e-10):
        raise ValueError("basis_suv_ecef must be right-handed and orthonormal")
    return basis


def position_local_m(position_ecef_km, origin_ecef_km, basis_suv_ecef):
    return 1000.0 * _basis(basis_suv_ecef).T @ (
        np.asarray(position_ecef_km) - np.asarray(origin_ecef_km)
    )


def position_ecef_km(position_suv_m, origin_ecef_km, basis_suv_ecef):
    return np.asarray(origin_ecef_km) + (
        _basis(basis_suv_ecef) @ np.asarray(position_suv_m)
    ) / 1000.0


def wavevector_local_rad_m(k_ecef_rad_km, basis_suv_ecef):
    return _basis(basis_suv_ecef).T @ np.asarray(k_ecef_rad_km) / 1000.0


def wavevector_ecef_rad_km(k_suv_rad_m, basis_suv_ecef):
    return 1000.0 * _basis(basis_suv_ecef) @ np.asarray(k_suv_rad_m)


def _fields(state, epsilon_suv, omega, ku, kv):
    eu, ev, hu, hv = state
    es = (-epsilon_suv[0, 1] * eu - epsilon_suv[0, 2] * ev
          + (kv * hu - ku * hv) / (omega * epsilon_0)) / epsilon_suv[0, 0]
    hs = (-kv * eu + ku * ev) / (omega * mu_0)
    return np.array([es, eu, ev]), np.array([hs, hu, hv])


def _flux(state):
    eu, ev, hu, hv = state
    return 0.5 * np.real(eu * np.conj(hv) - ev * np.conj(hu))


@dataclass(frozen=True)
class PortMode:
    label: str
    direction: int
    k_suv_rad_m: np.ndarray
    electric_suv: np.ndarray
    magnetic_suv: np.ndarray
    basis_suv_ecef: np.ndarray

    @property
    def state(self):
        return np.r_[self.electric_suv[1:], self.magnetic_suv[1:]]

    @property
    def power_W(self):
        return float(_flux(self.state))

    @property
    def electric_ecef(self):
        return self.basis_suv_ecef @ self.electric_suv

    @property
    def magnetic_ecef(self):
        return self.basis_suv_ecef @ self.magnetic_suv

    @property
    def k_ecef_rad_km(self):
        return wavevector_ecef_rad_km(self.k_suv_rad_m, self.basis_suv_ecef)


@dataclass(frozen=True)
class ModalPorts:
    forward: tuple[PortMode, PortMode]
    backward: tuple[PortMode, PortMode]


def _make_mode(label, k_suv, state, epsilon, omega, ku, kv, basis):
    power = float(_flux(state))
    if not np.isfinite(power) or abs(power) < 1e-12:
        raise ValueError("exterior mode has zero normal energy flux")
    state = state / np.sqrt(abs(power))
    electric, magnetic = _fields(state, epsilon, omega, ku, kv)
    phase = np.exp(-1j * np.angle(electric[np.argmax(np.abs(electric))]))
    return PortMode(label, 1 if power > 0 else -1, np.asarray(k_suv),
                    electric * phase, magnetic * phase, basis.copy())


def modal_ports(snapshot, basis_suv_ecef, ku_rad_m, kv_rad_m):
    """Return two forward and backward lossless, propagating 1 W modes."""
    basis = _basis(basis_suv_ecef)
    ku, kv = float(ku_rad_m), float(kv_rad_m)
    if snapshot.collision_s != 0 or not np.isfinite(ku) or not np.isfinite(kv):
        raise ValueError("external ports must be lossless with finite tangential k")
    epsilon = dielectric_tensor_local(snapshot.epsilon(), basis)
    omega = snapshot.omega
    k0 = omega / c
    modes = []
    if np.linalg.norm(snapshot.B_ecef_T) <= 1e-18 or snapshot.Ne_m3 == 0:
        if epsilon[0, 0].real <= snapshot.near_cutoff_n2:
            raise ValueError("external port is near cutoff; AH reconnection is unsafe")
        q2 = (k0 * k0 * epsilon[0, 0]).real - ku * ku - kv * kv
        if q2 <= 0:
            raise ValueError("exterior modes must propagate")
        q = np.sqrt(q2)
        kt = np.hypot(ku, kv)
        te = np.array([0.0, -kv / kt, ku / kt]) if kt else np.array([0., 1., 0.])
        for direction in (1, -1):
            k = np.array([direction * q, ku, kv])
            for label, electric in (("TE", te), ("TM", np.cross(k, te) / np.linalg.norm(k))):
                magnetic = np.cross(k, electric) / (omega * mu_0)
                state = np.r_[electric[1:], magnetic[1:]].astype(complex)
                modes.append(_make_mode(label, k, state, epsilon, omega, ku, kv, basis))
    else:
        matrix = berreman_matrix_uniform(epsilon, omega, ku, kv)
        scale = np.diag([1., 1., ETA0_OHM, ETA0_OHM])
        values, vectors = np.linalg.eig(scale @ matrix @ np.linalg.inv(scale))
        field_local = basis.T @ snapshot.B_ecef_T
        gyro_mhz = 27_992_489_872.333 * np.linalg.norm(field_local) / 1e6
        for i, value in enumerate(values):
            if abs(value.imag) > 1e-8 * k0:
                raise ValueError("exterior modes must propagate")
            k = np.array([value.real, ku, kv])
            roots = appleton_hartree_roots(
                omega, snapshot.Ne_m3, k, field_local / np.linalg.norm(field_local),
                gyro_mhz, plasma_coefficient=snapshot.plasma_coefficient,
            )
            n2 = float(np.dot(k, k) / (k0 * k0))
            label = min(roots, key=lambda key: abs(roots[key] - n2))
            if not np.isfinite(roots[label].real) or roots[label].real <= snapshot.near_cutoff_n2:
                raise ValueError("external port is near cutoff; AH reconnection is unsafe")
            state = np.linalg.solve(scale, vectors[:, i])
            modes.append(_make_mode(label, k, state, epsilon, omega, ku, kv, basis))
    forward = sorted((mode for mode in modes if mode.direction == 1), key=lambda mode: mode.label)
    backward = sorted((mode for mode in modes if mode.direction == -1), key=lambda mode: mode.label)
    if len(forward) != 2 or len(backward) != 2:
        raise ValueError("exterior must have two propagating modes in each direction")
    if len({mode.label for mode in forward}) != 2 or len({mode.label for mode in backward}) != 2:
        raise ValueError("exterior modes are degenerate or cannot be identified")
    return ModalPorts(tuple(forward), tuple(backward))


@dataclass(frozen=True)
class SlabResult:
    S: np.ndarray
    S0: np.ndarray
    left_ports: ModalPorts
    right_ports: ModalPorts
    epsilon_slab_suv: np.ndarray
    slab_matrix: np.ndarray
    thickness_m: float
    omega: float
    ku: float
    kv: float

    def field_suv(self, distance_m, incident):
        incident = np.asarray(incident, dtype=complex)
        if (incident.shape != (4,) or not np.all(np.isfinite(incident))
                or not 0 <= distance_m <= self.thickness_m):
            raise ValueError("invalid incident coefficients or slab position")
        outgoing = self.S @ incident
        state0 = (np.column_stack([mode.state for mode in self.left_ports.forward]) @ incident[:2]
                  + np.column_stack([mode.state for mode in self.left_ports.backward]) @ outgoing[:2])
        state = expm(1j * self.slab_matrix * distance_m) @ state0
        return _fields(state, self.epsilon_slab_suv, self.omega, self.ku, self.kv)

    def absorbed_power_W(self, incident, quadrature_order=64):
        nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
        loss = (self.epsilon_slab_suv - self.epsilon_slab_suv.conj().T) / (2j)
        integral = 0.0
        for node, weight in zip(nodes, weights):
            electric, _ = self.field_suv(self.thickness_m * (node + 1) / 2, incident)
            integral += weight * np.real(np.vdot(electric, loss @ electric))
        return float(self.omega * epsilon_0 * self.thickness_m * integral / 4)


def solve_slab(left, slab, right, background_slab, basis_suv_ecef,
               ku_rad_m, kv_rad_m, thickness_m):
    """Solve complex 4-port slab S and same-port background S0.

    Inputs are [left forward, right backward], outputs [left backward,
    right forward]. The common port and phase planes are at s=0 and s=d.
    """
    basis = _basis(basis_suv_ecef)
    media = (left, slab, right, background_slab)
    if any((m.frequency_hz != left.frequency_hz or m.time_s != left.time_s
            or m.plasma_coefficient != left.plasma_coefficient
            or m.earth_radius_km != left.earth_radius_km
            or m.near_cutoff_n2 != left.near_cutoff_n2
            or m.coordinate_model != left.coordinate_model) for m in media):
        raise ValueError("slabs and ports must share frequency, time, dispersion and Earth model")
    if not np.isfinite(thickness_m) or thickness_m <= 0:
        raise ValueError("thickness_m must be finite and positive")
    lp = modal_ports(left, basis, ku_rad_m, kv_rad_m)
    rp = modal_ports(right, basis, ku_rad_m, kv_rad_m)
    lf, lb = (np.column_stack([mode.state for mode in group]) for group in (lp.forward, lp.backward))
    rf, rb = (np.column_stack([mode.state for mode in group]) for group in (rp.forward, rp.backward))

    def solve_layer(medium):
        epsilon = dielectric_tensor_local(medium.epsilon(), basis)
        matrix = berreman_matrix_uniform(epsilon, medium.omega, ku_rad_m, kv_rad_m)
        # shortcut: use a stable scattering formulation for thick or strongly evanescent layers.
        transfer = expm(1j * matrix * thickness_m)
        scattering = np.linalg.solve(
            np.column_stack((transfer @ lb, -rf)),
            np.column_stack((-transfer @ lf, rb)),
        )
        return scattering, matrix, epsilon

    scattering, matrix, epsilon = solve_layer(slab)
    background, _, _ = solve_layer(background_slab)
    return SlabResult(scattering, background, lp, rp, epsilon, matrix,
                      float(thickness_m), left.omega, float(ku_rad_m), float(kv_rad_m))
