"""Bidirectional Maxwell scattering on a finite periodic transverse cell.

The input is one uniform-in-s layer with an arbitrary sampled two-dimensional
Ne field.  The exterior and its S1 ports come from a prepared S2 window.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import epsilon_0
from scipy.integrate import simpson
from scipy.linalg import expm
from scipy.sparse import csc_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply

from models.vector_pe_3d import (
    ETA0_OHM,
    _berreman_operator,
    cold_plasma_dielectric_tensor,
    dielectric_tensor_local,
)


def _periodic_k_matrix(count, period):
    modes = np.fft.fftfreq(count) * count
    grid = np.arange(count)
    fourier = np.exp(2j * np.pi * np.outer(grid, modes) / count) / np.sqrt(count)
    return fourier @ np.diag(2 * np.pi * modes / period) @ fourier.conj().T


def _port_matrix(window, side, direction):
    indices = window["indices"]
    periods = window["periods_m"]
    nu, nv = (len(np.unique(indices[:, j])) for j in range(2))
    u = np.arange(nu) * periods[0] / nu
    v = np.arange(nv) * periods[1] / nv
    uu, vv = np.meshgrid(u, v, indexing="ij")
    count = nu * nv
    matrix = np.empty((4 * count, 2 * count), dtype=complex)
    for i, slab in enumerate(window["slabs"]):
        ports = getattr(slab.left_ports if side == "left" else slab.right_ports, direction)
        ku, kv = window["tangent_k_rad_m"][i]
        harmonic = np.exp(1j * (ku * uu + kv * vv)).ravel() / np.sqrt(np.prod(periods))
        for j, port in enumerate(ports):
            for component, value in enumerate(port.state):
                matrix[component * count:(component + 1) * count, 2 * i + j] = value * harmonic
    return matrix


def _operator(ne_grid, medium, basis, ku, kv):
    epsilon = np.empty(ne_grid.shape + (3, 3), dtype=complex)
    for index in np.ndindex(ne_grid.shape):
        epsilon[index] = dielectric_tensor_local(cold_plasma_dielectric_tensor(
            medium.omega, ne_grid[index], medium.B_ecef_T, medium.collision_s,
            plasma_coefficient=medium.plasma_coefficient), basis)
    operator, electric_s = _berreman_operator(
        epsilon, medium.omega, ku, kv, np.ones(ne_grid.shape),
    )
    count = ne_grid.size
    scale = np.r_[np.ones(2 * count), np.full(2 * count, ETA0_OHM)]
    dense = operator.toarray()
    scaled = scale[:, None] * dense / scale[None, :]
    return epsilon, scaled, electric_s, scale


@dataclass(frozen=True)
class PeriodicSlabResult:
    S: np.ndarray
    S0: np.ndarray
    window: dict
    epsilon_suv: np.ndarray
    scaled_operator: np.ndarray
    electric_s_operators: tuple
    scale: np.ndarray
    left_forward_scaled: np.ndarray
    left_backward_scaled: np.ndarray
    matching_condition: float

    @property
    def count(self):
        return len(self.window["slabs"])

    @property
    def thickness_m(self):
        return self.window["slabs"][0].thickness_m

    def absorbed_power_W(self, incident, samples=33):
        """Integrate E† epsilon_loss E throughout the actual periodic layer."""
        incident = np.asarray(incident, dtype=complex)
        n = self.count
        if (incident.shape != (4 * n,) or not np.all(np.isfinite(incident))
                or not isinstance(samples, int) or samples < 3):
            raise ValueError("invalid incident coefficients or integration samples")
        outgoing = self.S @ incident
        initial = (self.left_forward_scaled @ incident[:2 * n]
                   + self.left_backward_scaled @ outgoing[:2 * n])
        matrix = csc_matrix(1j * self.scaled_operator)
        states = expm_multiply(matrix, initial, start=0, stop=self.thickness_m,
                               num=samples, endpoint=True, traceA=matrix.diagonal().sum())
        epsilon = self.epsilon_suv.reshape(n, 3, 3)
        loss = (epsilon - epsilon.swapaxes(-1, -2).conj()) / (2j)
        integrand = []
        for scaled_state in states:
            state = scaled_state / self.scale
            blocks = [state[j * n:(j + 1) * n] for j in range(4)]
            es = sum(op @ block for op, block in zip(self.electric_s_operators, blocks))
            electric = np.column_stack((es, blocks[0], blocks[1]))
            energy_density = np.einsum("ni,nij,nj->n", electric.conj(), loss, electric).real
            integrand.append(np.sum(energy_density) * np.prod(self.window["periods_m"]) / n)
        distance = np.linspace(0, self.thickness_m, samples)
        return float(self.window["left_medium"].omega * epsilon_0 / 2
                     * simpson(integrand, x=distance))


def solve_periodic_slab(window, layer, ne_grid_m3):
    """Return complete 4N-port S/S0 for the S2 window's periodic spectrum."""
    ne = np.asarray(ne_grid_m3, dtype=float)
    if (ne.ndim != 2 or ne.shape[0] != ne.shape[1] or ne.shape[0] < 3
            or ne.shape[0] % 2 != 1 or not np.all(np.isfinite(ne)) or np.any(ne < 0)):
        raise ValueError("Ne grid must be finite nonnegative on an odd square grid")
    size = ne.shape[0]
    count = size * size
    modes = np.arange(-(size // 2), size // 2 + 1)
    expected = np.array([(u, v) for u in modes for v in modes])
    if not np.array_equal(window["indices"], expected) or len(window["slabs"]) != count:
        raise ValueError("S2 window indices must match the full periodic grid")
    left = window["left_medium"]
    right = window["right_medium"]
    fields = ("frequency_hz", "time_s", "plasma_coefficient", "earth_radius_km",
              "coordinate_model", "near_cutoff_n2")
    if any(getattr(medium, field) != getattr(left, field)
           for medium in (layer, right) for field in fields):
        raise ValueError("periodic slab and exterior use inconsistent AH settings")
    basis = window["basis"]
    periods = window["periods_m"]
    ku = kron(csc_matrix(_periodic_k_matrix(size, periods[0])), eye(size, format="csc"), format="csc")
    kv = kron(eye(size, format="csc"), csc_matrix(_periodic_k_matrix(size, periods[1])), format="csc")
    epsilon, matrix, es_ops, scale = _operator(ne, layer, basis, ku, kv)
    _, matrix0, _, _ = _operator(np.full_like(ne, left.Ne_m3), left, basis, ku, kv)
    lf = scale[:, None] * _port_matrix(window, "left", "forward")
    lb = scale[:, None] * _port_matrix(window, "left", "backward")
    rf = scale[:, None] * _port_matrix(window, "right", "forward")
    rb = scale[:, None] * _port_matrix(window, "right", "backward")
    thickness = window["slabs"][0].thickness_m

    def scattering(operator):
        # shortcut: replace transfer matrices for thick or strongly evanescent layers.
        transfer = expm(1j * operator * thickness)
        match = np.column_stack((transfer @ lb, -rf))
        rhs = np.column_stack((-transfer @ lf, rb))
        return np.linalg.solve(match, rhs), float(np.linalg.cond(match))

    scattering_matrix, condition = scattering(matrix)
    background, _ = scattering(matrix0)
    if not np.all(np.isfinite(scattering_matrix)) or not np.isfinite(condition):
        raise ValueError("periodic slab matrix is non-finite or ill-conditioned")
    return PeriodicSlabResult(scattering_matrix, background, window, epsilon, matrix,
                              es_ops, scale, lf, lb, condition)


def apply_periodic_window(periodic, s2_window, s2_result):
    """Apply cross-bin S and S0 using the previously calculated S2 AH segments."""
    if s2_window is not periodic.window:
        raise ValueError("periodic and S2 results must share the same prepared window")
    if s2_result.get("window_context") is not s2_window:
        raise ValueError("S2 AH result was generated from a different prepared window")
    n = periodic.count
    source = np.asarray(s2_result["source"])
    receiver = np.asarray(s2_result["receiver"])
    before = np.asarray(s2_result["phase_before_rad"])
    after = np.asarray(s2_result["phase_after_rad"])
    if any(value.shape != (n, 2) for value in (source, receiver, before, after)):
        raise ValueError("S2 coefficients and AH phases must match periodic bins")
    rays = s2_result["rays"]
    if len(rays) != 4 * n:
        raise ValueError("S2 result does not contain both AH paths for each port")
    normal = s2_window["basis"][:, 0]
    for i, slab in enumerate(s2_window["slabs"]):
        groups = ((slab.left_ports.forward, s2_window["entrance_ecef_km"], 0),
                  (slab.right_ports.forward, s2_result["receiver_position_ecef_km"], 2))
        for ports, plane, offset in groups:
            for j, port in enumerate(ports):
                ray = rays[4 * i + offset + j]
                k_error = np.linalg.norm(ray["endpoint_k_ecef_rad_km"] - port.k_ecef_rad_km)
                k_error /= np.linalg.norm(port.k_ecef_rad_km)
                plane_error = abs(np.dot(ray["endpoint_ecef_km"] - plane, normal))
                if k_error > 1e-8 or plane_error > 1e-6:
                    raise ValueError("S2 AH paths do not match the periodic window ports and planes")
    at_window = (source * np.exp(1j * before)).reshape(2 * n)
    reflected = periodic.S[:2 * n, :2 * n] @ at_window
    forward = periodic.S[2 * n:, :2 * n] @ at_window
    background = periodic.S0[2 * n:, :2 * n] @ at_window
    transmitted = forward.reshape(n, 2) * np.exp(1j * after)
    transmitted0 = background.reshape(n, 2) * np.exp(1j * after)
    absorption = periodic.absorbed_power_W(np.r_[at_window, np.zeros(2 * n, complex)])
    rx = np.vdot(receiver, transmitted)
    rx0 = np.vdot(receiver, transmitted0)
    transmission_power = float(np.vdot(transmitted, transmitted).real)
    return {"incident_at_window": at_window, "reflected": reflected.reshape(n, 2),
            "transmitted": transmitted, "transmitted_background": transmitted0,
            "H_rx": rx, "H_rx_background": rx0,
            "incident_power_W": float(np.vdot(source, source).real),
            "reflected_power_W": float(np.vdot(reflected, reflected).real),
            "transmitted_power_W": transmission_power, "absorbed_power_W": absorption,
            "rx_power_W": float(abs(rx)**2), "uncaptured_power_W": transmission_power - float(abs(rx)**2)}
