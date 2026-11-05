"""CPU reference solver for a two-transverse-dimensional vector PE.

The solver is deliberately independent of the production propagation pipeline.
It marches the transverse Maxwell state ``[E_u, E_v, H_u, H_v]`` in a
ray-centred Bishop frame using a [1/1] Pade (Crank-Nicolson) step.  The model is
one-way in the sense that it is an initial-value march along ``s``; it is not a
replacement for FDTD when appreciable backward power is present.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from typing import Callable

import numpy as np
from scipy.constants import c, epsilon_0, mu_0
from scipy.interpolate import CubicHermiteSpline
from scipy.sparse import bmat, csc_matrix, diags, eye, kron
from scipy.sparse.linalg import splu

from config import K_FP


ETA0_OHM = float(np.sqrt(mu_0 / epsilon_0))


def _unit(vector, *, name: str) -> np.ndarray:
    value = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(value))
    if value.shape != (3,) or not np.all(np.isfinite(value)) or norm <= 1e-15:
        raise ValueError(f"{name} must be a finite non-zero length-3 vector")
    return value / norm


def _stable_hash(value) -> str:
    def convert(item):
        if isinstance(item, np.ndarray):
            array = np.ascontiguousarray(item)
            return {
                "shape": list(array.shape),
                "dtype": str(array.dtype),
                "sha256": hashlib.sha256(array.tobytes()).hexdigest(),
            }
        if isinstance(item, (np.floating, np.integer)):
            return item.item()
        if isinstance(item, complex):
            return [float(item.real), float(item.imag)]
        if isinstance(item, dict):
            return {str(k): convert(v) for k, v in sorted(item.items())}
        if isinstance(item, (list, tuple)):
            return [convert(v) for v in item]
        return item

    payload = json.dumps(
        convert(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def cold_plasma_dielectric_tensor(
    omega_rad_s: float,
    Ne_m3: float,
    B_ecef_T,
    collision_s: float = 0.0,
    *,
    plasma_coefficient: float | None = None,
) -> np.ndarray:
    """Return the complex cold-plasma dielectric tensor in ECEF coordinates.

    The time convention is exp(-i*omega*t).  With this convention a passive
    collisional plasma has a positive-semidefinite anti-Hermitian loss part.
    """
    omega = float(omega_rad_s)
    Ne = float(Ne_m3)
    nu = float(collision_s)
    field = np.asarray(B_ecef_T, dtype=float)
    if (
        not np.isfinite(omega)
        or omega <= 0.0
        or not np.isfinite(Ne)
        or Ne < 0.0
        or not np.isfinite(nu)
        or nu < 0.0
        or field.shape != (3,)
        or not np.all(np.isfinite(field))
    ):
        raise ValueError("invalid cold-plasma tensor input")

    coefficient = K_FP**2 if plasma_coefficient is None else float(plasma_coefficient)
    if not np.isfinite(coefficient) or coefficient <= 0.0:
        raise ValueError("plasma_coefficient must be finite and positive")
    frequency_hz = omega / (2.0 * np.pi)
    X = coefficient * Ne / (frequency_hz * frequency_hz)
    U = 1.0 + 1j * nu / omega
    magnitude = float(np.linalg.norm(field))
    if magnitude <= 1e-18:
        scalar = 1.0 - X / U
        return np.eye(3, dtype=complex) * scalar

    # Electron gyrofrequency magnitude in Hz/T, consistent with AH3D.
    gyro_hz = 27_992_489_872.333 * magnitude
    Y = gyro_hz / frequency_hz
    R = 1.0 - X / (U - Y)
    L = 1.0 - X / (U + Y)
    S = 0.5 * (R + L)
    D = 0.5 * (R - L)
    P = 1.0 - X / U
    tensor_b = np.array(
        [[S, -1j * D, 0.0], [1j * D, S, 0.0], [0.0, 0.0, P]],
        dtype=complex,
    )

    b_hat = field / magnitude
    reference = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(reference, b_hat))) > 0.9:
        reference = np.array([1.0, 0.0, 0.0])
    x_hat = _unit(np.cross(reference, b_hat), name="field transverse axis")
    y_hat = _unit(np.cross(b_hat, x_hat), name="field transverse axis")
    rotation = np.column_stack([x_hat, y_hat, b_hat])
    return rotation @ tensor_b @ rotation.T


def dielectric_tensor_local(epsilon_ecef, basis_suv) -> np.ndarray:
    """Rotate an ECEF dielectric tensor into the local (s,u,v) basis."""
    epsilon = np.asarray(epsilon_ecef, dtype=complex)
    basis = np.asarray(basis_suv, dtype=float)
    if epsilon.shape != (3, 3) or basis.shape != (3, 3):
        raise ValueError("epsilon_ecef and basis_suv must be 3x3")
    return basis.T @ epsilon @ basis


@dataclass(frozen=True)
class RayCenteredFrame3D:
    s_m: np.ndarray
    centre_ecef_km: np.ndarray
    basis_suv_ecef: np.ndarray
    curvature_u_m_inv: np.ndarray
    curvature_v_m_inv: np.ndarray
    jacobian: np.ndarray
    u_m: np.ndarray
    v_m: np.ndarray
    frame_hash: str

    @classmethod
    def from_ah3d(cls, ray: dict, s_m, u_m, v_m) -> "RayCenteredFrame3D":
        points = np.asarray(ray.get("trajectory_ecef_km", []), dtype=float)
        if points.ndim != 2 or points.shape[0] < 2 or points.shape[1] < 3:
            raise ValueError("AH3D ray requires trajectory_ecef_km")
        points = points[:, :3]
        path_km = np.asarray(ray.get("trajectory_geom_path_km", []), dtype=float)
        if path_km.shape != (points.shape[0],):
            path_km = np.concatenate(
                [[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
            )
        keep = np.concatenate([[True], np.diff(path_km) > 1e-10])
        points = points[keep]
        path_km = path_km[keep]
        if points.shape[0] < 2 or not np.all(np.diff(path_km) > 0.0):
            raise ValueError("AH3D trajectory path must be strictly increasing")

        stored = np.asarray(
            ray.get("trajectory_group_tangent_ecef", []), dtype=float
        )
        if stored.shape == (keep.size, 3):
            stored = stored[keep]
        else:
            stored = np.gradient(points, path_km, axis=0)
        tangents = np.vstack([_unit(row, name="AH3D tangent") for row in stored])

        s_query_m = np.asarray(s_m, dtype=float)
        u_axis = np.asarray(u_m, dtype=float)
        v_axis = np.asarray(v_m, dtype=float)
        if (
            s_query_m.ndim != 1
            or s_query_m.size < 2
            or not np.all(np.diff(s_query_m) > 0.0)
            or s_query_m[0] < path_km[0] * 1000.0 - 1e-8
            or s_query_m[-1] > path_km[-1] * 1000.0 + 1e-8
        ):
            raise ValueError("s_m must be an increasing subset of the AH3D path")
        if (
            u_axis.ndim != 1
            or v_axis.ndim != 1
            or u_axis.size < 3
            or v_axis.size < 3
            or not np.all(np.diff(u_axis) > 0.0)
            or not np.all(np.diff(v_axis) > 0.0)
        ):
            raise ValueError("u_m and v_m must be increasing one-dimensional axes")

        spline = CubicHermiteSpline(path_km, points, tangents, axis=0)
        s_query_km = s_query_m / 1000.0
        centre = np.asarray(spline(s_query_km), dtype=float)
        tangent = np.asarray(spline(s_query_km, 1), dtype=float)
        tangent = np.vstack([_unit(row, name="interpolated tangent") for row in tangent])

        radial = _unit(centre[0], name="initial radial direction")
        u0 = radial - np.dot(radial, tangent[0]) * tangent[0]
        if np.linalg.norm(u0) < 1e-12:
            reference = np.array([0.0, 0.0, 1.0])
            if abs(float(np.dot(reference, tangent[0]))) > 0.9:
                reference = np.array([1.0, 0.0, 0.0])
            u0 = reference - np.dot(reference, tangent[0]) * tangent[0]
        u_basis = np.empty_like(tangent)
        v_basis = np.empty_like(tangent)
        u_basis[0] = _unit(u0, name="initial Bishop u")
        v_basis[0] = _unit(np.cross(tangent[0], u_basis[0]), name="initial Bishop v")

        for index in range(1, tangent.shape[0]):
            previous = tangent[index - 1]
            current = tangent[index]
            axis = np.cross(previous, current)
            sine = float(np.linalg.norm(axis))
            cosine = float(np.clip(np.dot(previous, current), -1.0, 1.0))
            if sine < 1e-14:
                candidate = u_basis[index - 1]
            else:
                axis /= sine
                angle = float(np.arctan2(sine, cosine))
                old = u_basis[index - 1]
                candidate = (
                    old * np.cos(angle)
                    + np.cross(axis, old) * np.sin(angle)
                    + axis * np.dot(axis, old) * (1.0 - np.cos(angle))
                )
            candidate -= np.dot(candidate, current) * current
            u_basis[index] = _unit(candidate, name="transported Bishop u")
            v_basis[index] = _unit(
                np.cross(current, u_basis[index]), name="transported Bishop v"
            )

        dt_ds = np.gradient(tangent, s_query_m, axis=0, edge_order=2)
        curvature_u = np.einsum("ij,ij->i", dt_ds, u_basis)
        curvature_v = np.einsum("ij,ij->i", dt_ds, v_basis)
        uu, vv = np.meshgrid(u_axis, v_axis, indexing="ij")
        jacobian = (
            1.0
            - curvature_u[:, None, None] * uu[None, :, :]
            - curvature_v[:, None, None] * vv[None, :, :]
        )
        if not np.all(np.isfinite(jacobian)) or float(np.min(jacobian)) < 0.8:
            raise ValueError("coordinate_fold_risk: ray-centred Jacobian below 0.8")
        basis = np.stack([tangent, u_basis, v_basis], axis=2)
        payload = {
            "s_m": s_query_m,
            "centre_ecef_km": centre,
            "basis_suv_ecef": basis,
            "curvature_u_m_inv": curvature_u,
            "curvature_v_m_inv": curvature_v,
            "u_m": u_axis,
            "v_m": v_axis,
        }
        return cls(
            s_m=s_query_m,
            centre_ecef_km=centre,
            basis_suv_ecef=basis,
            curvature_u_m_inv=curvature_u,
            curvature_v_m_inv=curvature_v,
            jacobian=jacobian,
            u_m=u_axis,
            v_m=v_axis,
            frame_hash=_stable_hash(payload),
        )

    def ecef_grid_km(self, plane_index: int) -> np.ndarray:
        basis = self.basis_suv_ecef[int(plane_index)]
        uu, vv = np.meshgrid(self.u_m, self.v_m, indexing="ij")
        offsets_km = (
            uu[..., None] * basis[:, 1] + vv[..., None] * basis[:, 2]
        ) / 1000.0
        return self.centre_ecef_km[int(plane_index)] + offsets_km


@dataclass(frozen=True)
class VectorPE3DConfig:
    frequency_Hz: float
    pml_cells: int = 16
    pml_order: int = 3
    pml_sigma_max: float = 4.0
    pade_order: tuple[int, int] = (1, 1)
    max_angle_deg: float = 30.0
    max_out_of_cone_power: float = 0.01
    max_backward_power: float = 0.01
    linear_residual_tol: float = 1.0e-10
    min_forward_backward_gap_k0: float = 1.0e-3
    store_plane_indices: tuple[int, ...] = ()
    semantics: str = "research_reference_only"

    def validate(self, frame: RayCenteredFrame3D) -> None:
        if (
            not np.isfinite(self.frequency_Hz)
            or self.frequency_Hz <= 0.0
            or self.pade_order != (1, 1)
            or self.pml_cells < 2
            or self.pml_order < 1
            or not np.isfinite(self.pml_sigma_max)
            or self.pml_sigma_max <= 0.0
            or not (0.0 < self.max_angle_deg < 90.0)
            or not (0.0 <= self.max_out_of_cone_power < 1.0)
            or not (0.0 <= self.max_backward_power < 1.0)
            or self.linear_residual_tol <= 0.0
        ):
            raise ValueError("invalid VectorPE3DConfig")
        if 2 * self.pml_cells >= min(frame.u_m.size, frame.v_m.size):
            raise ValueError("PML occupies the complete transverse domain")
        wavelength = c / self.frequency_Hz
        du = float(np.min(np.diff(frame.u_m)))
        dv = float(np.min(np.diff(frame.v_m)))
        if self.pml_cells * min(du, dv) < 2.0 * wavelength:
            raise ValueError("PML physical thickness must be at least two wavelengths")


@dataclass
class VectorPE3DResult:
    final_field_transverse: np.ndarray
    final_E_s: np.ndarray
    final_H_s: np.ndarray
    poynting_s_W_m2: np.ndarray
    flux_history_W: np.ndarray
    linear_residual_max: float
    angular_out_of_cone_power_fraction: float
    backward_power_fraction: float
    collision_absorbed_W: float
    pml_or_numerical_loss_W: float
    energy_residual_fraction: float
    termination_status: str
    diagnostics: dict
    stored_fields: dict[int, np.ndarray]


def _first_derivative_matrix(axis: np.ndarray, stretch: np.ndarray) -> csc_matrix:
    size = axis.size
    spacing = np.diff(axis)
    if not np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-12):
        raise ValueError("B1 reference derivative requires a uniform transverse grid")
    h = float(spacing[0])
    matrix = diags(
        [-np.ones(size - 1), np.ones(size - 1)],
        [-1, 1],
        shape=(size, size),
        dtype=complex,
    ) / (2.0 * h)
    return (diags(1.0 / stretch, dtype=complex) @ matrix).tocsc()


def _pml_stretch(size: int, cells: int, order: int, sigma_max: float) -> np.ndarray:
    stretch = np.ones(size, dtype=complex)
    for index in range(cells):
        fraction = (cells - index) / float(cells)
        value = sigma_max * fraction**order
        stretch[index] = 1.0 + 1j * value
        stretch[-1 - index] = 1.0 + 1j * value
    return stretch


def _berreman_operator(
    epsilon_suv: np.ndarray,
    omega: float,
    K_u: csc_matrix,
    K_v: csc_matrix,
    jacobian: np.ndarray,
) -> tuple[csc_matrix, tuple[csc_matrix, ...]]:
    """Build the collocated four-component Maxwell spatial operator."""
    shape = epsilon_suv.shape[:2]
    count = int(np.prod(shape))
    epsilon = epsilon_suv.reshape(count, 3, 3)
    eps_ss = epsilon[:, 0, 0]
    if np.any(~np.isfinite(eps_ss)) or np.any(np.abs(eps_ss) < 1e-8):
        raise ValueError("longitudinal dielectric component is singular")

    def diagonal(values):
        return diags(np.asarray(values), dtype=complex).tocsc()

    inv_eps_ss = diagonal(1.0 / eps_ss)
    es_eu = diagonal(-epsilon[:, 0, 1] / eps_ss)
    es_ev = diagonal(-epsilon[:, 0, 2] / eps_ss)
    es_hu = inv_eps_ss @ K_v / (omega * epsilon_0)
    es_hv = -inv_eps_ss @ K_u / (omega * epsilon_0)
    hs_eu = -K_v / (omega * mu_0)
    hs_ev = K_u / (omega * mu_0)
    zero = csc_matrix((count, count), dtype=complex)
    identity = eye(count, dtype=complex, format="csc")

    row_eu = [
        K_u @ es_eu,
        K_u @ es_ev,
        K_u @ es_hu,
        K_u @ es_hv + omega * mu_0 * identity,
    ]
    row_ev = [
        K_v @ es_eu,
        K_v @ es_ev,
        K_v @ es_hu - omega * mu_0 * identity,
        K_v @ es_hv,
    ]

    eps_vu = diagonal(epsilon[:, 2, 1])
    eps_vv = diagonal(epsilon[:, 2, 2])
    eps_vs = diagonal(epsilon[:, 2, 0])
    eps_uu = diagonal(epsilon[:, 1, 1])
    eps_uv = diagonal(epsilon[:, 1, 2])
    eps_us = diagonal(epsilon[:, 1, 0])
    row_hu = [
        K_u @ hs_eu - omega * epsilon_0 * (eps_vu + eps_vs @ es_eu),
        K_u @ hs_ev - omega * epsilon_0 * (eps_vv + eps_vs @ es_ev),
        -omega * epsilon_0 * eps_vs @ es_hu,
        -omega * epsilon_0 * eps_vs @ es_hv,
    ]
    row_hv = [
        K_v @ hs_eu + omega * epsilon_0 * (eps_uu + eps_us @ es_eu),
        K_v @ hs_ev + omega * epsilon_0 * (eps_uv + eps_us @ es_ev),
        omega * epsilon_0 * eps_us @ es_hu,
        omega * epsilon_0 * eps_us @ es_hv,
    ]
    operator = bmat(
        [row_eu, row_ev, row_hu, row_hv], format="csc", dtype=complex
    )
    metric = diags(np.asarray(jacobian).reshape(count), dtype=complex)
    operator = bmat(
        [
            [metric if i == j else zero for j in range(4)]
            for i in range(4)
        ],
        format="csc",
    ) @ operator
    return operator, (es_eu, es_ev, es_hu, es_hv)


def berreman_matrix_uniform(
    epsilon_suv: np.ndarray,
    omega_rad_s: float,
    k_u_rad_m: float = 0.0,
    k_v_rad_m: float = 0.0,
) -> np.ndarray:
    """Dense uniform-medium matrix used for independent modal diagnostics."""
    ku = csc_matrix([[complex(k_u_rad_m)]])
    kv = csc_matrix([[complex(k_v_rad_m)]])
    operator, _ = _berreman_operator(
        np.asarray(epsilon_suv, dtype=complex).reshape(1, 1, 3, 3),
        float(omega_rad_s),
        ku,
        kv,
        np.ones((1, 1)),
    )
    return operator.toarray()


def _plane_flux(field: np.ndarray, jacobian, du: float, dv: float) -> tuple[np.ndarray, float]:
    eu, ev, hu, hv = field
    poynting_s = 0.5 * np.real(eu * np.conj(hv) - ev * np.conj(hu))
    flux = float(np.sum(poynting_s * jacobian) * du * dv)
    return poynting_s, flux


def _boundary_mask(shape: tuple[int, int]) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    mask[[0, -1], :] = True
    mask[:, [0, -1]] = True
    return mask


class VectorPE3D:
    """Independent CPU reference implementation of the B1 vector PE core."""

    def __init__(self, config: VectorPE3DConfig):
        self.config = config

    def _sample_epsilon(
        self,
        centre_ecef_km,
        basis_suv_ecef,
        u_m,
        v_m,
        ne_provider,
        geomagnetic_provider,
        earth_radius_km: float,
        collision_sampler: Callable[[np.ndarray], float],
    ) -> np.ndarray:
        basis = np.asarray(basis_suv_ecef, dtype=float)
        uu, vv = np.meshgrid(u_m, v_m, indexing="ij")
        grid = (
            np.asarray(centre_ecef_km, dtype=float)
            + (
                uu[..., None] * basis[:, 1]
                + vv[..., None] * basis[:, 2]
            )
            / 1000.0
        )
        output = np.empty(grid.shape[:2] + (3, 3), dtype=complex)
        omega = 2.0 * np.pi * self.config.frequency_Hz
        for index in np.ndindex(grid.shape[:2]):
            point = grid[index]
            Ne = float(ne_provider.get_Ne_ecef(point, earth_radius_km, t=0.0))
            B = np.asarray(
                geomagnetic_provider.get_B_ecef_T(point, t=0.0), dtype=float
            )
            collision = float(collision_sampler(point))
            epsilon_ecef = cold_plasma_dielectric_tensor(
                omega, Ne, B, collision_s=collision
            )
            output[index] = dielectric_tensor_local(epsilon_ecef, basis)
        return output

    def solve(
        self,
        frame: RayCenteredFrame3D,
        ne_provider,
        geomagnetic_provider,
        initial_field_transverse,
        *,
        earth_radius_km: float,
        collision_sampler: Callable[[np.ndarray], float] | None = None,
    ) -> VectorPE3DResult:
        self.config.validate(frame)
        collision_fn = collision_sampler or (lambda _point: 0.0)
        field = np.asarray(initial_field_transverse, dtype=complex)
        expected = (4, frame.u_m.size, frame.v_m.size)
        if field.shape != expected or not np.all(np.isfinite(field)):
            raise ValueError(f"initial_field_transverse must have shape {expected}")
        field = field.copy()
        du = float(np.diff(frame.u_m)[0])
        dv = float(np.diff(frame.v_m)[0])
        count = frame.u_m.size * frame.v_m.size
        stretch_u = _pml_stretch(
            frame.u_m.size,
            self.config.pml_cells,
            self.config.pml_order,
            self.config.pml_sigma_max,
        )
        stretch_v = _pml_stretch(
            frame.v_m.size,
            self.config.pml_cells,
            self.config.pml_order,
            self.config.pml_sigma_max,
        )
        derivative_u = _first_derivative_matrix(frame.u_m, stretch_u)
        derivative_v = _first_derivative_matrix(frame.v_m, stretch_v)
        K_u = (-1j * kron(derivative_u, eye(frame.v_m.size))).tocsc()
        K_v = (-1j * kron(eye(frame.u_m.size), derivative_v)).tocsc()
        identity = eye(4 * count, dtype=complex, format="csc")
        boundary = _boundary_mask(field.shape[1:])
        flux_history = []
        stored = {}
        residual_max = 0.0
        collision_absorbed = 0.0
        modal_gap_min = np.inf
        omega = 2.0 * np.pi * self.config.frequency_Hz
        carrier_k0 = omega / c
        last_es_ops = None
        last_epsilon = None
        cached_key = None
        cached_factor = None
        cached_right = None

        _, initial_flux = _plane_flux(
            field, frame.jacobian[0], du, dv
        )
        flux_history.append(initial_flux)
        if initial_flux <= 0.0:
            raise ValueError("initial field must carry positive s-directed power")

        for step in range(frame.s_m.size - 1):
            ds = float(frame.s_m[step + 1] - frame.s_m[step])
            midpoint_centre = 0.5 * (
                frame.centre_ecef_km[step] + frame.centre_ecef_km[step + 1]
            )
            tangent = _unit(
                frame.basis_suv_ecef[step, :, 0]
                + frame.basis_suv_ecef[step + 1, :, 0],
                name="midpoint tangent",
            )
            u_mid = (
                frame.basis_suv_ecef[step, :, 1]
                + frame.basis_suv_ecef[step + 1, :, 1]
            )
            u_mid -= np.dot(u_mid, tangent) * tangent
            u_mid = _unit(u_mid, name="midpoint Bishop u")
            v_mid = _unit(np.cross(tangent, u_mid), name="midpoint Bishop v")
            midpoint_basis = np.column_stack([tangent, u_mid, v_mid])
            epsilon = self._sample_epsilon(
                midpoint_centre,
                midpoint_basis,
                frame.u_m,
                frame.v_m,
                ne_provider,
                geomagnetic_provider,
                float(earth_radius_km),
                collision_fn,
            )
            jacobian = 0.5 * (
                frame.jacobian[step] + frame.jacobian[step + 1]
            )
            operator, es_ops = _berreman_operator(
                epsilon, omega, K_u, K_v, jacobian
            )
            centre_epsilon = epsilon[
                frame.u_m.size // 2, frame.v_m.size // 2
            ]
            modal_values = np.linalg.eigvals(
                berreman_matrix_uniform(centre_epsilon, omega)
            )
            positive = modal_values[modal_values.real > 0.0]
            negative = modal_values[modal_values.real < 0.0]
            if positive.size and negative.size:
                gap = min(
                    abs(p - n) for p in positive for n in negative
                ) / (omega / c)
                modal_gap_min = min(modal_gap_min, float(gap))
            if modal_gap_min < self.config.min_forward_backward_gap_k0:
                raise RuntimeError("forward/backward modal gap is unresolved")

            envelope_operator = operator - carrier_k0 * identity
            matrix_key = _stable_hash(
                {
                    "epsilon": epsilon,
                    "jacobian": jacobian,
                    "ds": ds,
                }
            )
            if matrix_key != cached_key:
                left = identity - 0.5j * ds * envelope_operator
                cached_right = identity + 0.5j * ds * envelope_operator
                cached_factor = splu(left)
                cached_key = matrix_key
            else:
                left = identity - 0.5j * ds * envelope_operator
            vector = field.reshape(4 * count)
            rhs = cached_right @ vector
            updated = cached_factor.solve(rhs)
            denominator = max(float(np.linalg.norm(rhs)), 1e-30)
            residual = float(np.linalg.norm(left @ updated - rhs) / denominator)
            residual_max = max(residual_max, residual)
            if not np.all(np.isfinite(updated)) or residual > self.config.linear_residual_tol:
                raise RuntimeError("vector PE sparse solve failed its residual gate")
            field = updated.reshape(expected)

            # Parallel-transport the transverse components to the next Bishop basis.
            old_basis = frame.basis_suv_ecef[step]
            new_basis = frame.basis_suv_ecef[step + 1]
            transform = new_basis[:, 1:3].T @ old_basis[:, 1:3]
            field[:2] = np.einsum("ab,bij->aij", transform, field[:2])
            field[2:] = np.einsum("ab,bij->aij", transform, field[2:])
            field[:, boundary] = 0.0

            vector_mid = 0.5 * (
                vector + field.reshape(4 * count)
            )
            es = sum(block @ vector_mid[i * count:(i + 1) * count]
                     for i, block in enumerate(es_ops))
            electric = np.stack(
                [
                    es.reshape(field.shape[1:]),
                    field[0],
                    field[1],
                ],
                axis=-1,
            ).reshape(count, 3)
            eps_flat = epsilon.reshape(count, 3, 3)
            absorption_density = np.empty(count, dtype=float)
            for index in range(count):
                loss_tensor = (
                    eps_flat[index] - eps_flat[index].conj().T
                ) / (2j)
                absorption_density[index] = max(
                    float(
                        0.5
                        * omega
                        * epsilon_0
                        * np.real(
                            np.vdot(
                                electric[index],
                                loss_tensor @ electric[index],
                            )
                        )
                    ),
                    0.0,
                )
            collision_absorbed += float(
                np.sum(absorption_density.reshape(field.shape[1:]) * jacobian)
                * du
                * dv
                * ds
            )
            _, flux = _plane_flux(
                field, frame.jacobian[step + 1], du, dv
            )
            flux_history.append(flux)
            if step + 1 in self.config.store_plane_indices:
                stored[step + 1] = field.copy()
            last_es_ops = es_ops
            last_epsilon = epsilon

        vector = field.reshape(4 * count)
        final_es = sum(
            block @ vector[i * count:(i + 1) * count]
            for i, block in enumerate(last_es_ops)
        ).reshape(field.shape[1:])
        final_hs = (
            K_u @ vector[count:2 * count] - K_v @ vector[:count]
        ) / (omega * mu_0)
        final_hs = final_hs.reshape(field.shape[1:])
        poynting, final_flux = _plane_flux(
            field, frame.jacobian[-1], du, dv
        )
        negative = float(
            np.sum(np.maximum(-poynting, 0.0) * frame.jacobian[-1]) * du * dv
        )
        backward_fraction = negative / max(abs(final_flux) + negative, 1e-30)

        spectrum_u = np.fft.fftshift(np.fft.fft2(field[0]))
        spectrum_v = np.fft.fftshift(np.fft.fft2(field[1]))
        spectrum_power = np.abs(spectrum_u) ** 2 + np.abs(spectrum_v) ** 2
        ku = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(frame.u_m.size, du))
        kv = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(frame.v_m.size, dv))
        kuu, kvv = np.meshgrid(ku, kv, indexing="ij")
        k0 = omega / c
        angle = np.degrees(np.arcsin(np.clip(np.hypot(kuu, kvv) / k0, 0.0, 1.0)))
        outside = angle > self.config.max_angle_deg
        outside_fraction = float(
            np.sum(spectrum_power[outside]) / max(np.sum(spectrum_power), 1e-30)
        )
        pml_or_numerical = max(
            initial_flux - final_flux - collision_absorbed, 0.0
        )
        residual_energy = abs(
            initial_flux - final_flux - collision_absorbed - pml_or_numerical
        ) / max(abs(initial_flux), 1e-30)
        status = "ok"
        if outside_fraction > self.config.max_out_of_cone_power:
            status = "wide_angle_validity_exceeded"
        if backward_fraction > self.config.max_backward_power:
            status = "backward_power_validity_exceeded"

        diagnostics = {
            "solver": "vector_pe_3d_cpu_reference",
            "solver_dimensionality": "full_3d_two_transverse_vector_pe",
            "propagation_direction": "one_way_forward",
            "validity": self.config.semantics,
            "calibration_status": "reference_uncalibrated",
            "time_convention": "exp(-i*omega*t)",
            "field_state_order": ["E_u", "E_v", "H_u", "H_v"],
            "frequency_Hz": float(self.config.frequency_Hz),
            "wavelength_m": float(c / self.config.frequency_Hz),
            "carrier_wavenumber_rad_m": float(carrier_k0),
            "field_representation": "complex_envelope_times_exp(i*k0*s)",
            "frame_hash": frame.frame_hash,
            "config_hash": _stable_hash(asdict(self.config)),
            "grid_shape": [int(frame.s_m.size), int(frame.u_m.size), int(frame.v_m.size)],
            "pml_cells": int(self.config.pml_cells),
            "pml_order": int(self.config.pml_order),
            "pml_sigma_max": float(self.config.pml_sigma_max),
            "min_jacobian": float(np.min(frame.jacobian)),
            "max_angle_deg": float(self.config.max_angle_deg),
            "forward_backward_gap_k0_min": float(modal_gap_min),
            "relative_flux_change": float(
                (final_flux - initial_flux) / initial_flux
            ),
            "last_epsilon_hash": _stable_hash(last_epsilon),
        }
        return VectorPE3DResult(
            final_field_transverse=field,
            final_E_s=final_es,
            final_H_s=final_hs,
            poynting_s_W_m2=poynting,
            flux_history_W=np.asarray(flux_history),
            linear_residual_max=float(residual_max),
            angular_out_of_cone_power_fraction=outside_fraction,
            backward_power_fraction=backward_fraction,
            collision_absorbed_W=float(collision_absorbed),
            pml_or_numerical_loss_W=float(pml_or_numerical),
            energy_residual_fraction=float(residual_energy),
            termination_status=status,
            diagnostics=diagnostics,
            stored_fields=stored,
        )


def gaussian_transverse_field(
    u_m,
    v_m,
    waist_m: float,
    *,
    polarization_uv=(1.0 + 0.0j, 0.0 + 0.0j),
    total_flux_W: float = 1.0,
) -> np.ndarray:
    """Return a vacuum forward Gaussian field with specified integrated flux."""
    u = np.asarray(u_m, dtype=float)
    v = np.asarray(v_m, dtype=float)
    if waist_m <= 0.0 or total_flux_W <= 0.0:
        raise ValueError("waist_m and total_flux_W must be positive")
    uu, vv = np.meshgrid(u, v, indexing="ij")
    envelope = np.exp(-(uu * uu + vv * vv) / float(waist_m) ** 2)
    polarization = np.asarray(polarization_uv, dtype=complex)
    norm = float(np.linalg.norm(polarization))
    if polarization.shape != (2,) or norm <= 0.0:
        raise ValueError("polarization_uv must be a non-zero complex pair")
    polarization /= norm
    field = np.zeros((4, u.size, v.size), dtype=complex)
    field[0] = polarization[0] * envelope
    field[1] = polarization[1] * envelope
    field[2] = -field[1] / ETA0_OHM
    field[3] = field[0] / ETA0_OHM
    du = float(np.diff(u)[0])
    dv = float(np.diff(v)[0])
    _, flux = _plane_flux(field, np.ones_like(envelope), du, dv)
    field *= np.sqrt(float(total_flux_W) / flux)
    return field


__all__ = [
    "ETA0_OHM",
    "RayCenteredFrame3D",
    "VectorPE3D",
    "VectorPE3DConfig",
    "VectorPE3DResult",
    "berreman_matrix_uniform",
    "cold_plasma_dielectric_tensor",
    "dielectric_tensor_local",
    "gaussian_transverse_field",
]
