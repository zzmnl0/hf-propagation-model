"""Run and freeze the B1 vector wide-angle PE reference acceptance cases."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import psutil
from scipy.constants import c
from scipy.linalg import eig

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import K_FP
from models.dispersion_ah import appleton_hartree_roots
from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from models.ne_provider_geocentric import VacuumGeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    integrate_fixed_path_ah3d_ecef,
)
from models.vector_pe_3d import (
    ETA0_OHM,
    RayCenteredFrame3D,
    VectorPE3D,
    VectorPE3DConfig,
    berreman_matrix_uniform,
    cold_plasma_dielectric_tensor,
    dielectric_tensor_local,
    gaussian_transverse_field,
)


OUT = ROOT / "output2" / "b1_vector_pe_3d"
FREQUENCY_HZ = 10.0e6
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
WAVELENGTH_M = c / FREQUENCY_HZ
EARTH_RADIUS_KM = 6371.0


class ConstantNeProvider:
    def __init__(self, Ne_m3: float):
        self.Ne_m3 = float(Ne_m3)

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        del ecef_km, earth_radius_km, t
        return self.Ne_m3

    def get_metadata(self):
        return {
            "name": "b1_constant_ne_provider",
            "coordinate_system": "geocentric_spherical",
            "units": "m^-3",
            "Ne_m3": self.Ne_m3,
            "earth_radius_km": EARTH_RADIUS_KM,
        }

    def domain_margin_ecef(self, ecef_km, earth_radius_km):
        del ecef_km
        if float(earth_radius_km) != EARTH_RADIUS_KM:
            raise ValueError("Earth radius mismatch")
        return float("inf")


class AlongTrackLinearNeProvider:
    def __init__(self, origin_ecef_km, tangent_ecef, Ne0_m3, slope_m3_per_m):
        self.origin = np.asarray(origin_ecef_km, dtype=float)
        self.tangent = np.asarray(tangent_ecef, dtype=float)
        self.tangent /= np.linalg.norm(self.tangent)
        self.Ne0 = float(Ne0_m3)
        self.slope = float(slope_m3_per_m)

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        del earth_radius_km, t
        ds_m = 1000.0 * float(
            np.dot(np.asarray(ecef_km, dtype=float) - self.origin, self.tangent)
        )
        return max(self.Ne0 + self.slope * ds_m, 0.0)

    def get_metadata(self):
        return {
            "name": "b1_along_track_linear_ne_provider",
            "coordinate_system": "geocentric_spherical",
            "units": "m^-3",
            "Ne0_m3": self.Ne0,
            "slope_m3_per_m": self.slope,
            "earth_radius_km": EARTH_RADIUS_KM,
        }

    def domain_margin_ecef(self, ecef_km, earth_radius_km):
        del ecef_km
        if float(earth_radius_km) != EARTH_RADIUS_KM:
            raise ValueError("Earth radius mismatch")
        return float("inf")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path, value) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _straight_ray(s_m, rotation=None, translation_km=None) -> dict:
    s = np.asarray(s_m, dtype=float)
    rotation = np.eye(3) if rotation is None else np.asarray(rotation, dtype=float)
    translation = (
        np.zeros(3)
        if translation_km is None
        else np.asarray(translation_km, dtype=float)
    )
    origin = np.array([EARTH_RADIUS_KM, 0.0, 0.0])
    tangent = np.array([0.0, 1.0, 0.0])
    points = origin + s[:, None] * tangent / 1000.0
    points = points @ rotation.T + translation
    tangents = np.tile(tangent @ rotation.T, (s.size, 1))
    return {
        "trajectory_ecef_km": points,
        "trajectory_geom_path_km": s / 1000.0,
        "trajectory_group_tangent_ecef": tangents,
    }


def _circular_ray(s_m, radius_m=100_000.0) -> dict:
    s = np.asarray(s_m, dtype=float)
    theta = s / float(radius_m)
    centre = np.array([EARTH_RADIUS_KM - radius_m / 1000.0, 0.0, 0.0])
    points = np.column_stack(
        [
            centre[0] + radius_m * np.cos(theta) / 1000.0,
            radius_m * np.sin(theta) / 1000.0,
            np.zeros_like(theta),
        ]
    )
    tangents = np.column_stack(
        [-np.sin(theta), np.cos(theta), np.zeros_like(theta)]
    )
    return {
        "trajectory_ecef_km": points,
        "trajectory_geom_path_km": s / 1000.0,
        "trajectory_group_tangent_ecef": tangents,
    }


def _tensor_n2_roots(epsilon, direction) -> np.ndarray:
    khat = np.asarray(direction, dtype=float)
    khat /= np.linalg.norm(khat)
    projector = np.eye(3) - np.outer(khat, khat)
    values = eig(np.asarray(epsilon, dtype=complex), projector, right=False)
    finite = values[np.isfinite(values)]
    finite = finite[np.abs(finite) < 1e6]
    if finite.size != 2:
        raise RuntimeError("cold-plasma wave equation did not yield two finite roots")
    return np.sort_complex(finite)


def _forward_mode(epsilon_suv, mode_index=0, omega=OMEGA):
    matrix = berreman_matrix_uniform(epsilon_suv, omega)
    values, vectors = np.linalg.eig(matrix)
    candidates = []
    for index, value in enumerate(values):
        vector = vectors[:, index]
        eu, ev, hu, hv = vector
        flux = 0.5 * np.real(eu * np.conj(hv) - ev * np.conj(hu))
        if value.real > 0.0 and flux > 0.0:
            vector = vector / np.sqrt(flux)
            candidates.append((float(value.real), vector))
    candidates.sort(key=lambda item: item[0])
    if len(candidates) != 2:
        raise RuntimeError("uniform tensor did not yield two forward modes")
    return candidates[int(mode_index)]


def _gaussian_width(field, u_m, v_m) -> tuple[float, np.ndarray]:
    intensity = np.abs(field[0]) ** 2 + np.abs(field[1]) ** 2
    uu, vv = np.meshgrid(u_m, v_m, indexing="ij")
    total = max(float(np.sum(intensity)), 1e-30)
    centroid = np.array(
        [np.sum(uu * intensity) / total, np.sum(vv * intensity) / total]
    )
    radius2 = (uu - centroid[0]) ** 2 + (vv - centroid[1]) ** 2
    width = float(np.sqrt(2.0 * np.sum(radius2 * intensity) / total))
    return width, centroid


def _run_tensor_and_dispersion():
    Ne = 1.0e11
    B = np.array([0.0, 0.0, 5.0e-5])
    vacuum = cold_plasma_dielectric_tensor(OMEGA, 0.0, B)
    zero_b = cold_plasma_dielectric_tensor(OMEGA, Ne, np.zeros(3))
    collisionless = cold_plasma_dielectric_tensor(OMEGA, Ne, B)
    collisional = cold_plasma_dielectric_tensor(OMEGA, Ne, B, 1.0e4)
    f_hz = FREQUENCY_HZ
    scalar_n2 = 1.0 - K_FP**2 * Ne / (f_hz * f_hz)

    tensor_rows = [
        {
            "metric": "vacuum_identity_max_abs",
            "value": float(np.max(np.abs(vacuum - np.eye(3)))),
            "limit": 1e-12,
        },
        {
            "metric": "B0_isotropic_max_abs",
            "value": float(np.max(np.abs(zero_b - scalar_n2 * np.eye(3)))),
            "limit": 1e-12,
        },
        {
            "metric": "collisionless_hermitian_max_abs",
            "value": float(np.max(np.abs(collisionless - collisionless.conj().T))),
            "limit": 1e-12,
        },
        {
            "metric": "collisional_passivity_negative_min_loss_eigenvalue",
            "value": float(
                max(
                    -np.min(
                        np.linalg.eigvalsh(
                            (collisional - collisional.conj().T) / (2j)
                        )
                    ),
                    0.0,
                )
            ),
            "limit": 1e-12,
        },
    ]

    ah_parallel = appleton_hartree_roots(
        OMEGA, Ne, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 1.0]),
        27_992.489872333 * np.linalg.norm(B), 0.0,
    )
    tensor_parallel = _tensor_n2_roots(collisionless, [0.0, 0.0, 1.0])
    parallel_expected = np.sort_complex(np.asarray(list(ah_parallel.values())))
    ah_perp = appleton_hartree_roots(
        OMEGA, Ne, np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]),
        27_992.489872333 * np.linalg.norm(B), 0.0,
    )
    tensor_perp = _tensor_n2_roots(collisionless, [1.0, 0.0, 0.0])
    perp_expected = np.sort_complex(np.asarray(list(ah_perp.values())))
    tensor_rows.extend(
        [
            {
                "metric": "parallel_AH_root_relative_error",
                "value": float(
                    np.max(
                        np.abs(tensor_parallel - parallel_expected)
                        / np.maximum(np.abs(parallel_expected), 1e-30)
                    )
                ),
                "limit": 1e-10,
            },
            {
                "metric": "perpendicular_AH_root_relative_error",
                "value": float(
                    np.max(
                        np.abs(tensor_perp - perp_expected)
                        / np.maximum(np.abs(perp_expected), 1e-30)
                    )
                ),
                "limit": 1e-10,
            },
        ]
    )

    angle_rows = []
    for angle_deg in (0.0, 10.0, 20.0, 30.0):
        ku = (OMEGA / c) * np.sin(np.radians(angle_deg))
        matrix = berreman_matrix_uniform(np.eye(3), OMEGA, ku, 0.0)
        values = np.linalg.eigvals(matrix)
        positive = np.sort(values.real[values.real > 0.0])
        exact = (OMEGA / c) * np.cos(np.radians(angle_deg))
        relative = float(np.max(np.abs(positive - exact) / exact))
        angle_rows.append(
            {
                "medium": "vacuum",
                "angle_deg": angle_deg,
                "relative_k_s_error": relative,
                "limit": 1e-3,
            }
        )
    return tensor_rows, angle_rows


def _run_geometry():
    s = np.linspace(0.0, 20.0 * WAVELENGTH_M, 21)
    uv = np.linspace(-2.0 * WAVELENGTH_M, 2.0 * WAVELENGTH_M, 17)
    base = RayCenteredFrame3D.from_ah3d(_straight_ray(s), s, uv, uv)
    angle = np.radians(37.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    rotated = RayCenteredFrame3D.from_ah3d(
        _straight_ray(s, rotation=rotation), s, uv, uv
    )
    translated = RayCenteredFrame3D.from_ah3d(
        _straight_ray(s, translation_km=[0.0, 100.0, 0.0]), s, uv, uv
    )
    curved = RayCenteredFrame3D.from_ah3d(
        _circular_ray(s), s, uv, uv
    )

    def orth_error(frame):
        gram = np.einsum(
            "nij,nik->njk", frame.basis_suv_ecef, frame.basis_suv_ecef
        )
        return float(np.max(np.abs(gram - np.eye(3))))

    return {
        "straight_orthonormal_error": orth_error(base),
        "rotated_orthonormal_error": orth_error(rotated),
        "translated_orthonormal_error": orth_error(translated),
        "curved_orthonormal_error": orth_error(curved),
        "straight_rotated_jacobian_max_abs": float(
            np.max(np.abs(base.jacobian - rotated.jacobian))
        ),
        "straight_translated_curvature_max_abs_m_inv": float(
            np.max(
                np.abs(
                    base.curvature_u_m_inv - translated.curvature_u_m_inv
                )
            )
        ),
        "curved_min_jacobian": float(np.min(curved.jacobian)),
        "curved_curvature_m_inv": float(
            np.median(
                np.hypot(curved.curvature_u_m_inv, curved.curvature_v_m_inv)
            )
        ),
        "curved_expected_curvature_m_inv": 1.0e-5,
    }, (base, rotated, curved)


def _run_gaussian(
    grid_points: int,
    ds_lambda: float,
    label: str,
    *,
    half_width_lambda: float = 8.0,
    waist_lambda: float = 2.0,
    length_lambda: float = 8.0,
    pml_cells_override: int | None = None,
):
    half_width = half_width_lambda * WAVELENGTH_M
    axis = np.linspace(-half_width, half_width, grid_points)
    length = length_lambda * WAVELENGTH_M
    n_steps = int(round(length_lambda / ds_lambda))
    s = np.linspace(0.0, length, n_steps + 1)
    frame = RayCenteredFrame3D.from_ah3d(_straight_ray(s), s, axis, axis)
    pml_cells = (
        int(pml_cells_override)
        if pml_cells_override is not None
        else max(5, int(np.ceil(2.0 * WAVELENGTH_M / np.diff(axis)[0])))
    )
    config = VectorPE3DConfig(
        frequency_Hz=FREQUENCY_HZ,
        pml_cells=pml_cells,
        pml_sigma_max=0.5,
        linear_residual_tol=1e-10,
        store_plane_indices=(n_steps // 2,),
    )
    waist = waist_lambda * WAVELENGTH_M
    initial = gaussian_transverse_field(axis, axis, waist)
    started = time.perf_counter()
    result = VectorPE3D(config).solve(
        frame,
        VacuumGeocentricNeProvider3D(EARTH_RADIUS_KM),
        ConstantGeomagneticField3D([0.0, 0.0, 1e-18]),
        initial,
        earth_radius_km=EARTH_RADIUS_KM,
    )
    runtime = time.perf_counter() - started
    width, centroid = _gaussian_width(result.final_field_transverse, axis, axis)
    z_rayleigh = np.pi * waist**2 / WAVELENGTH_M
    expected_width = waist * np.sqrt(1.0 + (length / z_rayleigh) ** 2)
    return {
        "label": label,
        "grid_points": grid_points,
        "du_lambda": float(np.diff(axis)[0] / WAVELENGTH_M),
        "ds_lambda": ds_lambda,
        "pml_cells": pml_cells,
        "half_width_lambda": half_width_lambda,
        "waist_lambda": waist_lambda,
        "length_lambda": length_lambda,
        "runtime_s": runtime,
        "width_m": width,
        "expected_width_m": expected_width,
        "width_relative_error": abs(width - expected_width) / expected_width,
        "centroid_error_lambda": float(np.linalg.norm(centroid) / WAVELENGTH_M),
        "flux_change_fraction": float(
            abs(result.flux_history_W[-1] - result.flux_history_W[0])
            / result.flux_history_W[0]
        ),
        "linear_residual_max": result.linear_residual_max,
        "out_of_cone_fraction": result.angular_out_of_cone_power_fraction,
        "backward_fraction": result.backward_power_fraction,
        "termination_status": result.termination_status,
        "axis_m": axis,
        "initial": initial,
        "result": result,
    }


def _run_uniform_magnetized_modes():
    Ne = 8.0e10
    B = np.array([0.0, 5.0e-5, 0.0])
    epsilon_ecef = cold_plasma_dielectric_tensor(OMEGA, Ne, B)
    basis = np.column_stack(
        [np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]),
         np.array([0.0, 0.0, -1.0])]
    )
    epsilon = dielectric_tensor_local(epsilon_ecef, basis)
    rows = []
    for mode in (0, 1):
        ks, vector = _forward_mode(epsilon, mode)
        ds = WAVELENGTH_M
        propagator = np.linalg.solve(
            np.eye(4)
            - 0.5j
            * ds
            * (berreman_matrix_uniform(epsilon, OMEGA) - OMEGA / c * np.eye(4)),
            np.eye(4)
            + 0.5j
            * ds
            * (berreman_matrix_uniform(epsilon, OMEGA) - OMEGA / c * np.eye(4)),
        )
        output = np.linalg.matrix_power(propagator, 8) @ vector
        coefficient = np.vdot(vector, output) / np.vdot(vector, vector)
        residual = output - coefficient * vector
        leakage = float(
            np.vdot(residual, residual).real
            / max(np.vdot(output, output).real, 1e-30)
        )
        exact_phase = np.exp(1j * (ks - OMEGA / c) * ds * 8)
        phase_error = float(abs(np.angle(coefficient / exact_phase)))
        rows.append(
            {
                "mode": mode,
                "k_s_rad_m": ks,
                "cross_polarization_power": leakage,
                "pade_phase_error_rad": phase_error,
            }
        )
    return rows


def _run_wkb_light():
    rows = []
    grid_points = 41
    axis = np.linspace(-10.0 * WAVELENGTH_M, 10.0 * WAVELENGTH_M, grid_points)
    length = 4.0 * WAVELENGTH_M
    s = np.linspace(0.0, length, 5)
    ray = _straight_ray(s)
    frame = RayCenteredFrame3D.from_ah3d(ray, s, axis, axis)
    origin = ray["trajectory_ecef_km"][0]
    tangent = ray["trajectory_group_tangent_ecef"][0]
    provider = AlongTrackLinearNeProvider(
        origin, tangent, 5.0e10, 5.0e10 / (400.0 * WAVELENGTH_M)
    )
    B_provider = ConstantGeomagneticField3D([0.0, 5.0e-5, 0.0])
    pml_cells = 5
    envelope = np.exp(
        -(
            np.meshgrid(axis, axis, indexing="ij")[0] ** 2
            + np.meshgrid(axis, axis, indexing="ij")[1] ** 2
        )
        / (4.0 * WAVELENGTH_M) ** 2
    )
    basis = frame.basis_suv_ecef[0]
    for mode in (0, 1):
        frequency_runs = []
        for frequency in FREQUENCY_HZ * np.array([1.0 - 1e-4, 1.0, 1.0 + 1e-4]):
            omega = 2.0 * np.pi * frequency
            epsilon0 = dielectric_tensor_local(
                cold_plasma_dielectric_tensor(
                    omega, provider.Ne0, [0.0, 5.0e-5, 0.0]
                ),
                basis,
            )
            _, modal = _forward_mode(epsilon0, mode, omega=omega)
            initial = modal[:, None, None] * envelope[None, :, :]
            config = VectorPE3DConfig(
                frequency_Hz=float(frequency),
                pml_cells=pml_cells,
                pml_sigma_max=0.4,
            )
            result = VectorPE3D(config).solve(
                frame,
                provider,
                B_provider,
                initial,
                earth_radius_km=EARTH_RADIUS_KM,
            )
            centre_ratio = (
                result.final_field_transverse[
                    0, grid_points // 2, grid_points // 2
                ]
                / initial[0, grid_points // 2, grid_points // 2]
            )
            sample_s = np.linspace(0.0, length, 257)
            ks_values = []
            for distance in sample_s:
                point = origin + tangent * distance / 1000.0
                epsilon_local = dielectric_tensor_local(
                    cold_plasma_dielectric_tensor(
                        omega,
                        provider.get_Ne_ecef(point, EARTH_RADIUS_KM),
                        [0.0, 5.0e-5, 0.0],
                    ),
                    basis,
                )
                ks_values.append(
                    _forward_mode(epsilon_local, mode, omega=omega)[0]
                )
            z_rayleigh = float(
                ks_values[0] * (4.0 * WAVELENGTH_M) ** 2 / 2.0
            )
            gouy_phase = float(np.arctan(length / z_rayleigh))
            expected_envelope_phase = float(
                np.trapezoid(
                    np.asarray(ks_values) - omega / c, sample_s
                )
            ) - gouy_phase
            frequency_runs.append(
                {
                    "frequency_Hz": float(frequency),
                    "omega": float(omega),
                    "initial": initial,
                    "result": result,
                    "centre_ratio": centre_ratio,
                    "ks_values": np.asarray(ks_values),
                    "gouy_phase": gouy_phase,
                    "expected_envelope_phase": expected_envelope_phase,
                }
            )
        central = frequency_runs[1]
        initial = central["initial"]
        result = central["result"]
        _, initial_centroid = _gaussian_width(initial, axis, axis)
        _, centroid = _gaussian_width(result.final_field_transverse, axis, axis)
        phase_error = float(
            abs(
                np.angle(
                    central["centre_ratio"]
                    * np.exp(-1j * central["expected_envelope_phase"])
                )
            )
        )
        physical_actual = np.unwrap(
            [
                np.angle(run["centre_ratio"])
                + run["omega"] / c * length
                for run in frequency_runs
            ]
        )
        physical_expected = np.unwrap(
            [
                run["expected_envelope_phase"] + run["omega"] / c * length
                for run in frequency_runs
            ]
        )
        delta_omega = frequency_runs[2]["omega"] - frequency_runs[0]["omega"]
        group_delay_s = float(
            (physical_actual[2] - physical_actual[0]) / delta_omega
        )
        expected_group_delay_s = float(
            (physical_expected[2] - physical_expected[0]) / delta_omega
        )
        group_delay_relative_error = float(
            abs(group_delay_s - expected_group_delay_s)
            / max(abs(expected_group_delay_s), 1e-30)
        )
        physical_phase = float(
            np.trapezoid(central["ks_values"], sample_s)
            - central["gouy_phase"]
        )
        launch_n2 = float((central["ks_values"][0] / (OMEGA / c)) ** 2)
        ah_roots = appleton_hartree_roots(
            OMEGA,
            provider.Ne0,
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            27_992.489872333 * 5.0e-5,
            0.0,
        )
        ah_branch = min(
            ("O", "X"), key=lambda branch: abs(ah_roots[branch].real - launch_n2)
        )
        ah_model = RefractiveIndexAH3DECEF(
            provider,
            earth_radius_km=EARTH_RADIUS_KM,
            freq_MHz=FREQUENCY_HZ / 1e6,
            wave_mode=ah_branch,
            geomagnetic_provider=B_provider,
            integration_profile="reference",
        )
        ah_ray = integrate_fixed_path_ah3d_ecef(
            ah_model,
            tx_lat_deg=0.0,
            tx_lon_deg=0.0,
            tx_height_km=0.0,
            elevation_deg=0.0,
            bearing_deg=90.0,
            path_length_km=length / 1000.0,
            stop_on_physical_events=False,
        )
        ah_group_delay_s = float(ah_ray["tau_ms"]) / 1000.0
        ah_group_delay_relative_error = float(
            abs(group_delay_s - ah_group_delay_s)
            / max(abs(ah_group_delay_s), 1e-30)
        )
        rows.append(
            {
                "mode": mode,
                "ah_branch": ah_branch,
                "n2_min": float(
                    min(
                        (value / (OMEGA / c)) ** 2
                        for value in central["ks_values"]
                    )
                ),
                "centroid_shift_lambda": float(
                    np.linalg.norm(centroid - initial_centroid) / WAVELENGTH_M
                ),
                "phase_error_rad": phase_error,
                "phase_relative_error": phase_error / abs(physical_phase),
                "gouy_phase_rad": central["gouy_phase"],
                "group_delay_s": group_delay_s,
                "expected_group_delay_s": expected_group_delay_s,
                "group_delay_relative_error": group_delay_relative_error,
                "ah3d_group_delay_s": ah_group_delay_s,
                "ah3d_group_delay_relative_error": ah_group_delay_relative_error,
                "flux_change_fraction": float(
                    abs(result.flux_history_W[-1] - result.flux_history_W[0])
                    / result.flux_history_W[0]
                ),
                "termination_status": result.termination_status,
            }
        )
    return rows


def _run_pml_reflection():
    axis = np.linspace(-12.0 * WAVELENGTH_M, 12.0 * WAVELENGTH_M, 97)
    s = np.linspace(0.0, 32.0 * WAVELENGTH_M, 33)
    frame = RayCenteredFrame3D.from_ah3d(_straight_ray(s), s, axis, axis)
    field = gaussian_transverse_field(
        axis, axis, 3.0 * WAVELENGTH_M, polarization_uv=(0.0, 1.0)
    )
    angle_deg = 15.0
    phase = np.exp(
        1j
        * (2.0 * np.pi / WAVELENGTH_M)
        * np.sin(np.radians(angle_deg))
        * axis
    )[:, None]
    field *= phase
    field[2] *= np.cos(np.radians(angle_deg))
    result = VectorPE3D(
        VectorPE3DConfig(
            frequency_Hz=FREQUENCY_HZ,
            pml_cells=16,
            pml_sigma_max=1.0,
        )
    ).solve(
        frame,
        VacuumGeocentricNeProvider3D(EARTH_RADIUS_KM),
        ConstantGeomagneticField3D([0.0, 0.0, 1e-18]),
        field,
        earth_radius_km=EARTH_RADIUS_KM,
    )
    spectrum = np.fft.fftshift(
        np.fft.fft2(result.final_field_transverse[1])
    )
    ku = 2.0 * np.pi * np.fft.fftshift(
        np.fft.fftfreq(axis.size, np.diff(axis)[0])
    )
    reflected = float(np.sum(np.abs(spectrum[ku < 0.0, :]) ** 2))
    incident = float(np.sum(np.abs(np.fft.fft2(field[1])) ** 2))
    reflection_dB = float(10.0 * np.log10(max(reflected / incident, 1e-30)))
    return {
        "angle_deg": angle_deg,
        "reflection_dB": reflection_dB,
        "limit_dB": -40.0,
        "out_of_cone_fraction": result.angular_out_of_cone_power_fraction,
        "termination_status": result.termination_status,
        "input_flux_W": float(result.flux_history_W[0]),
        "output_flux_W": float(result.flux_history_W[-1]),
        "axis_m": axis,
        "result": result,
    }


def _run_collision_energy():
    axis = np.linspace(-6.0 * WAVELENGTH_M, 6.0 * WAVELENGTH_M, 25)
    s = np.linspace(0.0, 4.0 * WAVELENGTH_M, 9)
    frame = RayCenteredFrame3D.from_ah3d(_straight_ray(s), s, axis, axis)
    Ne = 5.0e10
    collision_s = 1.0e5
    epsilon = dielectric_tensor_local(
        cold_plasma_dielectric_tensor(OMEGA, Ne, np.zeros(3), collision_s),
        frame.basis_suv_ecef[0],
    )
    _, modal = _forward_mode(epsilon, 0)
    uu, vv = np.meshgrid(axis, axis, indexing="ij")
    envelope = np.exp(
        -(uu * uu + vv * vv) / (2.0 * WAVELENGTH_M) ** 2
    )
    field = modal[:, None, None] * envelope[None, :, :]
    result = VectorPE3D(
        VectorPE3DConfig(
            frequency_Hz=FREQUENCY_HZ,
            pml_cells=5,
            pml_sigma_max=0.4,
        )
    ).solve(
        frame,
        ConstantNeProvider(Ne),
        ConstantGeomagneticField3D([0.0, 0.0, 0.0]),
        field,
        earth_radius_km=EARTH_RADIUS_KM,
        collision_sampler=lambda _point: collision_s,
    )
    flux_loss = float(result.flux_history_W[0] - result.flux_history_W[-1])
    balance = float(
        abs(flux_loss - result.collision_absorbed_W)
        / max(abs(result.flux_history_W[0]), 1e-30)
    )
    return {
        "collision_s": collision_s,
        "input_flux_W": float(result.flux_history_W[0]),
        "output_flux_W": float(result.flux_history_W[-1]),
        "flux_loss_W": flux_loss,
        "integrated_collision_absorbed_W": result.collision_absorbed_W,
        "balance_error_fraction": balance,
        "limit": 0.01,
        "termination_status": result.termination_status,
    }


def _plot_outputs(tensor_rows, angle_rows, geometry, frames, gaussian_rows,
                  magnetized_rows, wkb_rows, pml_row, collision_row):
    plt.figure(figsize=(8, 4.8))
    labels = [row["metric"] for row in tensor_rows]
    values = [max(abs(row["value"]), 1e-18) for row in tensor_rows]
    limits = [max(abs(row["limit"]), 1e-18) for row in tensor_rows]
    x = np.arange(len(labels))
    plt.semilogy(x, values, "o-", label="measured")
    plt.semilogy(x, limits, "k--", label="limit")
    plt.xticks(x, labels, rotation=25, ha="right")
    plt.ylabel("absolute / relative error")
    plt.title("B1 cold-plasma tensor acceptance")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "b1_tensor_eigen_acceptance.png", dpi=180)
    plt.close()

    plt.figure(figsize=(6.5, 4.5))
    plt.plot(
        [row["angle_deg"] for row in angle_rows],
        [row["relative_k_s_error"] for row in angle_rows],
        "o-",
    )
    plt.axhline(1e-3, color="k", linestyle="--", label="1e-3 gate")
    plt.yscale("log")
    plt.xlabel("local propagation angle (deg)")
    plt.ylabel("relative longitudinal wavenumber error")
    plt.title("B1 continuum vector wide-angle dispersion")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "b1_wide_angle_dispersion.png", dpi=180)
    plt.close()

    base, rotated, curved = frames
    fig = plt.figure(figsize=(10, 4.2))
    ax1 = fig.add_subplot(121, projection="3d")
    for frame, label in ((base, "straight"), (rotated, "rotated"), (curved, "curved")):
        points = frame.centre_ecef_km - frame.centre_ecef_km[0]
        ax1.plot(points[:, 0], points[:, 1], points[:, 2], label=label)
    ax1.set_title("Bishop centre curves")
    ax1.set_xlabel("ECEF dx (km)")
    ax1.set_ylabel("ECEF dy (km)")
    ax1.set_zlabel("ECEF dz (km)")
    ax1.legend()
    ax2 = fig.add_subplot(122)
    ax2.plot(curved.s_m / WAVELENGTH_M, np.min(curved.jacobian, axis=(1, 2)))
    ax2.axhline(0.8, color="k", linestyle="--")
    ax2.set_xlabel("s / wavelength")
    ax2.set_ylabel("minimum Jacobian")
    ax2.set_title("Curved-frame Jacobian")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "b1_bishop_frame_jacobian.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(
        len(gaussian_rows), 3, figsize=(11, 3.5 * len(gaussian_rows)), squeeze=False
    )
    for row_index, row in enumerate(gaussian_rows):
        axis_lambda = row["axis_m"] / WAVELENGTH_M
        initial_i = np.sum(np.abs(row["initial"][:2]) ** 2, axis=0)
        final_i = np.sum(
            np.abs(row["result"].final_field_transverse[:2]) ** 2, axis=0
        )
        axes[row_index, 0].imshow(
            initial_i.T,
            origin="lower",
            extent=[axis_lambda[0], axis_lambda[-1], axis_lambda[0], axis_lambda[-1]],
        )
        axes[row_index, 0].set_title(f"{row['label']} input intensity")
        axes[row_index, 1].imshow(
            final_i.T,
            origin="lower",
            extent=[axis_lambda[0], axis_lambda[-1], axis_lambda[0], axis_lambda[-1]],
        )
        axes[row_index, 1].set_title(
            f"output; width err={row['width_relative_error']:.3g}"
        )
        axes[row_index, 2].plot(
            np.arange(row["result"].flux_history_W.size) * row["ds_lambda"],
            row["result"].flux_history_W,
        )
        axes[row_index, 2].set_title("integrated Poynting flux")
        axes[row_index, 2].set_xlabel("s / wavelength")
        axes[row_index, 2].set_ylabel("W")
        axes[row_index, 2].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "b1_gaussian_field_and_flux.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].bar(
        [f"M{row['mode']}" for row in magnetized_rows],
        [row["cross_polarization_power"] for row in magnetized_rows],
    )
    axes[0].axhline(1e-4, color="k", linestyle="--")
    axes[0].set_yscale("log")
    axes[0].set_title("Uniform magnetized mode leakage")
    axes[0].set_ylabel("power fraction")
    axes[1].bar(
        [f"M{row['mode']}" for row in wkb_rows],
        [row["phase_error_rad"] for row in wkb_rows],
    )
    axes[1].axhline(0.005 * 2.0 * np.pi, color="k", linestyle="--")
    axes[1].set_title("Smooth-medium WKB phase error")
    axes[1].set_ylabel("rad")
    fig.tight_layout()
    fig.savefig(OUT / "b1_wkb_mode_acceptance.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].bar(["measured"], [pml_row["reflection_dB"]])
    axes[0].axhline(pml_row["limit_dB"], color="k", linestyle="--")
    axes[0].set_ylabel("reflected power (dB)")
    axes[0].set_title("15 deg transverse PML reflection")
    axes[1].bar(
        ["flux loss", "collision integral"],
        [collision_row["flux_loss_W"], collision_row["integrated_collision_absorbed_W"]],
    )
    axes[1].set_ylabel("W")
    axes[1].set_title(
        f"collision balance error={collision_row['balance_error_fraction']:.3g}"
    )
    fig.tight_layout()
    fig.savefig(OUT / "b1_pml_and_energy_budget.png", dpi=180)
    plt.close(fig)

def _plot_acceptance_matrix_from_metrics(metrics):
    gates = metrics["gates"]
    names = [
        "tensor", "angle", "geometry", "gaussian", "mode", "wkb",
        "pml", "energy", "full lambda/8",
    ]
    keys = [
        "tensor",
        "wide_angle_30deg",
        "bishop_geometry",
        "planned_w0_6lambda_25lambda_full_grid",
        "uniform_mode_leakage",
        "light_wkb",
        "measured_pml_reflection_minus40dB",
        "collision_energy_budget",
        "planned_w0_6lambda_25lambda_lambda8_refined",
    ]
    values = [1.0 if gates.get(key) else 0.0 for key in keys]
    if (
        not gates.get("planned_w0_6lambda_25lambda_lambda8_refined")
        and metrics.get("resource_assessment", {}).get("resource_blocked")
    ):
        values[-1] = 0.5
    plt.figure(figsize=(8, 3.8))
    plt.imshow(np.asarray(values)[None, :], vmin=0, vmax=1, cmap="RdYlGn")
    plt.yticks([0], ["gate"])
    plt.xticks(np.arange(len(names)), names, rotation=25, ha="right")
    for index, value in enumerate(values):
        text = "PASS" if value == 1.0 else ("BLOCKED" if value == 0.5 else "FAIL")
        plt.text(index, 0, text, ha="center", va="center")
    plt.title("B1 reference acceptance matrix")
    plt.tight_layout()
    plt.savefig(OUT / "b1_acceptance_matrix.png", dpi=180)
    plt.close()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    tensor_rows, angle_rows = _run_tensor_and_dispersion()
    geometry, frames = _run_geometry()
    magnetized_rows = _run_uniform_magnetized_modes()
    gaussian_rows = [
        _run_gaussian(65, 1.0, "lambda4_base"),
        _run_gaussian(65, 0.5, "lambda4_half_step"),
        _run_gaussian(129, 1.0, "lambda8_refined"),
        _run_gaussian(
            81, 1.0, "domain_expanded_25pct_each_side",
            half_width_lambda=10.0,
        ),
        _run_gaussian(
            65, 1.0, "pml_thickened", pml_cells_override=13,
        ),
        _run_gaussian(
            185,
            1.0,
            "planned_w0_6lambda_25lambda_base",
            half_width_lambda=23.0,
            waist_lambda=6.0,
            length_lambda=25.0,
            pml_cells_override=16,
        ),
    ]
    wkb_rows = _run_wkb_light()
    pml_row = _run_pml_reflection()
    collision_row = _run_collision_energy()

    tensor_pass = all(row["value"] <= row["limit"] for row in tensor_rows)
    angle_pass = all(
        row["relative_k_s_error"] <= row["limit"] for row in angle_rows
    )
    geometry_pass = (
        max(geometry[key] for key in geometry if "orthonormal" in key) <= 1e-10
        and geometry["curved_min_jacobian"] >= 0.8
        and abs(
            geometry["curved_curvature_m_inv"]
            - geometry["curved_expected_curvature_m_inv"]
        )
        <= 2e-7
    )
    mode_pass = all(
        row["cross_polarization_power"] <= 1e-4
        and row["pade_phase_error_rad"] <= 0.01
        for row in magnetized_rows
    )
    light_gaussian_pass = all(
        row["width_relative_error"] <= (
            0.01 if row["label"] == "lambda8_refined" else 0.02
        )
        and row["centroid_error_lambda"] <= 0.1
        and row["flux_change_fraction"] <= 0.005
        and row["linear_residual_max"] <= 1e-10
        and row["termination_status"] == "ok"
        for row in gaussian_rows
    )
    step_change = float(
        abs(
            gaussian_rows[0]["width_m"] - gaussian_rows[1]["width_m"]
        )
        / gaussian_rows[1]["width_m"]
    )
    step_pass = step_change <= 0.01
    base_field = gaussian_rows[0]["result"].final_field_transverse
    refined_on_base = gaussian_rows[2]["result"].final_field_transverse[
        :, ::2, ::2
    ]
    grid_field_change = float(
        np.linalg.norm(base_field - refined_on_base)
        / max(np.linalg.norm(refined_on_base), 1e-30)
    )
    grid_flux_change = float(
        abs(
            gaussian_rows[0]["result"].flux_history_W[-1]
            - gaussian_rows[2]["result"].flux_history_W[-1]
        )
        / gaussian_rows[2]["result"].flux_history_W[-1]
    )
    grid_half_gate = grid_field_change <= 0.02 and grid_flux_change <= 0.005
    expanded_on_base = gaussian_rows[3]["result"].final_field_transverse[
        :, 8:73, 8:73
    ]
    domain_field_change = float(
        np.linalg.norm(base_field - expanded_on_base)
        / max(np.linalg.norm(expanded_on_base), 1e-30)
    )
    domain_gate = domain_field_change <= 0.01
    pml_field = gaussian_rows[4]["result"].final_field_transverse
    core = np.s_[:, 13:-13, 13:-13]
    pml_field_change = float(
        np.linalg.norm(base_field[core] - pml_field[core])
        / max(np.linalg.norm(pml_field[core]), 1e-30)
    )
    pml_flux_change = float(
        abs(
            gaussian_rows[0]["result"].flux_history_W[-1]
            - gaussian_rows[4]["result"].flux_history_W[-1]
        )
        / gaussian_rows[4]["result"].flux_history_W[-1]
    )
    pml_thickness_gate = (
        pml_field_change <= 0.005 and pml_flux_change <= 0.005
    )
    wkb_pass = all(
        row["n2_min"] >= 0.2
        and row["centroid_shift_lambda"] <= 0.5
        and row["phase_relative_error"] <= 0.005
        and row["group_delay_relative_error"] <= 0.01
        and row["ah3d_group_delay_relative_error"] <= 0.01
        and row["termination_status"] == "ok"
        for row in wkb_rows
    )
    full_row = gaussian_rows[5]
    full_gaussian_gate = (
        full_row["width_relative_error"] <= 0.02
        and full_row["centroid_error_lambda"] <= 0.1
        and full_row["flux_change_fraction"] <= 0.005
        and full_row["linear_residual_max"] <= 1e-10
        and full_row["termination_status"] == "ok"
    )
    three_frequency_group_delay_gate = all(
        row["group_delay_relative_error"] <= 0.01 for row in wkb_rows
    )
    pml_reflection_gate = (
        pml_row["reflection_dB"] <= pml_row["limit_dB"]
        and pml_row["termination_status"] == "ok"
    )
    collision_energy_gate = (
        collision_row["balance_error_fraction"] <= collision_row["limit"]
        and collision_row["termination_status"] == "ok"
    )
    full_refined_points = 369
    full_refined_unknowns = 4 * full_refined_points**2
    available_memory_bytes = int(psutil.virtual_memory().available)
    estimated_refined_lu_memory_bytes = 12_000_000_000
    resource_blocked = available_memory_bytes < estimated_refined_lu_memory_bytes
    full_refined_gate = False
    farr_status = "farr_optional_evidence_resource_deferred"
    all_nonresource = all(
        [
            tensor_pass,
            angle_pass,
            geometry_pass,
            mode_pass,
            light_gaussian_pass,
            step_pass,
            wkb_pass,
            full_gaussian_gate,
            grid_half_gate,
            domain_gate,
            pml_thickness_gate,
            three_frequency_group_delay_gate,
            pml_reflection_gate,
            collision_energy_gate,
        ]
    )
    if all_nonresource and full_refined_gate:
        final_status = "passed_b1_vector_wide_angle_pe_reference"
    elif all_nonresource and resource_blocked:
        final_status = "blocked_by_cpu_memory_or_solver_stability"
    else:
        final_status = "failed_b1_vector_pe_acceptance"

    serial_gaussian = []
    for row in gaussian_rows:
        serial_gaussian.append(
            {
                key: value
                for key, value in row.items()
                if key not in {"axis_m", "initial", "result"}
            }
        )
    metrics = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "final_status": final_status,
        "farr_status": farr_status,
        "semantics": {
            "solver_dimensionality": "full_3d_two_transverse_vector_pe",
            "propagation_direction": "one_way_forward",
            "validity": "research_reference_only",
            "calibration_status": "reference_uncalibrated",
            "external_c1_status": "failed_c1_chapman_high_low_external_comparison",
            "external_pharlap_status": (
                "external_pharlap_comparison_paused_due_effective_medium_inconsistency"
            ),
        },
        "tensor": tensor_rows,
        "wide_angle": angle_rows,
        "geometry": geometry,
        "magnetized_modes": magnetized_rows,
        "gaussian_light": serial_gaussian,
        "wkb_light": wkb_rows,
        "pml_reflection": {
            key: value for key, value in pml_row.items()
            if key not in {"axis_m", "result"}
        },
        "collision_energy": collision_row,
        "step_half_width_change": step_change,
        "grid_lambda4_lambda8_field_l2": grid_field_change,
        "grid_lambda4_lambda8_flux_change": grid_flux_change,
        "domain_expansion_field_l2": domain_field_change,
        "pml_thickness_field_l2": pml_field_change,
        "pml_thickness_flux_change": pml_flux_change,
        "resource_assessment": {
            "full_refined_grid_points": full_refined_points,
            "full_refined_complex_unknowns": full_refined_unknowns,
            "available_memory_bytes": available_memory_bytes,
            "estimated_refined_lu_memory_bytes": estimated_refined_lu_memory_bytes,
            "resource_blocked": resource_blocked,
        },
        "gates": {
            "tensor": tensor_pass,
            "wide_angle_30deg": angle_pass,
            "bishop_geometry": geometry_pass,
            "uniform_mode_leakage": mode_pass,
            "light_gaussian": light_gaussian_pass,
            "step_halving": step_pass,
            "light_wkb": wkb_pass,
            "planned_w0_6lambda_25lambda_full_grid": full_gaussian_gate,
            "planned_w0_6lambda_25lambda_lambda8_refined": full_refined_gate,
            "planned_grid_lambda4_lambda8": grid_half_gate,
            "planned_domain_expansion_50pct": domain_gate,
            "planned_pml_thickness_change": pml_thickness_gate,
            "planned_three_frequency_group_delay": three_frequency_group_delay_gate,
            "measured_pml_reflection_minus40dB": pml_reflection_gate,
            "collision_energy_budget": collision_energy_gate,
        },
        "runtime_s": time.perf_counter() - started,
    }
    _json(OUT / "b1_metrics.json", metrics)
    with (OUT / "b1_tensor_metrics.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "limit"])
        writer.writeheader()
        writer.writerows(tensor_rows)
    np.savez_compressed(
        OUT / "b1_light_fields.npz",
        wavelength_m=WAVELENGTH_M,
        base_initial=gaussian_rows[0]["initial"],
        base_final=gaussian_rows[0]["result"].final_field_transverse,
        base_flux=gaussian_rows[0]["result"].flux_history_W,
        half_step_final=gaussian_rows[1]["result"].final_field_transverse,
        half_step_flux=gaussian_rows[1]["result"].flux_history_W,
        refined_final=gaussian_rows[2]["result"].final_field_transverse,
        refined_flux=gaussian_rows[2]["result"].flux_history_W,
    )
    _plot_outputs(
        tensor_rows,
        angle_rows,
        geometry,
        frames,
        gaussian_rows,
        magnetized_rows,
        wkb_rows,
        pml_row,
        collision_row,
    )
    _plot_acceptance_matrix_from_metrics(metrics)
    source_paths = [
        ROOT / "models" / "vector_pe_3d.py",
        ROOT / "models" / "dispersion_ah.py",
        ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        ROOT / "models" / "ne_provider_geocentric.py",
        ROOT / "models" / "geomagnetic_field_3d.py",
        ROOT / "tests" / "_run_b1_vector_pe_3d.py",
        ROOT / "tests" / "_verify_b1_vector_pe_3d.py",
        ROOT / "word" / "b1_vector_wide_angle_pe_3d_implementation_plan.txt",
        ROOT / "word" / "a4_c1g_spherical_high_ray_equation_and_parameter_attribution_report.txt",
    ]
    artifact_paths = sorted(
        path for path in OUT.iterdir() if path.name != "b1_manifest.json"
    )
    manifest = {
        "status": final_status,
        "sources": {
            str(path.relative_to(ROOT)): _sha256(path) for path in source_paths
        },
        "artifacts": {
            str(path.relative_to(ROOT)): _sha256(path) for path in artifact_paths
        },
    }
    _json(OUT / "b1_manifest.json", manifest)
    print(final_status)
    print(f"runtime_s={metrics['runtime_s']:.3f}")


if __name__ == "__main__":
    main()
