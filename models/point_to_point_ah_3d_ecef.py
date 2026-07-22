"""Seeded point-to-point solver for strict spherical-ECEF AH3D rays."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from utils import rx_enu_residual
from .ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF, trace_ray_ah3d_ecef


class P2PECEFMappingError(RuntimeError):
    """Raised when the local endpoint map is not a stable ground branch."""


def _trajectory_rms_km(left: dict, right: dict) -> float:
    def resample(ray):
        points = np.asarray(ray.get("trajectory_ecef_km", []), dtype=float)
        path = np.asarray(ray.get("trajectory_geom_path_km", []), dtype=float)
        if points.ndim != 2 or points.shape[0] < 2 or path.shape != (len(points),):
            return None
        if path[-1] <= 0.0:
            return None
        q = np.linspace(0.0, path[-1], 41)
        return np.column_stack([
            np.interp(q, path, points[:, axis]) for axis in range(3)])

    a, b = resample(left), resample(right)
    if a is None or b is None:
        return float("inf")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def solve_p2p_ray_ah3d_ecef(
    model: RefractiveIndexAH3DECEF,
    *,
    tx_geocentric,
    rx_geocentric,
    elevation_seed_deg: float,
    bearing_seed_deg: float,
    h_stop_km: float,
    max_path_km: float,
    elevation_bounds_deg=(2.0, 88.0),
    bearing_half_width_deg: float = 5.0,
    d_elevation_deg: float = 0.05,
    d_bearing_deg: float = 0.10,
    jacobian_scales=(0.5, 1.0, 2.0),
    jacobian_scale_rtol: float = 1.0e-2,
    jacobian_condition_max: float = 1.0e10,
    tol_km: float = 5.0,
    max_nfev: int = 20,
    profile: str = "reference",
    medium_time_s: float = 0.0,
    seed_semantics: str = "new_coordinate_seed",
) -> dict | None:
    """Solve the RX ENU endpoint map using a multiscale centered Jacobian."""
    if not isinstance(model, RefractiveIndexAH3DECEF):
        raise TypeError("model must be RefractiveIndexAH3DECEF")
    tx = np.asarray(tx_geocentric, dtype=float)
    rx = np.asarray(rx_geocentric, dtype=float)
    if tx.shape != (3,) or rx.shape != (3,) or not np.all(np.isfinite([tx, rx])):
        raise ValueError("TX/RX must be finite geocentric [lat,lon,height]")
    scales = tuple(float(value) for value in jacobian_scales)
    if len(scales) < 3 or any(value <= 0.0 for value in scales):
        raise ValueError("jacobian_scales requires at least three positive values")
    if any(b <= a for a, b in zip(scales, scales[1:])):
        raise ValueError("jacobian_scales must be strictly increasing")
    cache = {}
    jacobian_diagnostics = []

    def trace(theta):
        key = (float(theta[0]), float(theta[1]))
        if key not in cache:
            cache[key] = trace_ray_ah3d_ecef(
                model,
                tx_lat_deg=float(tx[0]), tx_lon_deg=float(tx[1]),
                tx_height_km=float(tx[2]), elevation_deg=key[0],
                bearing_deg=key[1] % 360.0, h_stop_km=h_stop_km,
                max_path_km=max_path_km, medium_time_s=medium_time_s,
                profile=profile)
        return cache[key]

    def endpoint_residual(ray):
        if ray.get("termination_status") != "ground":
            raise P2PECEFMappingError(
                f"endpoint_map_not_ground:{ray.get('termination_status')}")
        terminal = np.asarray(ray["trajectory_ecef_km"], dtype=float)[-1]
        residual = rx_enu_residual(
            terminal, float(rx[0]), float(rx[1]), float(rx[2]),
            model.earth_radius_km)
        if not np.all(np.isfinite(residual)):
            raise P2PECEFMappingError("endpoint_map_nonfinite")
        return residual

    def residual(theta):
        return endpoint_residual(trace(theta))

    def central(theta, elevation_step, bearing_step):
        ep = endpoint_residual(trace((theta[0] + elevation_step, theta[1])))
        em = endpoint_residual(trace((theta[0] - elevation_step, theta[1])))
        bp = endpoint_residual(trace((theta[0], theta[1] + bearing_step)))
        bm = endpoint_residual(trace((theta[0], theta[1] - bearing_step)))
        return np.column_stack((
            (ep - em) / (2.0 * elevation_step),
            (bp - bm) / (2.0 * bearing_step),
        ))

    def jacobian(theta):
        estimates = [central(
            theta, d_elevation_deg * scale, d_bearing_deg * scale)
            for scale in scales]
        changes = [
            float(np.linalg.norm(estimates[index] - estimates[index - 1])
                  / max(np.linalg.norm(estimates[index]), 1.0e-12))
            for index in range(1, len(estimates))]
        selected = int(np.argmin(changes)) + 1
        value = estimates[selected]
        condition = float(np.linalg.cond(value))
        diagnostics = {
            "scales": list(scales),
            "relative_changes": changes,
            "selected_scale": scales[selected],
            "stable": bool(changes[selected - 1] <= jacobian_scale_rtol),
            "condition": condition,
        }
        jacobian_diagnostics.append(diagnostics)
        if not diagnostics["stable"]:
            raise P2PECEFMappingError("jacobian_reference_unstable")
        if (not np.all(np.isfinite(value)) or not np.isfinite(condition)
                or condition > jacobian_condition_max):
            raise P2PECEFMappingError("jacobian_reference_ill_conditioned")
        return value

    bearing_center = float(bearing_seed_deg)
    lower = np.array([
        float(elevation_bounds_deg[0]),
        bearing_center - float(bearing_half_width_deg),
    ])
    upper = np.array([
        float(elevation_bounds_deg[1]),
        bearing_center + float(bearing_half_width_deg),
    ])
    initial = np.clip(
        np.array([float(elevation_seed_deg), bearing_center]),
        lower + 1.0e-10, upper - 1.0e-10)
    try:
        solved = least_squares(
            residual, initial, jac=jacobian, bounds=(lower, upper),
            method="trf", x_scale="jac", max_nfev=int(max_nfev),
            ftol=1.0e-10, xtol=1.0e-10, gtol=1.0e-10)
        final = dict(trace(solved.x))
        final_residual = endpoint_residual(final)
    except (P2PECEFMappingError, ValueError, FloatingPointError,
            np.linalg.LinAlgError):
        return None
    landing_error = float(np.linalg.norm(final_residual))
    if landing_error >= float(tol_km):
        return None
    final.update({
        "p2p_coordinate_model": "geocentric_TX_RX_ECEF_state_RX_ENU_residual",
        "p2p_residual_semantics": "terminal_minus_RX_in_RX_ENU_East_North_km",
        "p2p_residual_east_north_km": final_residual.tolist(),
        "p2p_landing_error_km": landing_error,
        "p2p_solver": "scipy_trf_bounded",
        "p2p_jacobian_method": "multiscale_centered_RX_ENU",
        "p2p_jacobian_reference": jacobian_diagnostics,
        "p2p_nfev": int(solved.nfev),
        "p2p_njev": int(solved.njev or 0),
        "p2p_trace_profile": str(profile),
        "p2p_seed_semantics": str(seed_semantics),
        "p2p_search_scope": "seeded_local_roots_not_full_domain_mode_audit",
        "rx_geocentric_lat_lon_height": rx.tolist(),
        "convergence_status": "passed",
    })
    return final


def find_seeded_rays_ah3d_ecef(
    model: RefractiveIndexAH3DECEF,
    *,
    tx_geocentric,
    rx_geocentric,
    seeds,
    h_stop_km: float,
    max_path_km: float,
    dedup_launch_deg: float = 0.02,
    dedup_trajectory_rms_km: float = 1.0,
    dedup_apex_km: float = 1.0,
    dedup_tau_ms: float = 0.05,
    **solve_kwargs,
) -> list:
    """Solve and deduplicate explicitly seeded roots without completeness claims."""
    candidates = []
    for index, seed in enumerate(seeds):
        ray = solve_p2p_ray_ah3d_ecef(
            model, tx_geocentric=tx_geocentric,
            rx_geocentric=rx_geocentric,
            elevation_seed_deg=float(seed["elevation_deg"]),
            bearing_seed_deg=float(seed["bearing_deg"]),
            h_stop_km=h_stop_km, max_path_km=max_path_km,
            seed_semantics=str(seed.get(
                "semantics", "legacy_coordinate_reference_seed_only")),
            **solve_kwargs)
        if ray is not None:
            ray["p2p_seed_id"] = str(seed.get("seed_id", index))
            candidates.append(ray)

    unique = []
    for candidate in sorted(candidates, key=lambda item: item["tau_ms"]):
        duplicate = any(
            candidate["branch_id"] == retained["branch_id"]
            and candidate["termination_status"] == retained["termination_status"]
            and candidate["turning_point_count"] == retained["turning_point_count"]
            and abs(candidate["elevation_tx_deg"]
                    - retained["elevation_tx_deg"]) <= dedup_launch_deg
            and abs((candidate["bearing_tx_deg"]
                     - retained["bearing_tx_deg"] + 180.0) % 360.0 - 180.0)
            <= dedup_launch_deg
            and _trajectory_rms_km(candidate, retained)
            <= dedup_trajectory_rms_km
            and abs(candidate["apex_height_km"]
                    - retained["apex_height_km"]) <= dedup_apex_km
            and abs(candidate["tau_ms"] - retained["tau_ms"])
            <= dedup_tau_ms
            for retained in unique)
        if not duplicate:
            unique.append(candidate)
    for ray in unique:
        ray["p2p_dedup_semantics"] = (
            "branch_status_turning_launch_ECEF_trajectory_apex_group_delay")
        ray["p2p_seeded_candidate_count"] = len(candidates)
        ray["p2p_unique_root_count"] = len(unique)
    return unique


__all__ = [
    "P2PECEFMappingError",
    "solve_p2p_ray_ah3d_ecef",
    "find_seeded_rays_ah3d_ecef",
]
