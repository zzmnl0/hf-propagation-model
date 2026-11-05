"""AH ray geometry and action between a physical port and an ECEF plane.

The action is a geometric-optics line integral.  Polarization, amplitude and
turning-point phase require separate field transport.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import DOP853
from scipy.optimize import brentq

from models.ray_tracer_ah_3d_ecef import (
    ah3d_ecef_rhs,
    hamiltonian_residual_ah3d_ecef,
    project_dispersion_ah3d_ecef,
)
from .config import ONEWAY_L2_SETTINGS


class _ProviderDomainBoundary(RuntimeError):
    """A stage or its AH finite-difference stencil is outside the Ne domain."""


def trace_to_plane(model, position_ecef_km, wavevector_ecef_rad_km,
                   plane_point_ecef_km, plane_normal_ecef, *, max_path_km,
                   direction=0, medium_time_s=0.0, profile="L2",
                   max_step_km=None):
    """Trace from an on-shell AH port to the first directed plane crossing.

    Distances and wavevectors use km and rad/km, group time uses seconds, and
    action is in radians.  A starting point on the target plane is ignored until
    the ray leaves it and returns.  ``direction`` is the crossing sign of the
    signed plane distance and can be -1, 0, or +1.  The isolated route uses
    its independent L2 integration profile.
    """
    position = np.asarray(position_ecef_km, dtype=float)
    wavevector = np.asarray(wavevector_ecef_rad_km, dtype=float)
    plane_point = np.asarray(plane_point_ecef_km, dtype=float)
    normal = np.asarray(plane_normal_ecef, dtype=float)
    vectors = (position, wavevector, plane_point, normal)
    if any(value.shape != (3,) or not np.all(np.isfinite(value)) for value in vectors):
        raise ValueError("port and plane must contain finite ECEF three-vectors")
    normal_length = float(np.linalg.norm(normal))
    if normal_length == 0 or np.linalg.norm(wavevector) == 0:
        raise ValueError("plane normal and wavevector must be nonzero")
    normal = normal / normal_length
    length = float(max_path_km)
    time_s = float(medium_time_s)
    if not np.isfinite(length) or length <= 0 or not np.isfinite(time_s):
        raise ValueError("max_path_km must be positive and medium_time_s finite")
    if direction not in (-1, 0, 1):
        raise ValueError("direction must be -1, 0, or 1")
    if model.provider_domain_margin(position) < 0:
        raise ValueError("initial port lies outside the Ne provider domain")

    initial = np.r_[position, wavevector, 0.0, 0.0]
    if hamiltonian_residual_ah3d_ecef(initial[:7], model, t=time_s) > 1e-8:
        raise ValueError("initial wavevector is not on the selected AH dispersion shell")
    if profile != "L2":
        raise ValueError(f"unknown one-way integration profile: {profile}")
    settings = dict(ONEWAY_L2_SETTINGS)
    if any(not np.isfinite(value) or value <= 0 for value in settings.values()):
        raise ValueError("one-way integration settings must be finite and positive")
    max_step = float(settings["ds_max_km"] if max_step_km is None else max_step_km)
    if not np.isfinite(max_step) or max_step <= 0:
        raise ValueError("max_step_km must be finite and positive")
    atol_k = float(settings.get("atol_k", model.k0 * settings.get("atol_p", 0.0)))
    atol = np.array([settings["atol_r_km"]] * 3 + [atol_k] * 3 + [1e-12, 1e-8])

    def rhs(path_km, state):
        center_margin = model.provider_domain_margin(state[:3])
        if np.isnan(center_margin):
            raise ValueError("Ne provider returned a non-finite domain margin")
        if center_margin < 0:
            raise _ProviderDomainBoundary("AH stage is outside the Ne provider domain")
        if np.isfinite(center_margin):
            for scale in model.spatial_fd_scales:
                step = model.spatial_fd_base_km * scale
                for axis in range(3):
                    for sign in (-1, 1):
                        sample = state[:3].copy()
                        sample[axis] += sign * step
                        margin = model.provider_domain_margin(sample)
                        if np.isnan(margin):
                            raise ValueError("Ne provider returned a non-finite domain margin")
                        if margin < 0:
                            raise _ProviderDomainBoundary(
                                "AH spatial-difference stencil crossed the Ne provider domain")
        ray_rhs = ah3d_ecef_rhs(path_km, state[:7], model, medium_time_s=time_s)
        return np.r_[ray_rhs, np.dot(state[3:6], ray_rhs[:3])]

    def plane_distance(state):
        return float(np.dot(state[:3] - plane_point, normal))

    try:
        solver = DOP853(rhs, 0.0, initial, length, rtol=float(settings["rtol"]),
                        atol=atol, first_step=min(float(settings["ds_min_km"]), length, max_step),
                        max_step=max_step)
    except _ProviderDomainBoundary as exc:
        raise ValueError("initial AH port lacks a valid Ne derivative stencil") from exc
    states = [initial.copy()]
    paths = [0.0]
    previous_distance = plane_distance(initial)
    termination = "max_path_reached"
    max_residual = hamiltonian_residual_ah3d_ecef(initial[:7], model, t=time_s)
    boundary_retries = 0
    minimum_step = max(1e-10, float(settings["atol_r_km"]) * 0.1)

    while solver.status == "running":
        previous_path = float(solver.t)
        previous_state = solver.y.copy()
        attempted_step = min(float(solver.h_abs), float(solver.max_step), length - previous_path)
        try:
            solver.step()
            if solver.status == "failed" or not np.all(np.isfinite(solver.y)):
                raise RuntimeError("AH integration failed before reaching the plane")
            dense = solver.dense_output()
            projected, _, ok = project_dispersion_ah3d_ecef(solver.y[:7], model, t=time_s)
            if not ok:
                raise RuntimeError("AH dispersion projection failed")
            solver.y[:7] = projected
            solver.f = rhs(solver.t, solver.y)
        except _ProviderDomainBoundary as exc:
            boundary_retries += 1
            max_step = min(max_step, attempted_step / 2)
            if boundary_retries > 30 or max_step < minimum_step:
                raise RuntimeError(
                    "AH integration cannot advance within the Ne provider derivative domain") from exc
            solver = DOP853(rhs, previous_path, previous_state, length,
                            rtol=float(settings["rtol"]), atol=atol,
                            first_step=min(max_step, length - previous_path), max_step=max_step)
            continue
        current_distance = plane_distance(solver.y)
        positive = previous_distance < 0 <= current_distance
        negative = previous_distance > 0 >= current_distance
        crossing = (direction >= 0 and positive) or (direction <= 0 and negative)
        if crossing:
            event_path = brentq(lambda distance: plane_distance(dense(distance)),
                                previous_path, float(solver.t), xtol=1e-12)
            event_state = dense(event_path)
            event_state[:7], _, ok = project_dispersion_ah3d_ecef(event_state[:7], model, t=time_s)
            if not ok:
                raise RuntimeError("AH dispersion projection failed at the receiver plane")
            if model.provider_domain_margin(event_state[:3]) < 0:
                raise RuntimeError("receiver crossing lies outside the Ne provider domain")
            states.append(event_state.copy())
            paths.append(event_path)
            max_residual = max(max_residual, hamiltonian_residual_ah3d_ecef(
                event_state[:7], model, t=time_s))
            termination = "receiver_plane_crossed"
            break
        if model.provider_domain_margin(solver.y[:3]) < 0:
            raise RuntimeError("AH trajectory left the Ne provider domain")
        states.append(solver.y.copy())
        paths.append(float(solver.t))
        max_residual = max(max_residual, hamiltonian_residual_ah3d_ecef(
            solver.y[:7], model, t=time_s))
        previous_distance = current_distance

    # shortcut: action omits turning and polarization phases; add field transport before coherent use.
    trajectory = np.asarray(states)
    return {
        "termination_status": termination,
        "path_km": np.asarray(paths),
        "trajectory_ecef_km": trajectory[:, :3],
        "wavevector_ecef_rad_km": trajectory[:, 3:6],
        "group_time_s": trajectory[:, 6],
        "action_rad": trajectory[:, 7],
        "endpoint_ecef_km": trajectory[-1, :3].copy(),
        "endpoint_k_ecef_rad_km": trajectory[-1, 3:6].copy(),
        "endpoint_k_normal_rad_km": float(np.dot(trajectory[-1, 3:6], normal)),
        "terminal_plane_distance_km": plane_distance(trajectory[-1]),
        "max_normalized_dispersion_residual": max_residual,
        "medium_time_s": time_s,
        "provider_boundary_retries": boundary_retries,
    }
