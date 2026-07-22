"""Local ray-frame utilities for update2.1.

The first Ray-PE/SSF coupling pass uses a 2.5-D ``s-u`` cross-section around a
3-D ray.  This module only builds geometry and samples electron density; it does
not modify propagation power, delay, or spectra.
"""

from __future__ import annotations

import numpy as np

from utils import geocentric_to_ecef


def _trajectory_positions(ray: dict) -> np.ndarray:
    traj = ray.get("trajectory_ecef_km")
    if traj is None or len(traj) < 2:
        raise ValueError(
            "formal local frame requires trajectory_ecef_km with at least two states")
    arr = np.asarray(traj, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError("trajectory_ecef_km must contain three position columns")
    return arr[:, :3]


def _frame_origin(frame: dict) -> np.ndarray:
    return np.asarray(frame["origin_ecef_km"], dtype=float)


def _unit(v: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < 1e-12:
        return fallback.astype(float) / max(float(np.linalg.norm(fallback)), 1e-30)
    return v.astype(float) / n


def make_ray_local_frame(ray: dict, sample_index: int | None = None) -> dict:
    """Return a right-handed local frame attached to a 3-D ray.

    The frame axes are in km-coordinate space:
    ``s_hat`` follows the ray tangent, ``v_hat`` is approximately horizontal
    cross-track, and ``u_hat = v_hat x s_hat`` completes the local propagation
    plane.  For near-vertical rays the fallback basis remains well-defined.
    """
    pts = _trajectory_positions(ray)
    n_pts = pts.shape[0]
    if sample_index is None:
        sample_index = n_pts // 2
    sample_index = int(np.clip(sample_index, 0, n_pts - 1))

    stored_tangents = np.asarray(
        ray.get("trajectory_group_tangent_ecef", []), dtype=float)
    if (stored_tangents.shape == pts.shape
            and np.all(np.isfinite(stored_tangents))):
        tangent = stored_tangents[sample_index]
    elif sample_index == 0:
        tangent = pts[1] - pts[0]
    elif sample_index == n_pts - 1:
        tangent = pts[-1] - pts[-2]
    else:
        tangent = pts[sample_index + 1] - pts[sample_index - 1]

    s_hat = _unit(tangent, np.array([1.0, 0.0, 0.0]))
    origin = pts[sample_index].astype(float)
    radial_up = _unit(origin, np.array([0.0, 0.0, 1.0]))
    u_hat = radial_up - np.dot(radial_up, s_hat) * s_hat
    if np.linalg.norm(u_hat) < 1e-10:
        u_hat = np.cross(np.array([0.0, 1.0, 0.0]), s_hat)
    u_hat = _unit(u_hat, np.array([0.0, 0.0, 1.0]))
    v_hat = np.cross(s_hat, u_hat)
    if np.linalg.norm(v_hat) < 1e-10:
        v_hat = np.cross(s_hat, np.array([0.0, 1.0, 0.0]))
    v_hat = _unit(v_hat, np.array([0.0, 1.0, 0.0]))
    u_hat = _unit(np.cross(v_hat, s_hat), radial_up)
    v_hat = _unit(np.cross(s_hat, u_hat), v_hat)

    frame = {
        "s_hat": s_hat,
        "u_hat": u_hat,
        "v_hat": v_hat,
        "radial_up_ecef": radial_up,
        "sample_index": sample_index,
        "basis_matrix": np.vstack([s_hat, u_hat, v_hat]),
        "frame_coordinate_model": "ECEF_ray_tangent_radial_up",
        "origin_ecef_km": origin,
    }
    return frame


def estimate_ray_curvature_km_inv(ray: dict, sample_index: int | None = None) -> float:
    """Estimate local ray curvature |dt/ds| in km^-1 from trajectory points."""
    pts = _trajectory_positions(ray)
    n_pts = pts.shape[0]
    if n_pts < 3:
        return 0.0
    if sample_index is None:
        sample_index = n_pts // 2
    i = int(np.clip(sample_index, 1, n_pts - 2))
    ds1 = max(float(np.linalg.norm(pts[i] - pts[i - 1])), 1e-12)
    ds2 = max(float(np.linalg.norm(pts[i + 1] - pts[i])), 1e-12)
    t1 = _unit(pts[i] - pts[i - 1], np.array([1.0, 0.0, 0.0]))
    t2 = _unit(pts[i + 1] - pts[i], t1)
    return float(np.linalg.norm(t2 - t1) / (0.5 * (ds1 + ds2)))


def make_boundary_field_descriptor(ray: dict,
                                   frame: dict | None = None,
                                   sample_index: int | None = None,
                                   Pr_W: float | None = None,
                                   k0_per_km: float | None = None,
                                   n3d=None,
                                   t: float = 0.0) -> dict:
    """Describe the local ray field used to launch a PE/SSF section.

    The descriptor carries geometry and phase information only. Callers decide
    whether the PE result is diagnostic or applied to power.
    """
    pts = _trajectory_positions(ray)
    if frame is None:
        frame = make_ray_local_frame(ray, sample_index=sample_index)
    idx = int(frame.get("sample_index", pts.shape[0] // 2))
    if idx == 0:
        tangent = pts[1] - pts[0]
    elif idx >= pts.shape[0] - 1:
        tangent = pts[-1] - pts[-2]
    else:
        tangent = pts[idx + 1] - pts[idx - 1]
    origin = _frame_origin(frame)
    k_hat = _unit(tangent, np.array([1.0, 0.0, 0.0]))
    up = np.asarray(frame["radial_up_ecef"], dtype=float)
    beta = float(np.degrees(np.arcsin(np.clip(np.dot(k_hat, up), -1.0, 1.0))))
    amp = float(np.sqrt(max(float(Pr_W if Pr_W is not None else ray.get("Pr_W", 0.0)), 0.0)))
    s_ref = 0.0
    geom = ray.get("geom_path_km")
    if geom is not None:
        s_ref = float(geom)
    elif idx > 0:
        s_ref = float(np.sum(np.linalg.norm(np.diff(pts[:idx + 1], axis=0), axis=1)))
    phase = float((k0_per_km or 0.0) * s_ref)
    grad_local = np.zeros(3, dtype=float)
    if n3d is not None:
        h = 1.0
        for axis, key in enumerate(("s_hat", "u_hat", "v_hat")):
            direction = np.asarray(frame[key], dtype=float)
            rp = origin + h * direction
            rm = origin - h * direction
            if not hasattr(n3d, "get_Ne_ecef"):
                raise TypeError("formal local frame sampling requires get_Ne_ecef")
            radius = float(ray["earth_radius_km"])
            nep = float(n3d.get_Ne_ecef(rp, radius, t=t))
            nem = float(n3d.get_Ne_ecef(rm, radius, t=t))
            grad_local[axis] = (nep - nem) / (2.0 * h)
    return {
        "amplitude_sqrt_W": amp,
        "phase_rad": phase,
        "k_hat": k_hat,
        "beta_inc_deg": beta,
        "curvature_km_inv": estimate_ray_curvature_km_inv(ray, idx),
        "polarization": str(ray.get("wave_mode", ray.get("polarization", "scalar"))),
        "sample_index": idx,
        "Ne_gradient_s_m3_per_km": float(grad_local[0]),
        "Ne_gradient_u_m3_per_km": float(grad_local[1]),
        "Ne_gradient_v_m3_per_km": float(grad_local[2]),
        "boundary_field_model": "ray_tangent_scalar_envelope",
        "boundary_coordinate_model": frame["frame_coordinate_model"],
        "theory_reference": "Knepp 1983 MPS; Lehtinen and Inan 2008 full-wave/FDTD interface concept",
    }


def local_to_global(frame: dict, s_km: float = 0.0, u_km: float = 0.0, v_km: float = 0.0) -> np.ndarray:
    """Map local ``(s,u,v)`` offsets in km to global ``(x,y,z)`` km."""
    return (
        _frame_origin(frame)
        + float(s_km) * np.asarray(frame["s_hat"], dtype=float)
        + float(u_km) * np.asarray(frame["u_hat"], dtype=float)
        + float(v_km) * np.asarray(frame["v_hat"], dtype=float)
    )


def sample_Ne_local_cross_section(
    n3d,
    ray: dict,
    s_offsets_km=None,
    u_offsets_km=None,
    sample_index: int | None = None,
    t: float = 0.0,
) -> dict:
    """Sample ``Ne`` on a local 2.5-D ``s-u`` cross-section.

    Returns a dictionary containing ``s_km``, ``u_km``, ``Ne_m3``, global
    coordinate grids ``x_km/y_km/z_km``, and the ``frame`` used for sampling.
    """
    if s_offsets_km is None:
        s_offsets_km = np.linspace(-50.0, 50.0, 51)
    if u_offsets_km is None:
        u_offsets_km = np.linspace(-25.0, 25.0, 51)
    s_arr = np.asarray(s_offsets_km, dtype=float)
    u_arr = np.asarray(u_offsets_km, dtype=float)
    frame = make_ray_local_frame(ray, sample_index=sample_index)

    Ne = np.empty((s_arr.size, u_arr.size), dtype=float)
    x_grid = np.empty_like(Ne)
    y_grid = np.empty_like(Ne)
    z_grid = np.empty_like(Ne)
    for i, s_val in enumerate(s_arr):
        for j, u_val in enumerate(u_arr):
            x_val, y_val, z_val = local_to_global(frame, s_val, u_val, 0.0)
            x_grid[i, j] = x_val
            y_grid[i, j] = y_val
            z_grid[i, j] = z_val
            if not hasattr(n3d, "get_Ne_ecef"):
                raise TypeError("formal local frame sampling requires get_Ne_ecef")
            Ne[i, j] = max(float(n3d.get_Ne_ecef(
                (x_val, y_val, z_val), float(ray["earth_radius_km"]), t=t)), 0.0)

    result = {
        "s_km": s_arr,
        "u_km": u_arr,
        "Ne_m3": Ne,
        "frame": frame,
        "coordinate_model": frame["frame_coordinate_model"],
    }
    result.update({"ecef_x_km": x_grid, "ecef_y_km": y_grid,
                   "ecef_z_km": z_grid})
    return result



def extract_irregularity_entry(ray: dict,
                               irregularity_params: dict,
                               q_boundary: float = 1.0) -> dict | None:
    """Find first trajectory point inside a field-aligned ellipsoidal irregularity.

    The boundary is q <= q_boundary, with q = |r_perp|^2/sigma_perp^2 +
    r_parallel^2/sigma_para^2. This matches the Gaussian field-aligned bubble
    geometry used by RefractiveIndex3D; it returns geometry diagnostics only.
    """
    pts = _trajectory_positions(ray)
    if not irregularity_params or not irregularity_params.get('enable', False):
        return None
    if "center_ecef_km" in irregularity_params:
        center = np.asarray(irregularity_params["center_ecef_km"], dtype=float)
    elif "center_geocentric_lat_lon_height" in irregularity_params:
        center = geocentric_to_ecef(
            *irregularity_params["center_geocentric_lat_lon_height"],
            float(ray["earth_radius_km"]))
    else:
        raise ValueError(
            "ECEF irregularity entry requires center_ecef_km or "
            "center_geocentric_lat_lon_height")
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise ValueError("irregularity ECEF center must be finite length-3")
    s_perp = max(float(irregularity_params.get('sigma_perp_km', irregularity_params.get('Lx_km', 50.0))), 1e-9)
    s_para = max(float(irregularity_params.get('sigma_para_km', irregularity_params.get('Lz_km', 100.0))), 1e-9)
    b_hat = np.asarray(irregularity_params.get("b_ecef_hat", []), dtype=float)
    if b_hat.shape != (3,) or not np.all(np.isfinite(b_hat)):
        raise ValueError("ECEF irregularity entry requires finite b_ecef_hat")
    b_hat = _unit(b_hat, np.array([0.0, 0.0, 1.0]))

    q_vals = []
    for pt in pts:
        r = pt - center
        r_para = float(np.dot(r, b_hat))
        r_perp = r - r_para * b_hat
        q_vals.append(float(np.dot(r_perp, r_perp) / (s_perp * s_perp) + (r_para * r_para) / (s_para * s_para)))
    q_vals = np.asarray(q_vals, dtype=float)
    inside = np.where(q_vals <= float(q_boundary))[0]
    if inside.size == 0:
        return None
    idx = int(inside[0])
    if idx == 0:
        tangent = pts[1] - pts[0]
    else:
        tangent = pts[idx] - pts[idx - 1]
    up = _unit(pts[idx], np.array([0.0, 0.0, 1.0]))
    beta = float(np.degrees(np.arcsin(np.clip(
        np.dot(_unit(tangent, up), up), -1.0, 1.0))))
    frame = make_ray_local_frame(ray, sample_index=idx)
    descriptor = make_boundary_field_descriptor(ray, frame=frame, sample_index=idx,
                                                Pr_W=ray.get('Pr_W'), k0_per_km=ray.get('k0_per_km'))
    result = {
        'index': idx,
        'q': float(q_vals[idx]),
        'beta_inc_deg': beta,
        'frame': frame,
        'boundary_field': descriptor,
        'coordinate_model': frame['frame_coordinate_model'],
    }
    result['entry_ecef_km'] = pts[idx].tolist()
    return result


def make_fdtd_adapter_descriptor(ray: dict,
                                 frame: dict | None = None,
                                 sample_index: int | None = None,
                                 domain_half_width_km: float = 50.0,
                                 domain_length_km: float = 100.0,
                                 grid_step_km: float = 1.0) -> dict:
    """Describe, but do not run, a local FDTD/FARR coupling domain."""
    if frame is None:
        frame = make_ray_local_frame(ray, sample_index=sample_index)
    boundary = make_boundary_field_descriptor(ray, frame=frame, sample_index=frame["sample_index"])
    return {
        "fdtd_adapter_status": "descriptor_only",
        "fdtd_solver": "farr_adapter_placeholder",
        "solver_dimensionality": "fullwave_domain_descriptor",
        "validity": "interface_only_not_executed",
        "fdtd_domain": {
            "origin_ecef_km": _frame_origin(frame).tolist(),
            "s_length_km": float(domain_length_km),
            "u_half_width_km": float(domain_half_width_km),
            "v_half_width_km": float(domain_half_width_km),
            "grid_step_km": float(grid_step_km),
        },
        "boundary_field": boundary,
        "receiver_surface": {
            "normal_axis": "s",
            "offset_km": float(0.5 * domain_length_km),
        },
        "field_to_pd": {
            "method": "field_time_series_to_delay_power",
            "implemented": False,
        },
        "theory_reference": "Yee 1966 FDTD; Lehtinen and Inan 2008 full-wave method; FARR local full-wave reference",
    }
