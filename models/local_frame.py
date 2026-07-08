"""Local ray-frame utilities for update2.1.

The first Ray-PE/SSF coupling pass uses a 2.5-D ``s-u`` cross-section around a
3-D ray.  This module only builds geometry and samples electron density; it does
not modify propagation power, delay, or spectra.
"""

from __future__ import annotations

import numpy as np


def _trajectory_positions(ray: dict) -> np.ndarray:
    traj = ray.get("trajectory_3d")
    if traj is None or len(traj) < 2:
        raise ValueError("ray must contain trajectory_3d with at least two states")
    arr = np.asarray(traj, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError("trajectory_3d must be an array-like sequence with x,y,z columns")
    return arr[:, :3]


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

    if sample_index == 0:
        tangent = pts[1] - pts[0]
    elif sample_index == n_pts - 1:
        tangent = pts[-1] - pts[-2]
    else:
        tangent = pts[sample_index + 1] - pts[sample_index - 1]

    s_hat = _unit(tangent, np.array([1.0, 0.0, 0.0]))
    z_hat = np.array([0.0, 0.0, 1.0])
    v_hat = np.cross(s_hat, z_hat)
    if np.linalg.norm(v_hat) < 1e-10:
        v_hat = np.cross(s_hat, np.array([0.0, 1.0, 0.0]))
    v_hat = _unit(v_hat, np.array([0.0, 1.0, 0.0]))
    u_hat = _unit(np.cross(v_hat, s_hat), z_hat)
    v_hat = _unit(np.cross(s_hat, u_hat), v_hat)

    origin = pts[sample_index].astype(float)
    return {
        "origin_km": origin,
        "s_hat": s_hat,
        "u_hat": u_hat,
        "v_hat": v_hat,
        "sample_index": sample_index,
        "basis_matrix": np.vstack([s_hat, u_hat, v_hat]),
    }


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
                                   k0_per_km: float | None = None) -> dict:
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
    k_hat = _unit(tangent, np.array([1.0, 0.0, 0.0]))
    beta = float(np.degrees(np.arctan2(k_hat[2], max(np.hypot(k_hat[0], k_hat[1]), 1e-12))))
    amp = float(np.sqrt(max(float(Pr_W if Pr_W is not None else ray.get("Pr_W", 0.0)), 0.0)))
    s_ref = 0.0
    geom = ray.get("geom_path_km")
    if geom is not None:
        s_ref = float(geom)
    elif idx > 0:
        s_ref = float(np.sum(np.linalg.norm(np.diff(pts[:idx + 1], axis=0), axis=1)))
    phase = float((k0_per_km or 0.0) * s_ref)
    return {
        "amplitude_sqrt_W": amp,
        "phase_rad": phase,
        "k_hat": k_hat,
        "beta_inc_deg": beta,
        "curvature_km_inv": estimate_ray_curvature_km_inv(ray, idx),
        "polarization": str(ray.get("wave_mode", ray.get("polarization", "scalar"))),
        "sample_index": idx,
    }


def local_to_global(frame: dict, s_km: float = 0.0, u_km: float = 0.0, v_km: float = 0.0) -> np.ndarray:
    """Map local ``(s,u,v)`` offsets in km to global ``(x,y,z)`` km."""
    return (
        np.asarray(frame["origin_km"], dtype=float)
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
            Ne[i, j] = max(float(n3d.Ne_at(x_val, y_val, z_val, t=t)), 0.0)

    return {
        "s_km": s_arr,
        "u_km": u_arr,
        "Ne_m3": Ne,
        "x_km": x_grid,
        "y_km": y_grid,
        "z_km": z_grid,
        "frame": frame,
    }



def _geomag_b_hat(params: dict | None = None) -> np.ndarray:
    """Approximate link-coordinate geomagnetic unit vector."""
    try:
        from config import GEOMAG, LINK_BEARING_DEG
        g = GEOMAG if params is None else {**GEOMAG, **params}
        dip = np.radians(float(g.get('dip_deg', 48.7)))
        decl = np.radians(float(g.get('decl_deg', -5.5)))
        link = np.radians(float(g.get('link_bearing_deg', LINK_BEARING_DEG)))
        bx = np.cos(dip) * np.cos(decl - link)
        by = np.cos(dip) * np.sin(decl - link)
        bz = -np.sin(dip)
        return _unit(np.array([bx, by, bz], dtype=float), np.array([1.0, 0.0, 0.0]))
    except Exception:
        return np.array([1.0, 0.0, 0.0], dtype=float)


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
    center = np.array([
        float(irregularity_params.get('x0_km', irregularity_params.get('x_km', 0.0))),
        float(irregularity_params.get('y0_km', 0.0)),
        float(irregularity_params.get('z0_km', 300.0)),
    ], dtype=float)
    s_perp = max(float(irregularity_params.get('sigma_perp_km', irregularity_params.get('Lx_km', 50.0))), 1e-9)
    s_para = max(float(irregularity_params.get('sigma_para_km', irregularity_params.get('Lz_km', 100.0))), 1e-9)
    b_hat = _geomag_b_hat(irregularity_params)

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
    beta = float(np.degrees(np.arctan2(tangent[2], max(np.hypot(tangent[0], tangent[1]), 1e-12))))
    frame = make_ray_local_frame(ray, sample_index=idx)
    descriptor = make_boundary_field_descriptor(ray, frame=frame, sample_index=idx,
                                                Pr_W=ray.get('Pr_W'), k0_per_km=ray.get('k0_per_km'))
    return {
        'index': idx,
        'x_km': float(pts[idx, 0]),
        'y_km': float(pts[idx, 1]),
        'z_km': float(pts[idx, 2]),
        'q': float(q_vals[idx]),
        'beta_inc_deg': beta,
        'frame': frame,
        'boundary_field': descriptor,
    }
