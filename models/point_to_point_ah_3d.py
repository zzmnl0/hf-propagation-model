"""Point-to-point solver for 3-D Appleton-Hartree O/X rays."""

from __future__ import annotations

import numpy as np

from config import FREQ_MHZ, P2P_3D, RT
from .point_to_point_3d import classify_mode_3d, annotate_mode_causes_3d, _seed_betas_3d
from .ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


def _shoot_ah_3d(
    tx: tuple,
    rx: tuple,
    n_ah: RefractiveIndexAH3D,
    freq_MHz: float,
    beta_init: float,
    az_init: float = 0.0,
    max_iter: int = 12,
    tol_km: float = 5.0,
    d_beta: float = 0.05,
    d_az: float = 0.10,
    jacobian: str = "forward",
    t: float = 0.0,
):
    beta = float(beta_init)
    az = float(az_init)
    x_rx = float(rx[0])
    y_rx = 0.0
    jacobian = str(jacobian).lower()
    if jacobian not in ("forward", "central"):
        jacobian = "forward"
    prev_err = float("inf")

    for k in range(max_iter):
        ray_c = trace_single_ray_ah_3d(tx, beta, az, n_ah, freq_MHz, t=t)
        x_c = ray_c["x_land_km"]
        y_c = ray_c["y_land_km"]

        if jacobian == "central":
            ray_bp = trace_single_ray_ah_3d(tx, beta + d_beta, az, n_ah, freq_MHz, t=t)
            ray_bm = trace_single_ray_ah_3d(tx, beta - d_beta, az, n_ah, freq_MHz, t=t)
            ray_ap = trace_single_ray_ah_3d(tx, beta, az + d_az, n_ah, freq_MHz, t=t)
            ray_am = trace_single_ray_ah_3d(tx, beta, az - d_az, n_ah, freq_MHz, t=t)
            dx_db = (ray_bp["x_land_km"] - ray_bm["x_land_km"]) / (2.0 * d_beta)
            dy_db = (ray_bp["y_land_km"] - ray_bm["y_land_km"]) / (2.0 * d_beta)
            dx_da = (ray_ap["x_land_km"] - ray_am["x_land_km"]) / (2.0 * d_az)
            dy_da = (ray_ap["y_land_km"] - ray_am["y_land_km"]) / (2.0 * d_az)
        else:
            ray_b = trace_single_ray_ah_3d(tx, beta + d_beta, az, n_ah, freq_MHz, t=t)
            ray_a = trace_single_ray_ah_3d(tx, beta, az + d_az, n_ah, freq_MHz, t=t)
            dx_db = (ray_b["x_land_km"] - x_c) / d_beta
            dy_db = (ray_b["y_land_km"] - y_c) / d_beta
            dx_da = (ray_a["x_land_km"] - x_c) / d_az
            dy_da = (ray_a["y_land_km"] - y_c) / d_az

        J = np.array([[dx_db, dx_da], [dy_db, dy_da]], dtype=float)
        F = np.array([x_rx - x_c, y_rx - y_c], dtype=float)
        err = float(np.linalg.norm(F))
        if err < tol_km:
            ray_c["beta_deg_tx"] = beta
            ray_c["az_deg_tx"] = az
            ray_c["elevation_tx_deg"] = beta
            ray_c["elevation_rx_deg"] = ray_c["beta_deg_rx"]
            ray_c["label"] = classify_mode_3d(ray_c)
            return ray_c

        if k >= 3 and err > prev_err:
            break
        prev_err = err

        det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
        if abs(det) < 1e-12:
            break
        J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]], dtype=float) / det
        delta = J_inv @ F
        beta += float(delta[0])
        az += float(delta[1])
        beta = float(np.clip(beta, 2.0, 88.0))
    return None


def _deduplicate(candidates: list, clust_h_km: float, clust_tau_ms: float) -> list:
    unique = []
    for c in candidates:
        is_dup = any(
            abs(c["h_reflect_km"] - u["h_reflect_km"]) < clust_h_km
            and abs(c["tau_ms"] - u["tau_ms"]) < clust_tau_ms
            and c.get("wave_mode") == u.get("wave_mode")
            for u in unique
        )
        if not is_dup:
            unique.append(c)
    return unique


def find_all_rays_ah_3d(
    tx: tuple,
    rx: tuple,
    n3d,
    freq_MHz: float = FREQ_MHZ,
    p2p_params: dict = P2P_3D,
    t: float = 0.0,
    wave_mode: str = "both",
    geomag: dict | None = None,
    rays_fan: list | None = None,
) -> list:
    """Find AH3D P2P rays for O, X, or both branches."""
    n_init = int(p2p_params.get("n_init", 6))
    clust_h = float(p2p_params.get("clust_h_km", 10.0))
    clust_tau = float(p2p_params.get("clust_tau_ms", 0.05))
    tol_km = float(p2p_params.get("tol_km", 5.0))
    max_iter = int(p2p_params.get("max_iter", 12))
    d_beta = float(p2p_params.get("d_beta", 0.05))
    d_az = float(p2p_params.get("d_az", 0.10))
    jacobian = p2p_params.get("jacobian", "forward")

    wm = str(wave_mode).upper()
    modes = ["O", "X"] if wm == "BOTH" else [wm]
    betas, seed_source = _seed_betas_3d(rx, p2p_params, rays_fan)

    candidates = []
    for mode in modes:
        n_ah = RefractiveIndexAH3D(n3d, freq_MHz, wave_mode=mode, geomag=geomag)
        for beta in betas:
            ray = _shoot_ah_3d(
                tx,
                rx,
                n_ah,
                freq_MHz,
                beta,
                0.0,
                max_iter=max_iter,
                tol_km=tol_km,
                d_beta=d_beta,
                d_az=d_az,
                jacobian=jacobian,
                t=t,
            )
            if ray is None:
                continue
            if abs(ray["x_land_km"] - float(rx[0])) < tol_km and abs(ray["y_land_km"]) < tol_km:
                if ray["h_reflect_km"] < RT["z_stop_km"] and ray["group_path_km"] > 0.1:
                    ray["label"] = classify_mode_3d(ray)
                    ray["p2p_seed_source"] = seed_source
                    candidates.append(ray)

    unique = _deduplicate(candidates, clust_h, clust_tau)
    unique.sort(key=lambda r: r["tau_ms"])
    annotate_mode_causes_3d(unique, p2p_params, seed_source)
    return unique
