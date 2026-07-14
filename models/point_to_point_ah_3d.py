"""Point-to-point solver for 3-D Appleton-Hartree O/X rays."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

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
    jacobian: str = "variational",
    t: float = 0.0,
    formal_profile: str = "research",
    jacobian_rtol: float = 1e-3,
):
    x_rx = float(rx[0])
    y_rx = 0.0
    jacobian = str(jacobian).lower()
    if jacobian not in ("variational", "central", "forward"):
        jacobian = "variational"
    cache = {}
    fallback_count = 0
    variational_checked = False
    variational_valid = True
    variational_error = np.nan

    def trace(theta, sensitivity=False):
        key = (float(theta[0]), float(theta[1]), bool(sensitivity))
        if key not in cache:
            cache[key] = trace_single_ray_ah_3d(
                tx, key[0], key[1], n_ah, freq_MHz, t=t,
                with_sensitivity=sensitivity)
        return cache[key]

    def residual(theta):
        ray = trace(theta, sensitivity=(jacobian == "variational"))
        return np.array([ray["x_land_km"] - x_rx,
                         ray["y_land_km"] - y_rx])

    def central_jac(theta):
        bp = trace((theta[0] + d_beta, theta[1]))
        bm = trace((theta[0] - d_beta, theta[1]))
        ap = trace((theta[0], theta[1] + d_az))
        am = trace((theta[0], theta[1] - d_az))
        return np.array([
            [(bp["x_land_km"] - bm["x_land_km"]) / (2.0 * d_beta),
             (ap["x_land_km"] - am["x_land_km"]) / (2.0 * d_az)],
            [(bp["y_land_km"] - bm["y_land_km"]) / (2.0 * d_beta),
             (ap["y_land_km"] - am["y_land_km"]) / (2.0 * d_az)],
        ])

    def jac(theta):
        nonlocal fallback_count, variational_checked, variational_valid, variational_error
        if jacobian == "forward":
            center = trace(theta)
            bp = trace((theta[0] + d_beta, theta[1]))
            ap = trace((theta[0], theta[1] + d_az))
            return np.array([
                [(bp["x_land_km"] - center["x_land_km"]) / d_beta,
                 (ap["x_land_km"] - center["x_land_km"]) / d_az],
                [(bp["y_land_km"] - center["y_land_km"]) / d_beta,
                 (ap["y_land_km"] - center["y_land_km"]) / d_az],
            ])
        if jacobian == "central":
            return central_jac(theta)
        if variational_checked and not variational_valid:
            return central_jac(theta)
        ray = trace(theta, sensitivity=True)
        value = np.asarray(ray.get("landing_jacobian_deg"), dtype=float)
        condition = float(ray.get("landing_jacobian_condition", np.inf))
        if value.shape != (2, 2) or not np.all(np.isfinite(value)) or condition > 1e10:
            fallback_count += 1
            variational_checked = True
            variational_valid = False
            return central_jac(theta)
        if not variational_checked:
            reference = central_jac(theta)
            variational_error = float(np.linalg.norm(value - reference) /
                                      max(np.linalg.norm(reference), 1e-12))
            variational_checked = True
            variational_valid = variational_error <= jacobian_rtol
            if not variational_valid:
                fallback_count += 1
                return reference
        return value

    try:
        solved = least_squares(
            residual, np.array([beta_init, az_init], dtype=float), jac=jac,
            bounds=([2.0, -180.0], [88.0, 180.0]), method="trf",
            x_scale="jac", max_nfev=max_iter, ftol=1e-9, xtol=1e-9, gtol=1e-9)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return None
    if float(np.linalg.norm(solved.fun)) >= tol_km:
        return None

    beta, az = map(float, solved.x)
    previous_profile = n_ah.integration_profile
    n_ah.integration_profile = formal_profile
    try:
        ray_final = trace_single_ray_ah_3d(tx, beta, az, n_ah, freq_MHz, t=t)
    finally:
        n_ah.integration_profile = previous_profile
    landing_error = float(np.hypot(ray_final["x_land_km"] - x_rx,
                                   ray_final["y_land_km"] - y_rx))
    if ray_final.get("convergence_status") != "passed" or landing_error >= tol_km:
        return None
    ray_final.update({
        "beta_deg_tx": beta, "az_deg_tx": az,
        "elevation_tx_deg": beta,
        "elevation_rx_deg": ray_final["beta_deg_rx"],
        "label": classify_mode_3d(ray_final),
        "p2p_search_profile": "smoke",
        "p2p_final_profile": formal_profile,
        "p2p_solver": "scipy_trf",
        "p2p_jacobian_method": jacobian,
        "p2p_jacobian_fallback_count": fallback_count,
        "p2p_jacobian_validation_relative": float(variational_error),
        "p2p_jacobian_effective": ("variational" if variational_valid
                                    else "central_reference_fallback"),
        "p2p_landing_error_km": landing_error,
        "p2p_nfev": int(solved.nfev),
        "p2p_njev": int(solved.njev or 0),
    })
    return ray_final


def _candidate_seeds(tx, rx, n_ah, freq_MHz, betas, rays_fan, t, nearest=3):
    """Find full-domain crossing/fold candidates from a coarse landing map."""
    branch = n_ah.wave_mode
    fan = [r for r in (rays_fan or []) if r.get("wave_mode") == branch]
    if len(fan) < 3:
        fan = [trace_single_ray_ah_3d(tx, beta, 0.0, n_ah, freq_MHz, t=t)
               for beta in betas]
    fan = sorted((r for r in fan if r.get("branch_status") == "ground"),
                 key=lambda r: float(r.get("beta_deg_tx", 0.0)))
    if not fan:
        return [], 0
    target = float(rx[0])
    errors = np.array([float(r["x_land_km"]) - target for r in fan])
    seeds = []
    for left, right, e_left, e_right in zip(fan[:-1], fan[1:], errors[:-1], errors[1:]):
        if e_left == 0.0 or e_left * e_right < 0.0:
            seeds.append((0.5 * (float(left["beta_deg_tx"]) +
                                 float(right["beta_deg_tx"])), 0.0))
    for i in range(1, len(fan) - 1):
        if abs(errors[i]) <= abs(errors[i - 1]) and abs(errors[i]) <= abs(errors[i + 1]):
            seeds.append((float(fan[i]["beta_deg_tx"]),
                          float(fan[i].get("az_deg_tx", 0.0))))
    order = np.argsort(np.abs(errors))[:min(nearest, len(fan))]
    seeds.extend((float(fan[i]["beta_deg_tx"]),
                  float(fan[i].get("az_deg_tx", 0.0))) for i in order)
    unique = []
    for seed in seeds:
        if not any(abs(seed[0] - old[0]) < 0.02 and abs(seed[1] - old[1]) < 0.02
                   for old in unique):
            unique.append(seed)
    return unique, len(fan)


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
    continuation_modes: list | None = None,
) -> list:
    """Find AH3D P2P rays for O, X, or both branches."""
    n_init = int(p2p_params.get("n_init", 6))
    clust_h = float(p2p_params.get("clust_h_km", 10.0))
    clust_tau = float(p2p_params.get("clust_tau_ms", 0.05))
    tol_km = float(p2p_params.get("tol_km", 5.0))
    max_iter = int(p2p_params.get("max_iter", 12))
    d_beta = float(p2p_params.get("d_beta", 0.05))
    d_az = float(p2p_params.get("d_az", 0.10))
    jacobian = p2p_params.get("jacobian", "variational")
    jacobian_rtol = float(p2p_params.get("jacobian_rtol", 1e-3))

    wm = str(wave_mode).upper()
    modes = ["O", "X"] if wm == "BOTH" else [wm]
    betas, seed_source = _seed_betas_3d(rx, p2p_params, None)

    candidates = []
    rejected_status = {}
    for mode in modes:
        n_ah = RefractiveIndexAH3D(
            n3d, freq_MHz, wave_mode=mode, geomag=geomag,
            integration_profile="smoke",
        )
        seeds, audit_count = _candidate_seeds(
            tx, rx, n_ah, freq_MHz, betas, rays_fan, t,
            nearest=int(p2p_params.get("candidate_nearest", 3)))
        continued = [
            (float(r.get("beta_deg_tx", r.get("elevation_tx_deg", 0.0))),
             float(r.get("az_deg_tx", 0.0)))
            for r in (continuation_modes or []) if r.get("wave_mode") == mode
        ]
        for seed in continued + seeds:
            beta, az = seed
            ray = _shoot_ah_3d(
                tx,
                rx,
                n_ah,
                freq_MHz,
                beta,
                az,
                max_iter=max_iter,
                tol_km=tol_km,
                d_beta=d_beta,
                d_az=d_az,
                jacobian=jacobian,
                jacobian_rtol=jacobian_rtol,
                t=t,
            )
            if ray is None:
                continue
            status = ray.get("branch_status", "ok")
            if status != "ground" or ray.get("convergence_status") != "passed":
                rejected_status[status] = rejected_status.get(status, 0) + 1
                continue
            if abs(ray["x_land_km"] - float(rx[0])) < tol_km and abs(ray["y_land_km"]) < tol_km:
                if ray["h_reflect_km"] < RT["z_stop_km"] and ray["group_path_km"] > 0.1:
                    ray["label"] = classify_mode_3d(ray)
                    ray["p2p_seed_source"] = seed_source
                    ray["freq_MHz"] = float(freq_MHz)
                    ray["ox_solver"] = "ah3d_hamiltonian"
                    ray["p2p_candidate_audit_rays"] = audit_count
                    ray["p2p_full_domain_audit"] = True
                    ray["p2p_continuation_used"] = seed in continued
                    candidates.append(ray)

    unique = _deduplicate(candidates, clust_h, clust_tau)
    unique.sort(key=lambda r: r["tau_ms"])
    annotate_mode_causes_3d(unique, p2p_params, seed_source)
    for ray in unique:
        ray["ah3d_rejected_status_counts"] = dict(rejected_status)
    return unique
