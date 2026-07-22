"""Point-to-point solver for 3-D Appleton-Hartree O/X rays."""

from __future__ import annotations

import time

import numpy as np
from scipy.optimize import least_squares

from config import FREQ_MHZ, P2P_3D, RT
from .point_to_point_3d import classify_mode_3d, annotate_mode_causes_3d, _seed_betas_3d
from .ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


class _JacobianReferenceError(RuntimeError):
    """Reject a local root whose finite-difference endpoint map is unreliable."""


def _trace_cached(tx, n_ah, freq_MHz, beta, az, t, sensitivity,
                  profile, cache, stats):
    key = (float(beta), float(az), bool(sensitivity), str(profile), float(t))
    if key in cache:
        stats["cache_hits"] += 1
        stats["avoided_ray_calls"] += 1
        return cache[key]
    stats["cache_misses"] += 1
    if stats.get("print_progress"):
        print("  p0r_ray_call={} branch={} profile={} beta={:.4f} az={:.4f} sensitivity={}".format(
            stats["cache_misses"], n_ah.wave_mode, profile, key[0], key[1], key[2]),
            flush=True)
    previous_profile = n_ah.integration_profile
    n_ah.integration_profile = profile
    started = time.perf_counter()
    try:
        ray = trace_single_ray_ah_3d(
            tx, key[0], key[1], n_ah, freq_MHz, t=t,
            with_sensitivity=sensitivity)
    finally:
        n_ah.integration_profile = previous_profile
    stats["profile_wall_s"][profile] = (
        stats["profile_wall_s"].get(profile, 0.0)
        + time.perf_counter() - started)
    stats["profile_nfev"][profile] = (
        stats["profile_nfev"].get(profile, 0) + int(ray.get("nfev", 0)))
    stats["profile_accepted_steps"][profile] = (
        stats["profile_accepted_steps"].get(profile, 0)
        + int(ray.get("accepted_steps", 0)))
    if profile == "research" and ray.get("convergence_reference_profile") == "reference":
        stats["reference_ray_calls"] += 1
        stats["reference_nfev"] += int(ray.get("convergence_reference_nfev", 0))
        stats["reference_accepted_steps"] += int(
            ray.get("convergence_reference_accepted_steps", 0))
    cache[key] = ray
    return ray


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
    jacobian: str = "central_reference_multiscale",
    t: float = 0.0,
    formal_profile: str = "research",
    jacobian_rtol: float = 1e-3,
    jacobian_scales: tuple = (0.5, 1.0, 2.0),
    jacobian_min_scale: float = 0.125,
    jacobian_scale_rtol: float = 1e-2,
    jacobian_condition_max: float = 1e10,
    trace_cache: dict | None = None,
    trace_stats: dict | None = None,
    defer_formal: bool = False,
    seed_info: dict | None = None,
    angle_bounds: tuple | None = None,
):
    x_rx = float(rx[0])
    y_rx = 0.0
    jacobian = str(jacobian).lower()
    if jacobian == "central":
        jacobian = "central_reference_multiscale"
    if jacobian not in ("variational", "central_reference_multiscale", "forward"):
        jacobian = "central_reference_multiscale"
    scales = tuple(sorted(float(scale) for scale in jacobian_scales
                          if np.isfinite(scale) and float(scale) > 0.0))
    if len(scales) < 3:
        scales = (0.5, 1.0, 2.0)
    cache = {} if trace_cache is None else trace_cache
    stats = trace_stats if trace_stats is not None else {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {},
        "profile_accepted_steps": {}, "reference_ray_calls": 0,
        "reference_nfev": 0, "reference_accepted_steps": 0,
    }
    fallback_count = 0
    variational_checked = False
    variational_valid = True
    variational_error = np.nan
    variational_column_errors = [np.nan, np.nan]
    fallback_reasons = []
    reference_diagnostics = {"evaluations": []}

    def trace(theta, sensitivity=False):
        return _trace_cached(
            tx, n_ah, freq_MHz, theta[0], theta[1], t, sensitivity,
            "smoke", cache, stats)

    def residual(theta):
        ray = trace(theta, sensitivity=(jacobian == "variational"))
        return np.array([ray["x_land_km"] - x_rx,
                         ray["y_land_km"] - y_rx])

    def central_jac(theta, beta_step=d_beta, az_step=d_az):
        bp = trace((theta[0] + beta_step, theta[1]))
        bm = trace((theta[0] - beta_step, theta[1]))
        ap = trace((theta[0], theta[1] + az_step))
        am = trace((theta[0], theta[1] - az_step))
        rays = (bp, bm, ap, am)
        if any(ray.get("branch_status") != "ground" for ray in rays):
            raise _JacobianReferenceError("jacobian_perturbation_not_ground")
        if any(not np.all(np.isfinite(
                [ray.get("x_land_km", np.nan), ray.get("y_land_km", np.nan)]))
               for ray in rays):
            raise _JacobianReferenceError("jacobian_perturbation_nonfinite")
        return np.array([
            [(bp["x_land_km"] - bm["x_land_km"]) / (2.0 * beta_step),
             (ap["x_land_km"] - am["x_land_km"]) / (2.0 * az_step)],
            [(bp["y_land_km"] - bm["y_land_km"]) / (2.0 * beta_step),
             (ap["y_land_km"] - am["y_land_km"]) / (2.0 * az_step)],
        ])

    def reference_jacobian(theta):
        active_scales = list(scales)
        refs = [central_jac(theta, d_beta * scale, d_az * scale)
                for scale in active_scales]
        changes = [float(np.linalg.norm(refs[i] - refs[i - 1]) /
                         max(np.linalg.norm(refs[i]), 1e-12))
                   for i in range(1, len(refs))]
        while min(changes) > jacobian_scale_rtol:
            finer = 0.5 * active_scales[0]
            if finer < jacobian_min_scale - 1e-15:
                break
            active_scales.insert(0, finer)
            refs.insert(0, central_jac(
                theta, d_beta * finer, d_az * finer))
            changes.insert(0, float(np.linalg.norm(refs[1] - refs[0]) /
                                    max(np.linalg.norm(refs[1]), 1e-12)))
        index = int(np.argmin(changes))
        reference = refs[index]
        stable = bool(changes[index] <= jacobian_scale_rtol)
        condition = float(np.linalg.cond(reference))
        meta = {
            "scales": list(active_scales),
            "relative_changes": changes,
            "selected_scale": active_scales[index],
            "stable": stable,
            "scale_rtol": float(jacobian_scale_rtol),
            "condition": condition,
            "condition_max": float(jacobian_condition_max),
        }
        reference_diagnostics["evaluations"].append(meta)
        reference_diagnostics.update(meta)
        if not stable:
            raise _JacobianReferenceError("jacobian_reference_unstable")
        if (not np.all(np.isfinite(reference)) or not np.isfinite(condition)
                or condition > jacobian_condition_max):
            raise _JacobianReferenceError("jacobian_reference_ill_conditioned")
        return reference

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
        if jacobian == "central_reference_multiscale":
            return reference_jacobian(theta)
        if variational_checked and not variational_valid:
            scale = float(reference_diagnostics.get("selected_scale", 1.0))
            return central_jac(theta, d_beta * scale, d_az * scale)
        ray = trace(theta, sensitivity=True)
        value = np.asarray(ray.get("landing_jacobian_deg"), dtype=float)
        condition = float(ray.get("landing_jacobian_condition", np.inf))
        if value.shape != (2, 2) or not np.all(np.isfinite(value)) or condition > 1e10:
            fallback_count += 1
            variational_checked = True
            variational_valid = False
            fallback_reasons.append(
                "nonfinite_or_condition" if value.shape != (2, 2)
                or not np.all(np.isfinite(value)) else "condition_gt_1e10")
            reference = reference_jacobian(theta)
            return reference
        if not variational_checked:
            reference = reference_jacobian(theta)
            variational_column_errors[:] = [
                float(np.linalg.norm(value[:, j] - reference[:, j]) /
                      max(np.linalg.norm(reference[:, j]), 1e-12))
                for j in range(2)]
            variational_error = float(np.linalg.norm(value - reference) /
                                      max(np.linalg.norm(reference), 1e-12))
            variational_checked = True
            variational_valid = variational_error <= jacobian_rtol
            if not variational_valid:
                fallback_count += 1
                fallback_reasons.append("relative_error_gt_threshold")
                return reference
        return value

    bounds = angle_bounds or ((2.0, 88.0), (-180.0, 180.0))
    lower = np.array([float(bounds[0][0]), float(bounds[1][0])])
    upper = np.array([float(bounds[0][1]), float(bounds[1][1])])
    initial = np.clip(np.array([beta_init, az_init], dtype=float),
                      lower + 1e-10, upper - 1e-10)
    try:
        solved = least_squares(
            residual, initial, jac=jac, bounds=(lower, upper), method="trf",
            x_scale="jac", max_nfev=max_iter, ftol=1e-9, xtol=1e-9, gtol=1e-9)
    except _JacobianReferenceError as exc:
        reason = str(exc)
        rejected = stats.setdefault("jacobian_rejected_local_roots", {})
        rejected[reason] = rejected.get(reason, 0) + 1
        return None
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return None
    if float(np.linalg.norm(solved.fun)) >= tol_km:
        return None

    beta, az = map(float, solved.x)
    if jacobian == "central_reference_multiscale":
        try:
            reference_jacobian((beta, az))
        except _JacobianReferenceError as exc:
            reason = str(exc)
            rejected = stats.setdefault("jacobian_rejected_local_roots", {})
            rejected[reason] = rejected.get(reason, 0) + 1
            return None
    smoke = dict(trace((beta, az)))
    smoke.update({
        "beta_deg_tx": beta, "az_deg_tx": az,
        "elevation_tx_deg": beta,
        "p2p_search_profile": "smoke",
        "p2p_final_profile": formal_profile,
        "p2p_solver": "scipy_trf",
        "p2p_jacobian_method": jacobian,
        "p2p_jacobian_fallback_count": fallback_count,
        "p2p_jacobian_validation_relative": float(variational_error),
        "p2p_jacobian_column_relative": list(variational_column_errors),
        "p2p_jacobian_reference": reference_diagnostics,
        "p2p_jacobian_fallback_reasons": list(fallback_reasons),
        "p2p_jacobian_effective": (
            "central_reference_multiscale" if jacobian == "central_reference_multiscale" else
            "forward" if jacobian == "forward" else
            "variational" if variational_valid else "central_reference_fallback"),
        "p2p_landing_error_km": float(np.hypot(
            smoke["x_land_km"] - x_rx, smoke["y_land_km"] - y_rx)),
        "p2p_nfev": int(solved.nfev),
        "p2p_njev": int(solved.njev or 0),
        "p2p_seed_info": dict(seed_info or {}),
    })
    if defer_formal:
        return smoke
    return _formalize_ah_3d(
        tx, rx, n_ah, freq_MHz, smoke, tol_km, t, formal_profile,
        cache, stats)


def _formalize_ah_3d(tx, rx, n_ah, freq_MHz, local, tol_km, t,
                     formal_profile, cache, stats):
    beta = float(local["beta_deg_tx"])
    az = float(local["az_deg_tx"])
    ray_final = dict(_trace_cached(
        tx, n_ah, freq_MHz, beta, az, t, False, formal_profile,
        cache, stats))
    landing_error = float(np.hypot(
        ray_final["x_land_km"] - float(rx[0]), ray_final["y_land_km"]))
    if ray_final.get("convergence_status") != "passed" or landing_error >= tol_km:
        return None
    diagnostics = {key: value for key, value in local.items()
                   if key.startswith("p2p_")}
    ray_final.update(diagnostics)
    ray_final.update({
        "beta_deg_tx": beta, "az_deg_tx": az,
        "elevation_tx_deg": beta,
        "elevation_rx_deg": ray_final["beta_deg_rx"],
        "label": classify_mode_3d(ray_final),
        "p2p_landing_error_km": landing_error,
    })
    return ray_final


def _candidate_seeds(tx, rx, n_ah, freq_MHz, betas, rays_fan, t, nearest=3,
                     trace_cache=None, trace_stats=None, diagnostics=None):
    """Find full-domain crossing/fold candidates from a coarse landing map."""
    branch = n_ah.wave_mode
    fan = [r for r in (rays_fan or []) if r.get("wave_mode") == branch]
    if len(fan) < 3:
        cache = {} if trace_cache is None else trace_cache
        stats = trace_stats if trace_stats is not None else {
            "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
            "profile_wall_s": {}, "profile_nfev": {},
            "profile_accepted_steps": {}, "reference_ray_calls": 0,
            "reference_nfev": 0, "reference_accepted_steps": 0,
        }
        fan = [_trace_cached(tx, n_ah, freq_MHz, beta, 0.0, t, False,
                             "smoke", cache, stats) for beta in betas]
    fan = sorted((r for r in fan if r.get("branch_status") == "ground"),
                 key=lambda r: float(r.get("beta_deg_tx", 0.0)))
    if not fan:
        if diagnostics is not None:
            diagnostics.update(candidate_raw_count=0, candidate_retained_count=0,
                               candidate_removed_as_covered=0,
                               candidate_mapping=[])
        return [], 0
    target = float(rx[0])
    errors = np.array([float(r["x_land_km"]) - target for r in fan])
    raw = []

    def add(provenance, beta, az, support):
        raw.append({"beta_deg": float(beta), "az_deg": float(az),
                    "provenance": provenance,
                    "support_cells": sorted(set(int(value) for value in support))})

    for i, (left, right, e_left, e_right) in enumerate(
            zip(fan[:-1], fan[1:], errors[:-1], errors[1:])):
        if e_left == 0.0 or e_left * e_right < 0.0:
            add("crossing", 0.5 * (float(left["beta_deg_tx"]) +
                                    float(right["beta_deg_tx"])), 0.0, [i])
    for i in range(1, len(fan) - 1):
        x_left = float(fan[i]["x_land_km"] - fan[i - 1]["x_land_km"])
        x_right = float(fan[i + 1]["x_land_km"] - fan[i]["x_land_km"])
        if x_left * x_right <= 0.0 and (x_left != 0.0 or x_right != 0.0):
            add("fold", fan[i]["beta_deg_tx"], fan[i].get("az_deg_tx", 0.0),
                [i - 1, i])
    order = np.argsort(np.abs(errors))[:min(nearest, len(fan))]
    for i in order:
        add("nearest", fan[i]["beta_deg_tx"], fan[i].get("az_deg_tx", 0.0),
            ([i - 1] if i > 0 else []) + ([i] if i < len(fan) - 1 else []))

    unique = []
    for seed in raw:
        duplicate = next((old for old in unique
                          if seed["beta_deg"] == old["beta_deg"]
                          and seed["az_deg"] == old["az_deg"]), None)
        if duplicate is None:
            seed = dict(seed)
            seed["candidate_id"] = f"{branch}-{len(unique)}"
            seed["also_provenance"] = []
            unique.append(seed)
        else:
            duplicate["also_provenance"].append(seed["provenance"])

    physical = [seed for seed in unique if seed["provenance"] in ("crossing", "fold")]
    for seed in unique:
        if seed["provenance"] != "nearest":
            continue
        covered = [candidate["candidate_id"] for candidate in physical
                   if set(seed["support_cells"]) & set(candidate["support_cells"])]
        seed["covered_by_candidate_ids"] = covered
    if diagnostics is not None:
        diagnostics.update(
            candidate_raw_count=len(raw),
            candidate_retained_count=len(unique),
            candidate_removed_as_covered=0,
            candidate_mapping=[dict(seed) for seed in unique],
        )
    return unique, len(fan)


def _same_local_root(left, right, d_beta, d_az):
    if (left.get("wave_mode") != right.get("wave_mode")
            or left.get("branch_status") != right.get("branch_status")
            or left.get("turning_point_count") != right.get("turning_point_count")):
        return False
    return (
        abs(left["beta_deg_tx"] - right["beta_deg_tx"]) <= 1e-3 * d_beta
        and abs(left["az_deg_tx"] - right["az_deg_tx"]) <= 1e-3 * d_az
        and np.hypot(left["x_land_km"] - right["x_land_km"],
                     left["y_land_km"] - right["y_land_km"]) <= 1e-3
        and abs(left["h_reflect_km"] - right["h_reflect_km"]) <= 1e-3
        and abs(left["tau_ms"] - right["tau_ms"]) <= 1e-6
    )


def _merge_local_roots(local_candidates, d_beta, d_az):
    groups = []
    for candidate in local_candidates:
        group = next((item for item in groups
                      if _same_local_root(candidate, item["representative"],
                                          d_beta, d_az)), None)
        if group is None:
            groups.append({"representative": candidate, "members": [candidate]})
            continue
        group["members"].append(candidate)
        if candidate["p2p_landing_error_km"] < group["representative"]["p2p_landing_error_km"]:
            group["representative"] = candidate
    return groups


def _az_distance(left, right):
    return abs((float(left) - float(right) + 180.0) % 360.0 - 180.0)


def _trajectory_rms_km(left, right):
    def resample(ray):
        points = np.asarray(ray.get("trajectory_3d", []), dtype=float)
        if points.ndim != 2 or len(points) < 2 or points.shape[1] < 3:
            return None
        xyz = points[:, :3]
        distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]
        if distance[-1] <= 0.0:
            return None
        q = np.linspace(0.0, distance[-1], 33)
        return np.column_stack([np.interp(q, distance, xyz[:, i]) for i in range(3)])
    a, b = resample(left), resample(right)
    if a is None or b is None:
        return np.inf
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def _deduplicate(candidates: list, clust_h_km: float, clust_tau_ms: float,
                 clust_launch_deg: float = 0.02,
                 clust_trajectory_rms_km: float = 1.0) -> list:
    unique = []
    for c in candidates:
        is_dup = any(
            c.get("branch_id", c.get("wave_mode")) ==
            u.get("branch_id", u.get("wave_mode"))
            and c.get("branch_status") == u.get("branch_status")
            and c.get("turning_point_count") == u.get("turning_point_count")
            and abs(c["beta_deg_tx"] - u["beta_deg_tx"]) < clust_launch_deg
            and _az_distance(c["az_deg_tx"], u["az_deg_tx"]) < clust_launch_deg
            and _trajectory_rms_km(c, u) < clust_trajectory_rms_km
            and abs(c["h_reflect_km"] - u["h_reflect_km"]) < clust_h_km
            and abs(c["tau_ms"] - u["tau_ms"]) < clust_tau_ms
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
    diagnostics: dict | None = None,
) -> list:
    """Find AH3D P2P rays for O, X, or both branches."""
    n_init = int(p2p_params.get("n_init", 6))
    clust_h = float(p2p_params.get("clust_h_km", 10.0))
    clust_tau = float(p2p_params.get("clust_tau_ms", 0.05))
    clust_launch = float(p2p_params.get("clust_launch_deg", 0.02))
    clust_trajectory = float(p2p_params.get("clust_trajectory_rms_km", 1.0))
    tol_km = float(p2p_params.get("tol_km", 5.0))
    max_iter = int(p2p_params.get("max_iter", 12))
    d_beta = float(p2p_params.get("d_beta", 0.05))
    d_az = float(p2p_params.get("d_az", 0.10))
    jacobian = p2p_params.get("jacobian", "central_reference_multiscale")
    jacobian_rtol = float(p2p_params.get("jacobian_rtol", 1e-3))
    jacobian_scales = tuple(p2p_params.get("jacobian_scales", (0.5, 1.0, 2.0)))
    jacobian_min_scale = float(p2p_params.get("jacobian_min_scale", 0.125))
    jacobian_scale_rtol = float(p2p_params.get("jacobian_scale_rtol", 1e-2))
    jacobian_condition_max = float(
        p2p_params.get("jacobian_condition_max", 1e10))

    wm = str(wave_mode).upper()
    modes = ["O", "X"] if wm == "BOTH" else [wm]
    betas, seed_source = _seed_betas_3d(rx, p2p_params, None)

    candidates = []
    rejected_status = {}
    run_stats = {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {},
        "profile_accepted_steps": {}, "reference_ray_calls": 0,
        "reference_nfev": 0, "reference_accepted_steps": 0,
        "candidate_fan_wall_s": 0.0,
        "local_solve_wall_s": 0.0, "formal_research_wall_s": 0.0,
        "formal_reference_wall_s": 0.0, "dedup_wall_s": 0.0,
        "local_root_count": 0, "formal_root_count": 0,
        "merged_root_count": 0, "branches": {},
        "jacobian_rejected_local_roots": {},
    }
    for mode in modes:
        n_ah = RefractiveIndexAH3D(
            n3d, freq_MHz, wave_mode=mode, geomag=geomag,
            integration_profile="smoke",
        )
        trace_cache = {}
        candidate_diagnostics = {}
        stage_started = time.perf_counter()
        if str(p2p_params.get("audit_mode", "internal_scan")) == "adaptive_2d":
            from .p2p_mode_audit_3d import adaptive_landing_audit
            seeds, audit_diagnostics = adaptive_landing_audit(
                tx, rx, n_ah, freq_MHz, p2p_params, t,
                lambda beta, az: _trace_cached(
                    tx, n_ah, freq_MHz, beta, az, t, False, "smoke",
                    trace_cache, run_stats))
            audit_count = int(audit_diagnostics["node_count"])
            candidate_diagnostics.update(
                **audit_diagnostics,
                candidate_raw_count=len(seeds),
                candidate_retained_count=len(seeds),
                candidate_removed_as_covered=0,
                candidate_mapping=[dict(seed) for seed in seeds],
            )
        else:
            seeds, audit_count = _candidate_seeds(
                tx, rx, n_ah, freq_MHz, betas, rays_fan, t,
                nearest=int(p2p_params.get("candidate_nearest", 3)),
                trace_cache=trace_cache, trace_stats=run_stats,
                diagnostics=candidate_diagnostics)
        run_stats["candidate_fan_wall_s"] += time.perf_counter() - stage_started
        continued = [
            {"beta_deg": float(r.get("beta_deg_tx", r.get("elevation_tx_deg", 0.0))),
             "az_deg": float(r.get("az_deg_tx", 0.0)),
             "provenance": "continuation", "support_cells": [],
             "candidate_id": f"{mode}-continuation-{i}",
             "also_provenance": [], "covered_by_candidate_ids": []}
            for i, r in enumerate(continuation_modes or []) if r.get("wave_mode") == mode
        ]
        local_candidates = []
        stage_started = time.perf_counter()
        for seed in continued + seeds:
            beta, az = seed["beta_deg"], seed["az_deg"]
            local = _shoot_ah_3d(
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
                jacobian_scales=jacobian_scales,
                jacobian_min_scale=jacobian_min_scale,
                jacobian_scale_rtol=jacobian_scale_rtol,
                jacobian_condition_max=jacobian_condition_max,
                t=t,
                trace_cache=trace_cache,
                trace_stats=run_stats,
                defer_formal=True,
                seed_info=seed,
                angle_bounds=seed.get("angle_bounds"),
            )
            if local is None:
                continue
            local_candidates.append(local)
        run_stats["local_solve_wall_s"] += time.perf_counter() - stage_started
        groups = _merge_local_roots(local_candidates, d_beta, d_az)
        run_stats["local_root_count"] += len(local_candidates)
        run_stats["formal_root_count"] += len(groups)
        run_stats["merged_root_count"] += len(local_candidates) - len(groups)

        branch_groups = []
        for root_id, group in enumerate(groups):
            member_ids = [member["p2p_seed_info"].get("candidate_id")
                          for member in group["members"]]
            group["representative"]["p2p_merged_seed_ids"] = member_ids
            before_reference = run_stats["formal_reference_wall_s"]
            ray = _formalize_ah_3d(
                tx, rx, n_ah, freq_MHz, group["representative"], tol_km,
                t, "research", trace_cache, run_stats)
            if ray is not None:
                reference_wall = float(ray.get("convergence_reference_wall_s", 0.0))
                run_stats["formal_reference_wall_s"] = before_reference + reference_wall
            branch_groups.append({"root_id": f"{mode}-root-{root_id}",
                                  "member_candidate_ids": member_ids,
                                  "formal_passed": ray is not None})
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
                    ray["p2p_continuation_used"] = (
                        ray["p2p_seed_info"].get("provenance") == "continuation")
                    ray["p2p_local_root_count"] = len(local_candidates)
                    ray["p2p_formal_root_count"] = len(groups)
                    ray["p2p_merged_root_count"] = len(local_candidates) - len(groups)
                    candidates.append(ray)

        candidate_by_id = {item["candidate_id"]: item
                           for item in candidate_diagnostics.get("candidate_mapping", [])}
        root_by_candidate = {candidate_id: group["root_id"]
                             for group in branch_groups
                             for candidate_id in group["member_candidate_ids"]}
        for item in candidate_by_id.values():
            if item.get("provenance") != "nearest":
                continue
            own_root = root_by_candidate.get(item["candidate_id"])
            covered_roots = {root_by_candidate.get(candidate_id)
                             for candidate_id in item.get("covered_by_candidate_ids", [])}
            covered_roots.discard(None)
            item["shadow_root_id"] = own_root if own_root in covered_roots else None
            item["shadow_duplicate_proved"] = item["shadow_root_id"] is not None
        candidate_diagnostics["candidate_mapping"] = list(candidate_by_id.values())
        run_stats["branches"][mode] = {
            **candidate_diagnostics,
            "audit_ray_count": audit_count,
            "root_groups": branch_groups,
        }

    research_total = float(run_stats["profile_wall_s"].get("research", 0.0))
    run_stats["formal_research_wall_s"] = max(
        research_total - run_stats["formal_reference_wall_s"], 0.0)
    stage_started = time.perf_counter()
    unique = _deduplicate(candidates, clust_h, clust_tau,
                          clust_launch, clust_trajectory)
    unique.sort(key=lambda r: r["tau_ms"])
    annotate_mode_causes_3d(unique, p2p_params, seed_source)
    for ray in unique:
        ray["ah3d_rejected_status_counts"] = dict(rejected_status)
    run_stats["dedup_wall_s"] = time.perf_counter() - stage_started
    if diagnostics is not None:
        diagnostics.clear()
        diagnostics.update(run_stats)
    return unique
