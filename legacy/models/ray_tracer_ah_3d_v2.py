"""Legacy local-coordinate adaptive AH3D integration (A2).

Retained only until the C0R-5 caller migration.  The strict spherical-ECEF
Hamiltonian core is in ``ray_tracer_ah_3d_ecef``.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.integrate import DOP853

from config import RE_KM, RT, UPDATE2
from utils import freq_to_k0, group_delay_ms
from models.dispersion_ah import C_KM_S
from .ray_tracer_ah_3d import (
    RefractiveIndexAH3D, ah3d_deriv, group_index_ah3d,
    hamiltonian_ah3d, hamiltonian_residual_ah3d,
)


def _profile(name: str | None) -> tuple[str, dict]:
    profiles = UPDATE2["AH3D_INTEGRATION"]
    selected = str(name or UPDATE2.get("AH3D_INTEGRATION_PROFILE", "research"))
    if selected not in profiles:
        raise ValueError(f"unknown AH3D integration profile: {selected}")
    return selected, dict(profiles[selected])


def _project_dispersion(state, n_ah, k0, t, tol=1e-10, max_iter=5):
    """Minimum correction along D_k; never rescales |k| without diagnostics."""
    out = np.asarray(state, dtype=float).copy()
    total = 0.0
    for _ in range(max_iter):
        D = hamiltonian_ah3d(out[:6], n_ah, k0, t=t)
        if abs(D) / max(k0 * k0, 1e-30) <= tol:
            return out, total, True
        k = out[3:6]
        h = max(1e-7, 1e-6 * float(np.linalg.norm(k)))
        grad = np.empty(3)
        for i in range(3):
            yp, ym = out[:6].copy(), out[:6].copy()
            yp[3 + i] += h
            ym[3 + i] -= h
            grad[i] = (hamiltonian_ah3d(yp, n_ah, k0, t=t)
                       - hamiltonian_ah3d(ym, n_ah, k0, t=t)) / (2.0 * h)
        norm2 = float(np.dot(grad, grad))
        if not np.isfinite(norm2) or norm2 < 1e-24:
            break
        correction = -D * grad / norm2
        out[3:6] += correction
        total += float(np.linalg.norm(correction))
    return out, total, False


def _resolution_check(n_ah, spatial_fd_km: float) -> tuple[bool, str]:
    n3d = n_ah.n3d
    dx = float(np.min(np.diff(n3d.x_km)))
    dz = float(np.min(np.diff(n3d.z_km)))
    provider = getattr(n3d, "ne_provider", None)
    scales = []
    if provider is not None:
        for comp in getattr(provider, "components", []):
            meta = comp.metadata()
            for key in ("lambda_h_km", "lambda_z_km", "sigma_z_km",
                        "sigma_perp_km", "sigma_para_km", "delta_h_km",
                        "outer_scale_km"):
                value = meta.get(key)
                if value is not None and float(value) > 0:
                    scales.append(float(value))
            scales.extend(float(v) for v in meta.get("lambda_h_km_list", [])
                          if float(v) > 0)
            for key in ("inner_scale_m", "L_perp_m", "L_para_m"):
                value = meta.get(key)
                if value is not None and float(value) > 0:
                    scales.append(float(value) / 1000.0)
    if scales and min(scales) < 4.0 * max(dx, dz, spatial_fd_km):
        return True, "component_scale_below_4_samples"
    return False, "grid_and_component_scale_checked" if scales else "grid_only_no_component_scale"


def trace_single_ray_ah_3d_v2(tx_pos, beta_deg, az_deg,
                              n_ah: RefractiveIndexAH3D, freq_MHz,
                              rt_params=RT, t=0.0, profile=None,
                              with_sensitivity=False):
    """Trace one real AH ray with adaptive steps and D=0 projection."""
    profile_name, cfg = _profile(profile)
    k0 = freq_to_k0(freq_MHz)
    x0, z0 = map(float, tx_pos[:2])

    def launch_state(beta_value, az_value):
        beta, az = np.radians([beta_value, az_value])
        launch = np.array([np.cos(beta) * np.cos(az),
                           np.cos(beta) * np.sin(az), np.sin(beta)])
        root = n_ah.root(x0, 0.0, z0, launch, t=t)
        if root.status not in ("ok", "near_cutoff"):
            raise ValueError(f"invalid launch dispersion: {root.status}")
        return np.array([x0, 0.0, z0,
                         *(np.sqrt(root.n2.real) * k0 * launch), 0.0]), root

    base0, root0 = launch_state(beta_deg, az_deg)
    atol_base = np.array([cfg["atol_r_km"]] * 3 + [cfg["atol_k"]] * 3 + [1e-10])
    if with_sensitivity:
        h_angle = 1e-4
        sensitivity0 = np.column_stack([
            (launch_state(beta_deg + h_angle, az_deg)[0]
             - launch_state(beta_deg - h_angle, az_deg)[0]) / (2.0 * h_angle),
            (launch_state(beta_deg, az_deg + h_angle)[0]
             - launch_state(beta_deg, az_deg - h_angle)[0]) / (2.0 * h_angle),
        ])
        y0 = np.r_[base0, sensitivity0.ravel()]
        atol = np.r_[atol_base, np.full(sensitivity0.size, cfg["atol_r_km"])]
    else:
        y0 = base0
        atol = atol_base

    def base_rhs(state):
        base = ah3d_deriv(state[:6], n_ah, k0, t=t)
        ng = group_index_ah3d(*state[:3], state[3:6], n_ah, t=t)
        return np.r_[base, ng / C_KM_S]

    def rhs(s, state):
        base = base_rhs(state[:7])
        if not with_sensitivity:
            return base
        sensitivity = state[7:].reshape(7, 2)
        tangent = np.empty_like(sensitivity)
        for j in range(2):
            # Directional central differences integrate dS/ds=(df/dy)S
            # without forming the full 7x7 Hessian of the AH Hamiltonian.
            eps = min(1e-4, max(1e-7, 1e-2 /
                      max(float(np.linalg.norm(sensitivity[:3, j])), 1.0)))
            tangent[:, j] = (
                base_rhs(state[:7] + eps * sensitivity[:, j])
                - base_rhs(state[:7] - eps * sensitivity[:, j])
            ) / (2.0 * eps)
        return np.r_[base, tangent.ravel()]

    max_path = float(rt_params.get(
        "max_path_km",
        float(rt_params.get("max_steps", 20000))
        * float(rt_params.get("ds_max_km", RT["ds_max_km"]))))
    max_accepted_steps = int(np.ceil(max_path / cfg["ds_min_km"]))
    solver = DOP853(rhs, 0.0, y0, max_path, rtol=cfg["rtol"], atol=atol,
                    first_step=cfg["ds_min_km"], max_step=cfg["ds_max_km"])
    traj = [base0[:6].copy()]
    accepted = 0
    ds_used = []
    projection_max = 0.0
    residual_max = hamiltonian_residual_ah3d(base0[:6], n_ah, k0, t=t)
    event_error_max = 0.0
    statuses = [root0.status]
    previous_kz = base0[5]
    turning_count = 0
    has_risen = False
    termination = "max_path"
    h_max_state = base0.copy()

    while solver.status == "running":
        s_prev, y_prev = solver.t, solver.y.copy()
        solver.step()
        if solver.status == "failed" or not np.all(np.isfinite(solver.y)):
            termination = "step_failed"
            break
        ds_used.append(float(solver.t - s_prev))
        dense = solver.dense_output()
        base_trial = solver.y[:7].copy()
        projected, correction, ok = _project_dispersion(base_trial, n_ah, k0, t)
        projection_max = max(projection_max, correction)
        if not ok:
            termination = "constraint_failed"
            break
        if with_sensitivity:
            sensitivity = solver.y[7:].reshape(7, 2)
            projected_sensitivity = np.empty_like(sensitivity)
            for j in range(2):
                eps = min(1e-4, max(1e-7, 1e-2 /
                          max(float(np.linalg.norm(sensitivity[:3, j])), 1.0)))
                plus, _, plus_ok = _project_dispersion(
                    base_trial + eps * sensitivity[:, j], n_ah, k0, t)
                minus, _, minus_ok = _project_dispersion(
                    base_trial - eps * sensitivity[:, j], n_ah, k0, t)
                if not (plus_ok and minus_ok):
                    ok = False
                    break
                projected_sensitivity[:, j] = (plus - minus) / (2.0 * eps)
            if not ok:
                termination = "constraint_failed"
                break
            solver.y[:] = np.r_[projected, projected_sensitivity.ravel()]
        else:
            solver.y[:] = projected
        solver.f = rhs(solver.t, solver.y)
        accepted += 1
        state = solver.y[:7]
        residual_max = max(residual_max,
                           hamiltonian_residual_ah3d(state[:6], n_ah, k0, t=t))
        root = n_ah.root(*state[:3], state[3:6], t=t)
        if root.status not in statuses:
            statuses.append(root.status)
        if root.status in ("invalid_dispersion", "resonance_excluded", "branch_lost"):
            termination = root.status
            break
        if previous_kz * state[5] < 0.0 and state[2] > 50.0:
            lo, hi = s_prev, solver.t
            flo = float(dense(lo)[5])
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                fm = float(dense(mid)[5])
                if flo * fm <= 0.0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            event_error_max = max(event_error_max, hi - lo)
            turning_count += 1
            if "turning_point" not in statuses:
                statuses.append("turning_point")
        previous_kz = state[5]
        traj.append(state[:6].copy())
        if state[2] > h_max_state[2]:
            h_max_state = state.copy()
        has_risen = has_risen or state[2] > 50.0
        if state[2] <= 0.0 and has_risen:
            frac = y_prev[2] / max(y_prev[2] - state[2], 1e-30)
            event_base = y_prev[:7] + np.clip(frac, 0.0, 1.0) * (state - y_prev[:7])
            if with_sensitivity:
                s_prev_mat = y_prev[7:].reshape(7, 2)
                s_now = solver.y[7:].reshape(7, 2)
                q = np.clip(frac, 0.0, 1.0)
                fixed_s = s_prev_mat + q * (s_now - s_prev_mat)
                event_rhs = base_rhs(event_base)
                if abs(event_rhs[2]) < 1e-12:
                    termination = "ground_event_tangent"
                    break
                # Continuous event-surface correction for g(y)=z=0:
                # S_event = S - f * (g_y S)/(g_y f).
                event_s = fixed_s - np.outer(
                    event_rhs, fixed_s[2] / event_rhs[2])
                solver.y[:] = np.r_[event_base, event_s.ravel()]
            else:
                solver.y[:] = event_base
            traj[-1] = event_base[:6].copy()
            termination = "ground"
            break
        if state[2] >= float(rt_params["z_stop_km"]):
            termination = "upper_boundary"
            break
        if accepted >= max_accepted_steps:
            termination = "max_steps"
            break

    final = solver.y[:7]
    kmag = max(float(np.linalg.norm(final[3:6])), 1e-30)
    h_ref = float(h_max_state[2])
    if n_ah.earth_flat:
        h_ref -= (h_max_state[0] ** 2 + h_max_state[1] ** 2) / (2.0 * RE_KM)
    under_resolved, resolution_reason = _resolution_check(n_ah, 0.5)
    converged = (termination == "ground" and turning_count > 0 and
                 residual_max <= 1e-6 and not under_resolved)
    # DOP853 uses 12 RHS stages per accepted step; excess stages diagnose retries.
    rejected = max(0, int(round((solver.nfev - 12 * accepted - 2) / 12)))
    geom_path = sum(float(np.linalg.norm(np.diff(np.asarray(traj)[:, :3], axis=0)[i]))
                    for i in range(max(len(traj) - 1, 0)))
    result = {
        "x_land_km": float(final[0]), "y_land_km": float(final[1]),
        "z_land_km": float(final[2]), "geom_path_km": geom_path,
        "group_path_km": float(C_KM_S * final[6]),
        "tau_ms": group_delay_ms(C_KM_S * final[6]),
        "h_reflect_km": h_ref, "reflect_x_km": float(h_max_state[0]),
        "reflect_y_km": float(h_max_state[1]),
        "beta_deg_tx": float(beta_deg), "az_deg_tx": float(az_deg),
        "elevation_tx_deg": float(beta_deg),
        "beta_deg_rx": float(np.degrees(np.arcsin(np.clip(abs(final[5]) / kmag, 0, 1)))),
        "elevation_rx_deg": float(np.degrees(np.arcsin(np.clip(abs(final[5]) / kmag, 0, 1)))),
        "az_deg_rx": float(np.degrees(np.arctan2(final[4], final[3]))),
        "azimuth_deflect_deg": float(np.degrees(np.arctan2(final[1], max(abs(final[0]), 1e-6)))),
        "trajectory_3d": traj, "wave_mode": n_ah.wave_mode,
        "branch_id": n_ah.wave_mode, "branch_status": termination,
        "branch_status_history": ",".join(statuses),
        "integration_profile": profile_name, "integrator": "DOP853",
        "accepted_steps": accepted, "rejected_steps_estimate": rejected,
        "nfev": int(solver.nfev),
        "ds_min_used": min(ds_used, default=0.0),
        "ds_max_used": max(ds_used, default=0.0),
        "hamiltonian_residual_max": float(residual_max),
        "dispersion_projection_correction_max": float(projection_max),
        "event_localization_error": float(event_error_max),
        "turning_point_count": turning_count,
        "under_resolved_Ne": under_resolved,
        "resolution_check": resolution_reason,
        "resolution_scope": ("grid_and_component_scale" if "component" in resolution_reason
                             else "grid_only"),
        "convergence_status": ("screening_only" if profile_name == "smoke"
                               else "passed" if converged else "convergence_failed"),
        "group_index_model": "dispersion_gradient_-Dk_Domega",
        "group_path_definition": "c_integral_dt",
        "ray_attenuation_semantics": "real_ray_complex_attenuation",
        "n2_floor_propagation": False, "ox_solver": "ah3d_v2_adaptive",
        "solver_dimensionality": "3d_ray", "freq_MHz": float(freq_MHz),
        "theory_reference": "Appleton-Hartree; Haselgrove; Jones-Stephenson; DOP853",
    }
    if with_sensitivity:
        landing_jacobian = solver.y[7:].reshape(7, 2)[:2]
        result.update({
            "landing_jacobian_deg": landing_jacobian.tolist(),
            "landing_jacobian_method": "integrated_tangent_linear_central_rhs",
            "landing_jacobian_condition": float(np.linalg.cond(landing_jacobian)),
            "landing_event_sensitivity_corrected": termination == "ground",
        })
    if profile_name == "research":
        reference_started = time.perf_counter()
        reference = trace_single_ray_ah_3d_v2(
            tx_pos, beta_deg, az_deg, n_ah, freq_MHz, rt_params, t=t,
            profile="reference", with_sensitivity=False,
        )
        reference_wall_s = time.perf_counter() - reference_started
        landing_delta = float(np.hypot(result["x_land_km"] - reference["x_land_km"],
                                       result["y_land_km"] - reference["y_land_km"]))
        height_delta = abs(result["h_reflect_km"] - reference["h_reflect_km"])
        delay_rel = abs(result["tau_ms"] - reference["tau_ms"]) / max(
            abs(reference["tau_ms"]), 1e-12)
        result.update({
            "convergence_landing_delta_km": landing_delta,
            "convergence_reflection_height_delta_km": height_delta,
            "convergence_group_delay_relative": delay_rel,
            "convergence_reference_profile": "reference",
            "convergence_reference_wall_s": reference_wall_s,
            "convergence_reference_nfev": int(reference.get("nfev", 0)),
            "convergence_reference_accepted_steps": int(
                reference.get("accepted_steps", 0)),
        })
        passed = (result["convergence_status"] == "passed"
                  and reference["convergence_status"] == "passed"
                  and landing_delta <= 1.0 and height_delta <= 0.5
                  and delay_rel <= 1e-3)
        result["convergence_status"] = "passed" if passed else "convergence_failed"
    return result
