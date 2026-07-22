"""Unified scatter-regime coupler for 3-D ray modes.

P3 keeps scatter coupling explicit and diagnostic by default. Born/Rytov
statistics annotate weak or moderate regimes; MPS is used only when explicitly
enabled and then remains a path-related phase-screen statistical kernel, not a
complete 3-D full-wave solution.
"""

from __future__ import annotations

import numpy as np

from .scatter_born import born_scatter
from .scatter_phase_screen import make_path_mps_delay_kernel
from .scatter_rytov import rytov_full
from .scatter_selector import scatter_phase_variance, select_scatter_method


def _estimate_ne0_from_ray(mode: dict, n3d, fallback: float) -> float:
    """Estimate Ne0 from samples along a 3-D trajectory."""
    traj = mode.get("trajectory_ecef_km")
    if traj is None or len(traj) == 0:
        raise ValueError("formal scatter diagnostics require trajectory_ecef_km")
    if not hasattr(n3d, "get_Ne_ecef"):
        raise TypeError("formal scatter sampling requires an ECEF NeProvider")
    samples = traj[::max(len(traj) // 64, 1)]
    vals = []
    for s in samples:
        try:
            vals.append(float(n3d.get_Ne_ecef(
                np.asarray(s)[:3], float(mode["earth_radius_km"]))))
        except Exception:
            pass
    if not vals:
        return float(fallback)
    return float(max(max(vals), 1.0))


def apply_scatter_regime_3d(modes: list,
                            n3d,
                            scatter_params: dict,
                            freq_MHz: float,
                            Pt_W: float,
                            Gt: float,
                            Gr: float,
                            enable_mps: bool = False) -> list:
    """Annotate modes with physically separated scatter diagnostics."""
    Cs_rel = float(scatter_params.get("Cs_rel", 0.01))
    p = float(scatter_params.get("p", 3.0))
    L_outer_km = float(scatter_params.get("L_outer_km", 10.0))
    DeltaL_km = float(scatter_params.get("DeltaL_km", 100.0))
    z_eff_km = float(scatter_params.get("z_eff_km", 300.0))
    Ne0_fallback = float(scatter_params.get("Ne0_m3", 1.0e11))

    for mode in modes:
        if str(mode.get("branch_status", "ok")) not in ("ok", "near_cutoff", ""):
            mode.update({
                "scatter_method": "not_applied",
                "scatter_validity": "invalid_or_lost_ray",
                "scatter_power_applied": False,
                "kernel_normalization": 0.0,
            })
            continue

        Ne0 = _estimate_ne0_from_ray(mode, n3d, Ne0_fallback)
        sigma_phi_sq = scatter_phase_variance(
            Ne0, Cs_rel, DeltaL_km, L_outer_km, freq_MHz)
        method = select_scatter_method(sigma_phi_sq)
        sigma_phi = float(np.sqrt(max(sigma_phi_sq, 0.0)))
        mode.update({
            "scatter_method": method,
            "sigma_phi_sq": float(sigma_phi_sq),
            "sigma_phi_rad": sigma_phi,
            "scatter_Ne0_m3": Ne0,
            "scatter_Cs_rel": Cs_rel,
            "scatter_L_outer_km": L_outer_km,
            "scatter_DeltaL_km": DeltaL_km,
            "scatter_power_applied": False,
            "scatter_solver_dimensionality": "statistical_diagnostic",
            "scatter_validity": "born_rytov_mps_regime_diagnostic",
            "scatter_path_coordinate_model": "spherical_ECEF",
            "theory_reference": "Rino 1979; Yeh and Liu 1982; Knepp 1983",
        })

        if method == "Born":
            mode.update(born_scatter(
                mode, Pt_W, Gt, Gr, freq_MHz,
                Cs_rel, p, L_outer_km, Ne0))
            mode["scatter_validity"] = "weak_scatter_born_statistical"
        elif method == "Rytov":
            mode.update(rytov_full(
                freq_MHz, Cs_rel, p,
                L_outer_km, Ne0, z_eff_km, DeltaL_km))
            mode["scatter_method"] = "Rytov"
            mode["scatter_validity"] = "moderate_scatter_rytov_statistical"
        else:
            mode["mps_required"] = True
            if enable_mps:
                diag = make_path_mps_delay_kernel(
                    mode, n3d, np.linspace(
                        max(0.0, float(mode.get("tau_ms", 0.0)) - 1.0),
                        float(mode.get("tau_ms", 0.0)) + 1.0,
                        101),
                    scatter_params=scatter_params,
                    phase_params=scatter_params.get("phase_screen", None),
                    freq_MHz=freq_MHz)
                mode.update(diag)
            else:
                mode["scatter_method"] = "MPS_required"
                mode["scatter_validity"] = "mps_required_not_applied"
                mode["kernel_normalization"] = 0.0

    return modes
