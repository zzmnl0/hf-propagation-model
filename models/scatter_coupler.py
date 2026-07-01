"""
Unified scatter-regime coupler for 3-D ray modes.

This module connects the existing Born, Rytov, and MPS selector formulas to
3-D mode dictionaries. It is diagnostic by default: it annotates each mode
with the selected scatter regime and theory outputs without changing Pr_W.
"""
import numpy as np

from .scatter_selector import scatter_phase_variance, select_scatter_method
from .scatter_born import born_scatter
from .scatter_rytov import rytov_full


def _estimate_ne0_from_ray(mode: dict, n3d, fallback: float) -> float:
    """Estimate representative Ne0 from samples along a 3-D trajectory."""
    traj = mode.get("trajectory_3d")
    if not traj:
        return float(fallback)

    samples = traj[::max(len(traj) // 64, 1)]
    vals = []
    for s in samples:
        try:
            vals.append(float(n3d.Ne_at(float(s[0]), float(s[1]), float(s[2]))))
        except Exception:
            continue
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
    """
    Annotate 3-D modes with Born/Rytov/MPS scatter-regime diagnostics.

    Fields added include:
      scatter_method, sigma_phi_sq, sigma_phi_rad, scatter_Ne0_m3,
      scatter_Cs_rel, scatter_L_outer_km, scatter_DeltaL_km,
      plus Born or Rytov outputs where applicable.
    """
    Cs_rel = float(scatter_params.get("Cs_rel", 0.01))
    p_spec = float(scatter_params.get("p", 3.0))
    L_outer_km = float(scatter_params.get("L_outer_km", 10.0))
    DeltaL_km = float(scatter_params.get("DeltaL_km", 100.0))
    z_eff_default = float(scatter_params.get("z_eff_km", 300.0))
    Ne0_default = float(scatter_params.get("Ne0_m3", 5.0e10))

    for mode in modes:
        Ne0 = _estimate_ne0_from_ray(mode, n3d, Ne0_default)
        sigma_phi_sq = scatter_phase_variance(
            Ne0, Cs_rel, DeltaL_km, L_outer_km, freq_MHz)
        method = select_scatter_method(sigma_phi_sq)

        mode.update({
            "scatter_method": method,
            "sigma_phi_sq": float(sigma_phi_sq),
            "sigma_phi_rad": float(np.sqrt(max(sigma_phi_sq, 0.0))),
            "scatter_Ne0_m3": float(Ne0),
            "scatter_Cs_rel": Cs_rel,
            "scatter_L_outer_km": L_outer_km,
            "scatter_DeltaL_km": DeltaL_km,
            "mps_required": False,
        })

        if method == "Born":
            mode.update(born_scatter(
                mode, Pt_W, Gt, Gr, freq_MHz,
                Cs_rel, p_spec, L_outer_km, Ne0))
        elif method == "Rytov":
            z_eff_km = float(mode.get("h_reflect_km", z_eff_default))
            mode.update(rytov_full(
                freq_MHz, Cs_rel, p_spec, L_outer_km,
                Ne0, z_eff_km, DeltaL_km))
            mode["scatter_method"] = "Rytov"
        else:
            mode["mps_required"] = True
            mode["scatter_method"] = "MPS" if enable_mps else "MPS_required"

    return modes
