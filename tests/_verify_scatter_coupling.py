"""
Stage C verification: unified 3-D scatter coupler diagnostics.

Run:
    conda run -n pytorch_cpu python tests/_verify_scatter_coupling.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.scatter_coupler import apply_scatter_regime_3d


def _mode():
    return {
        "label": "test",
        "group_path_km": 1000.0,
        "tau_ms": 1000.0 / cfg.C_KMS * 1e3,
        "h_reflect_km": 300.0,
    }


def _run_case(name, cs_rel, expected):
    params = {
        **cfg.SCATTER,
        "Cs_rel": cs_rel,
        "Ne0_m3": 5.0e10,
        "DeltaL_km": 100.0,
        "L_outer_km": 10.0,
    }
    modes = [_mode()]
    apply_scatter_regime_3d(
        modes, None, params, cfg.FREQ_MHZ, cfg.PT_W, cfg.GT, cfg.GR)
    m = modes[0]
    ok = m["scatter_method"] == expected
    common = all(k in m for k in (
        "sigma_phi_sq", "sigma_phi_rad", "scatter_Ne0_m3",
        "scatter_Cs_rel", "scatter_L_outer_km", "scatter_DeltaL_km"))
    if expected == "Born":
        common = common and "Pr_scatter_W" in m and "sigma_v_m1" in m
    if expected == "Rytov":
        common = common and "S4" in m and "sigma_tau_ms" in m
    if expected == "MPS_required":
        common = common and m.get("mps_required") is True
    ok = ok and common
    print("  {} Cs={} method={} sigma_phi_sq={:.3e} {}".format(
        name, cs_rel, m["scatter_method"], m["sigma_phi_sq"],
        "OK" if ok else "FAIL"))
    return 0 if ok else 1


def run():
    print("=" * 55)
    print("  Stage C verification: scatter coupler")
    print("=" * 55)
    n_fail = 0
    n_fail += _run_case("Born", 1.0e-3, "Born")
    n_fail += _run_case("Rytov", 1.0e-2, "Rytov")
    n_fail += _run_case("MPS", 1.0e-1, "MPS_required")
    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
