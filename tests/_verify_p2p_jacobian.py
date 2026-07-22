"""
AH3D verification: variational vs central-reference P2P Jacobian.

Run:
    conda run -n pytorch_cpu python tests/_verify_p2p_jacobian.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.point_to_point_ah_3d import find_all_rays_ah_3d
from models.scenario_provider_factory import ScenarioProviderFactory


def run():
    print("=" * 55)
    print("  AH3D verification: P2P Jacobian")
    print("=" * 55)

    provider = ScenarioProviderFactory(cfg).create(
        "baseline", seed=0, t_s=0.0, validation=True)
    ne_slice = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z, t=0.0, y_km=0.0)
    n3d = RefractiveIndex3D(
        ne_slice, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ, earth_flat=True,
        ne_provider=provider)

    base = {**cfg.P2P_3D, "n_init": 6, "max_iter": 12}
    modes_f = find_all_rays_ah_3d(
        cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ,
        p2p_params={**base, "jacobian": "variational"}, wave_mode="O")
    diagnostics = {}
    modes_c = find_all_rays_ah_3d(
        cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ,
        p2p_params={**base, "jacobian": "central_reference_multiscale"},
        wave_mode="O", diagnostics=diagnostics)

    n_fail = 0
    ok_count = bool(modes_f) and bool(modes_c)
    print("  variational_modes={} central_modes={} {}".format(
        len(modes_f), len(modes_c), "OK" if ok_count else "FAIL"))
    if not ok_count:
        return 1

    dtau = abs(modes_f[0]["tau_ms"] - modes_c[0]["tau_ms"])
    dh = abs(modes_f[0]["h_reflect_km"] - modes_c[0]["h_reflect_km"])
    ok_close = dtau < 0.2 and dh < 20.0
    print("  dtau={:.6f} ms dh={:.3f} km {}".format(
        dtau, dh, "OK" if ok_close else "FAIL"))
    if not ok_close:
        n_fail += 1

    checked = [m for m in modes_c if m.get("p2p_jacobian_reference")]
    if not checked:
        print("  jacobian_diagnostics=missing FAIL")
        n_fail += 1
    else:
        mode = checked[0]
        meta = mode["p2p_jacobian_reference"]
        ok_diag = (meta.get("stable") is True
                   and len(meta.get("relative_changes", [])) >= 2
                   and np.isfinite(meta.get("condition", np.nan))
                   and meta["condition"] <= cfg.P2P_3D["jacobian_condition_max"]
                   and mode.get("p2p_jacobian_effective")
                   == "central_reference_multiscale")
        print("  reference_scale={} condition={:.3e} {}".format(
            meta.get("selected_scale"), float(meta.get("condition", np.nan)),
            "OK" if ok_diag else "FAIL"))
        if not ok_diag:
            n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
