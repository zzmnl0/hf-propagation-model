"""
Stage D verification: 3-D P2P forward vs central Jacobian option.

Run:
    conda run -n pytorch_cpu python tests/_verify_p2p_jacobian.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.ionosphere_model import IonosphereModel
from models.ray_tracer_3d import RefractiveIndex3D
from models.point_to_point_3d import find_all_rays_3d


def run():
    print("=" * 55)
    print("  Stage D verification: P2P Jacobian")
    print("=" * 55)

    iono = IonosphereModel()
    Ne, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d = RefractiveIndex3D(Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ, earth_flat=True)

    base = {**cfg.P2P_3D, "n_init": 6, "max_iter": 12}
    modes_f = find_all_rays_3d(
        cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ,
        p2p_params={**base, "jacobian": "forward"})
    modes_c = find_all_rays_3d(
        cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ,
        p2p_params={**base, "jacobian": "central"})

    n_fail = 0
    ok_count = bool(modes_f) and bool(modes_c)
    print("  forward_modes={} central_modes={} {}".format(
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

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
