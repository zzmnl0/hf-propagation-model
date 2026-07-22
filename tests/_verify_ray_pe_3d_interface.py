import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import config as cfg
from legacy.models.hybrid_ray_pe import apply_ray_pe_coupling_3d
from legacy.models.pe_propagator import PEPropagator


def _mode():
    states = []
    for x in np.linspace(540.0, 660.0, 25):
        states.append(np.array([x, 0.0, 350.0, 1.0, 0.0, 0.0], dtype=float))
    return {"trajectory_3d": states, "Pr_W": 1e-9}


def run():
    print("=" * 55)
    print(" P2 verification: Ray-PE 3D interface")
    print("=" * 55)
    x = np.linspace(500.0, 700.0, 101)
    z = np.linspace(300.0, 400.0, 101)
    Ne = np.zeros((len(x), len(z)), dtype=float)
    pe = PEPropagator(freq_MHz=cfg.FREQ_MHZ, pe_params={**cfg.PE, "dx_km": 2.0, "dz_m": 1000.0, "n_pml": 0, "store_history": False})
    bubble = {**cfg.BUBBLE_3D, "enable": True, "x0_km": 600.0, "y0_km": 0.0, "z0_km": 350.0, "sigma_perp_km": 80.0, "sigma_para_km": 300.0}
    Pr_in = 1e-9
    Pr_diag, diag = apply_ray_pe_coupling_3d(pe, _mode(), Ne, x, z, Pr_in, pe.k0, bubble, {"enable": True, "apply_to_power": False})
    Pr_apply, diag_apply = apply_ray_pe_coupling_3d(pe, _mode(), Ne, x, z, Pr_in, pe.k0, bubble, {"enable": True, "apply_to_power": True})
    n_fail = 0
    ok_diag = diag.get("pe_coupled") and diag.get("pe_frame_status") == "local_2p5d" and "pe_y_entry_km" in diag
    ok_power = abs(Pr_diag - Pr_in) <= 1e-21
    ok_apply = abs(Pr_apply - Pr_in * diag_apply.get("pe_power_ratio", 1.0)) <= 1e-21
    print(" diagnostic_fields {}".format("OK" if ok_diag else "FAIL")); n_fail += 0 if ok_diag else 1
    print(" diagnostic_power_unchanged {}".format("OK" if ok_power else "FAIL")); n_fail += 0 if ok_power else 1
    print(" apply_to_power_ratio {}".format("OK" if ok_apply else "FAIL")); n_fail += 0 if ok_apply else 1
    print("{} fail(s)".format(n_fail))
    return n_fail

if __name__ == "__main__":
    sys.exit(run())
