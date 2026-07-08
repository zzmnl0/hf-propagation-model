import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import config as cfg
from models.hybrid_model import build_pd_spectrum
from models.pd_scatter_kernel import apply_scatter_kernel


def _modes(method="Rytov"):
    return [{"label": "m1", "tau_ms": 3.0, "delta_tau_ms": 0.02, "Pr_W": 2.0, "scatter_method": method, "sigma_tau_ms": 0.20}]


def run():
    print("=" * 55)
    print(" P2 verification: scatter P-D kernel")
    print("=" * 55)
    old = dict(cfg.UPDATE2["SCATTER_PD"])
    n_fail = 0
    try:
        cfg.UPDATE2["SCATTER_PD"].update({"apply_to_pd": False})
        tau0, pd0 = build_pd_spectrum(_modes())
        cfg.UPDATE2["SCATTER_PD"].update({"apply_to_pd": True})
        tau1, pd1 = build_pd_spectrum(_modes())
        ok_changed = abs(float(np.sum(pd1)) - 2.0) < 1e-9 and float(np.max(pd1)) < 2.0
        print(" apply_on_power_conserved {}".format("OK" if ok_changed else "FAIL")); n_fail += 0 if ok_changed else 1
        tau_mps, pd_mps = apply_scatter_kernel(_modes("MPS_required"), tau0)
        tau_base, pd_base = build_pd_spectrum(_modes("MPS_required"), tau_axis=tau0, method="gaussian")
        ok_mps = np.allclose(tau_mps, tau_base) and abs(float(np.sum(pd_mps)) - 2.0) < 1e-9
        print(" mps_required_no_fake_power {}".format("OK" if ok_mps else "FAIL")); n_fail += 0 if ok_mps else 1
        cfg.UPDATE2["SCATTER_PD"].update({"apply_to_pd": False})
        tau_off, pd_off = build_pd_spectrum(_modes())
        ok_off = np.allclose(tau_off, tau0) and np.allclose(pd_off, pd0)
        print(" default_off_legacy_pd {}".format("OK" if ok_off else "FAIL")); n_fail += 0 if ok_off else 1
    finally:
        cfg.UPDATE2["SCATTER_PD"].clear(); cfg.UPDATE2["SCATTER_PD"].update(old)
    print("{} fail(s)".format(n_fail))
    return n_fail

if __name__ == "__main__":
    sys.exit(run())
