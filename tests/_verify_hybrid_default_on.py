"""Verify reflection/scatter hybrid coupling is enabled by default."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import config as cfg
from models.hybrid_model import build_pd_spectrum


def run():
    print("=" * 55)
    print(" hybrid default-on verification")
    print("=" * 55)
    n_fail = 0
    ok_pe = cfg.UPDATE2["PE_COUPLING"].get("enable") is True and cfg.UPDATE2["PE_COUPLING"].get("apply_to_power") is True
    ok_pd = cfg.UPDATE2["SCATTER_PD"].get("apply_to_pd") is True
    print(" PE_COUPLING default_on {}".format("OK" if ok_pe else "FAIL")); n_fail += 0 if ok_pe else 1
    print(" SCATTER_PD default_on {}".format("OK" if ok_pd else "FAIL")); n_fail += 0 if ok_pd else 1
    modes = [{"label":"m", "tau_ms":3.0, "delta_tau_ms":0.02, "Pr_W":1.0, "scatter_method":"Rytov", "sigma_tau_ms":0.2}]
    tau, pd = build_pd_spectrum(modes)
    ok_power = abs(float(np.sum(pd)) - 1.0) < 1e-9
    ok_spread = float(np.max(pd)) < 1.0
    print(" default_scatter_pd_power_conserved {}".format("OK" if ok_power else "FAIL")); n_fail += 0 if ok_power else 1
    print(" default_scatter_pd_broadened {}".format("OK" if ok_spread else "FAIL")); n_fail += 0 if ok_spread else 1
    print("{} fail(s)".format(n_fail))
    return n_fail

if __name__ == "__main__":
    sys.exit(run())
