"""Verify reflection/scatter hybrid coupling is disabled by default."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import config as cfg
from models.hybrid_model import build_pd_spectrum


def run():
    print("=" * 55)
    print(" hybrid default-off verification")
    print("=" * 55)
    n_fail = 0
    ok_pe = cfg.UPDATE2["PE_COUPLING"].get("enable") is False and cfg.UPDATE2["PE_COUPLING"].get("apply_to_power") is False
    ok_pd = cfg.UPDATE2["SCATTER_PD"].get("apply_to_pd") is False
    print(" PE_COUPLING default_off {}".format("OK" if ok_pe else "FAIL")); n_fail += 0 if ok_pe else 1
    print(" SCATTER_PD default_off {}".format("OK" if ok_pd else "FAIL")); n_fail += 0 if ok_pd else 1
    modes = [{"label":"m", "tau_ms":3.0, "delta_tau_ms":0.02, "Pr_W":1.0, "scatter_method":"Rytov", "sigma_tau_ms":0.2}]
    tau, pd = build_pd_spectrum(modes)
    ok_power = abs(float(np.max(pd)) - 1.0) < 1e-9
    ok_spread = float(np.max(pd)) == 1.0
    print(" default_legacy_pd_peak_preserved {}".format("OK" if ok_power else "FAIL")); n_fail += 0 if ok_power else 1
    print(" default_scatter_pd_unchanged {}".format("OK" if ok_spread else "FAIL")); n_fail += 0 if ok_spread else 1
    print("{} fail(s)".format(n_fail))
    return n_fail

if __name__ == "__main__":
    sys.exit(run())


