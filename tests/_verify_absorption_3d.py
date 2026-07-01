"""
Stage B verification: 3-D D-layer absorption in main._add_power_3d().

Run:
    conda run -n pytorch_cpu python tests/_verify_absorption_3d.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import copy
import numpy as np

import config as cfg
from main import _add_power_3d
from utils import d_layer_absorption_dB


def _mode():
    return {
        "group_path_km": 1000.0,
        "beta_deg_tx": 30.0,
        "tau_ms": 1000.0 / cfg.C_KMS * 1e3,
        "h_reflect_km": 250.0,
    }


def _ok(label, condition):
    print("  {} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run():
    n_fail = 0
    old_abs = copy.deepcopy(cfg.ABSORPTION)
    try:
        print("=" * 55)
        print("  Stage B verification: 3-D absorption")
        print("=" * 55)

        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update({"enable": False, "A0": 500.0, "chi_deg": 60.0})
        clean = [_mode()]
        _add_power_3d(clean)
        clean_p = clean[0]["Pr_W"]
        n_fail += _ok("disabled abs_total_dB == 0:",
                      abs(clean[0]["abs_total_dB"]) < 1e-12)

        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update({"enable": True, "A0": 500.0, "chi_deg": 60.0})
        comm = [_mode()]
        _add_power_3d(comm)
        fH_L = cfg.GEOMAG["fH_MHz"] * np.cos(np.radians(cfg.GEOMAG["dip_deg"]))
        A = d_layer_absorption_dB(
            cfg.FREQ_MHZ, 30.0, A0=500.0, chi_deg=60.0, fH_L_MHz=fH_L)
        expected_comm = clean_p * 10.0**(-A / 10.0)
        print("\n[1] Communication one-way absorption")
        print("  A={:.6f} dB clean={:.6e} abs={:.6e}".format(
            A, clean_p, comm[0]["Pr_W"]))
        n_fail += _ok("Pr ratio matches one-way A:",
                      abs(comm[0]["Pr_W"] / expected_comm - 1.0) < 1e-10)
        n_fail += _ok("abs_total_dB == A:",
                      abs(comm[0]["abs_total_dB"] - A) < 1e-10)

        radar_clean = [_mode()]
        cfg.ABSORPTION["enable"] = False
        _add_power_3d(radar_clean, radar_mode=True, sigma_rcs=1.0)
        radar_abs = [_mode()]
        cfg.ABSORPTION["enable"] = True
        _add_power_3d(radar_abs, radar_mode=True, sigma_rcs=1.0)
        expected_radar = radar_clean[0]["Pr_W"] * 10.0**(-(2.0*A) / 10.0)
        print("\n[2] Radar two-way absorption")
        print("  2A={:.6f} dB clean={:.6e} abs={:.6e}".format(
            2.0*A, radar_clean[0]["Pr_W"], radar_abs[0]["Pr_W"]))
        n_fail += _ok("Pr ratio matches two-way A:",
                      abs(radar_abs[0]["Pr_W"] / expected_radar - 1.0) < 1e-10)
        n_fail += _ok("abs_total_dB == 2A:",
                      abs(radar_abs[0]["abs_total_dB"] - 2.0*A) < 1e-10)

        print("\n{} fail(s)\n".format(n_fail))
        return n_fail
    finally:
        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update(old_abs)


if __name__ == "__main__":
    sys.exit(run())
