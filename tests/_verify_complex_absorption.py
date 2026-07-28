"""Verification for update2.1 AH-collisional complex absorption.

Run:
    conda run -n pytorch_cpu python tests/_verify_complex_absorption.py
"""

from __future__ import annotations

import copy
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from main import _add_power_3d
from models.absorption_collision import (
    complex_absorption_path_dB,
    complex_refractive_index,
    complex_refractive_index_ah,
)
from utils import geodetic_to_ecef


class _ConstN3D:
    def __init__(self, ne):
        self.ne = float(ne)

    def get_Ne_ecef(self, ecef_km, earth_radius_km=None, t=0.0):
        return self.ne


def _mode(wave_mode="O"):
    traj = np.array([
        geodetic_to_ecef(30.0, 120.0, 60.0),
        geodetic_to_ecef(30.5, 120.0, 90.0),
        geodetic_to_ecef(31.0, 120.0, 60.0),
    ])
    return {
        "group_path_km": 220.0,
        "geom_path_km": 220.0,
        "tau_ms": 220.0 / cfg.C_KMS * 1e3,
        "beta_deg_tx": 20.0,
        "trajectory_ecef_km": traj,
        "earth_radius_km": cfg.RE_KM,
        "earth_model": "WGS84",
        "frequency_MHz": cfg.FREQ_MHZ,
        "branch_id": wave_mode,
        "b_ecef_hat": [0.0, 0.0, 1.0],
        "gyrofrequency_MHz": float(cfg.GEOMAG["fH_MHz"]),
        "wave_mode": wave_mode,
    }


def _ok(label, condition):
    print(" {} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run():
    n_fail = 0
    print("=" * 55)
    print(" update2.1 verification: AH complex absorption")
    print("=" * 55)

    path = [np.array([0.0, 0.0, 80.0]), np.array([100.0, 0.0, 80.0])]
    ne_sampler = lambda x, y, z: 1.0e10

    zero = complex_absorption_path_dB(
        path, cfg.FREQ_MHZ, ne_sampler,
        {"collision_model": "constant", "nu0_s": 0.0}, geomag=cfg.GEOMAG)
    n_fail += _ok("nu=0 gives zero loss", abs(zero["abs_total_dB"]) < 1e-12)

    one = complex_absorption_path_dB(
        path, cfg.FREQ_MHZ, ne_sampler,
        {"collision_model": "constant", "nu0_s": 5.0e6}, two_way=False,
        wave_mode="O", geomag=cfg.GEOMAG)
    two = complex_absorption_path_dB(
        path, cfg.FREQ_MHZ, ne_sampler,
        {"collision_model": "constant", "nu0_s": 5.0e6}, two_way=True,
        wave_mode="O", geomag=cfg.GEOMAG)
    n_fail += _ok("one-way positive", one["abs_total_dB"] > 0.0)
    n_fail += _ok("two-way doubles one-way", abs(two["abs_total_dB"] - 2.0 * one["abs_total_dB"]) < 1e-9)
    n_fail += _ok("AH index model recorded", one["abs_index_model"] == "appleton_hartree_collision")

    n_c = complex_refractive_index(cfg.FREQ_MHZ, 1.0e10, 5.0e6)
    n_fail += _ok("Y=0 complex index finite", np.isfinite(n_c.real) and np.isfinite(n_c.imag) and n_c.imag >= 0.0)
    n_ah_o = complex_refractive_index_ah(cfg.FREQ_MHZ, 1.0e10, 5.0e6, np.array([1.0, 0.0, 0.2]), wave_mode="O", geomag=cfg.GEOMAG)
    n_ah_x = complex_refractive_index_ah(cfg.FREQ_MHZ, 1.0e10, 5.0e6, np.array([1.0, 0.0, 0.2]), wave_mode="X", geomag=cfg.GEOMAG)
    n_fail += _ok("AH O/X complex indices finite", all(np.isfinite(v.real) and np.isfinite(v.imag) and v.imag >= 0.0 for v in (n_ah_o, n_ah_x)))

    old_abs = copy.deepcopy(cfg.ABSORPTION)
    old_complex = copy.deepcopy(cfg.UPDATE2["ABSORPTION_COMPLEX"])
    try:
        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update({"enable": False, "A0": 500.0, "chi_deg": 60.0})
        cfg.UPDATE2["ABSORPTION_COMPLEX"].clear()
        cfg.UPDATE2["ABSORPTION_COMPLEX"].update({"enable": False})
        clean = [_mode()]
        _add_power_3d(clean, n3d=_ConstN3D(1.0e10))

        cfg.UPDATE2["ABSORPTION_COMPLEX"].clear()
        cfg.UPDATE2["ABSORPTION_COMPLEX"].update({
            "enable": True,
            "collision_model": "constant",
            "nu0_s": 5.0e6,
        })
        comm = [_mode("O")]
        radar = [_mode("X")]
        _add_power_3d(comm, n3d=_ConstN3D(1.0e10), radar_mode=False)
        _add_power_3d(radar, n3d=_ConstN3D(1.0e10), radar_mode=True, sigma_rcs=1.0)
        unified = [_mode("O")]
        unified[0].update({"Pr_W": 1.0, "absorption_applied_once": True})
        _add_power_3d(unified, n3d=_ConstN3D(1.0e10))
        n_fail += _ok("explicit off no complex loss", clean[0]["abs_model"] == "none" and abs(clean[0]["abs_total_dB"]) < 1e-12)
        n_fail += _ok("enabled reduces comm power", comm[0]["Pr_W"] < clean[0]["Pr_W"] and comm[0]["abs_model"] == "appleton_hartree_collision_path")
        n_fail += _ok("radar absorption positive", radar[0]["abs_total_dB"] > 0.0)
        n_fail += _ok("unified absorption not overwritten", unified[0]["Pr_W"] == 1.0)
    finally:
        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update(old_abs)
        cfg.UPDATE2["ABSORPTION_COMPLEX"].clear()
        cfg.UPDATE2["ABSORPTION_COMPLEX"].update(old_complex)

    print("complex_absorption: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
