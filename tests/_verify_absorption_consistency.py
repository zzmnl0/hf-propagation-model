"""Verify P1 AH-collisional absorption consistency."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.absorption_collision import complex_absorption_path_dB


def _sampler(_x, _y, _z):
    return 1.0e10


def _ok(label, condition):
    print("{} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run() -> int:
    path = [np.array([0.0, 0.0, 80.0]), np.array([100.0, 0.0, 80.0])]
    base = {
        "collision_model": "constant",
        "nu0_s": 5.0e6,
    }
    n_fail = 0
    z0 = complex_absorption_path_dB(path, 10.0, _sampler, {**base, "nu0_s": 0.0})
    one = complex_absorption_path_dB(path, 10.0, _sampler, base, two_way=False)
    two = complex_absorption_path_dB(path, 10.0, _sampler, base, two_way=True)
    high = complex_absorption_path_dB(path, 20.0, _sampler, base, two_way=False)
    ox_o = complex_absorption_path_dB(path, 10.0, _sampler, base, wave_mode="O")
    ox_x = complex_absorption_path_dB(path, 10.0, _sampler, base, wave_mode="X")
    n_fail += _ok("nu_zero_abs_zero", abs(z0["abs_total_dB"]) < 1e-12)
    n_fail += _ok("two_way_double", abs(two["abs_total_dB"] - 2.0 * one["abs_total_dB"]) <= max(1e-12, one["abs_total_dB"] * 1e-9))
    n_fail += _ok("higher_freq_lower_abs", high["abs_total_dB"] < one["abs_total_dB"])
    n_fail += _ok("ox_finite_nonnegative", all(np.isfinite(v["abs_total_dB"]) and v["abs_total_dB"] >= 0.0 for v in (ox_o, ox_x)))
    print("{} fail(s)".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
