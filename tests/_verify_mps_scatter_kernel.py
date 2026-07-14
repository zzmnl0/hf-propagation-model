"""Verify P3 MPS scatter-kernel semantics."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.pd_scatter_kernel import apply_scatter_kernel
from models.scatter_phase_screen import make_path_mps_delay_kernel


class ConstN3D:
    def Ne_at(self, _x, _y, _z, t=0.0):
        return 1.0e11


def _mode():
    pts = []
    for x in np.linspace(0.0, 120.0, 25):
        pts.append(np.array([x, 0.0, 240.0, 1.0, 0.0, 0.0], dtype=float))
    return {
        "label": "m1",
        "tau_ms": 3.0,
        "Pr_W": 2.0,
        "trajectory_3d": pts,
        "scatter_DeltaL_km": 60.0,
    }


def _ok(label, cond):
    print("{} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def run():
    tau = np.linspace(2.5, 3.5, 101)
    phase = {"N_screen": 2, "n_pts": 32, "dx_m": 500.0, "seed": 7}
    scatter = {**cfg.SCATTER, "Cs_rel": 0.01, "DeltaL_km": 60.0}
    zero = {**cfg.SCATTER, "Cs_rel": 0.0, "sigma_phi_rad": 0.0, "DeltaL_km": 60.0}
    n_fail = 0

    r1 = make_path_mps_delay_kernel(_mode(), ConstN3D(), tau, scatter, phase, cfg.FREQ_MHZ)
    r2 = make_path_mps_delay_kernel(_mode(), ConstN3D(), tau, scatter, phase, cfg.FREQ_MHZ)
    ok = np.allclose(r1["mps_delay_kernel"], r2["mps_delay_kernel"])
    ok = ok and abs(float(np.sum(r1["mps_delay_kernel"])) - 1.0) < 1e-9
    ok = ok and r1["scatter_validity"] == "phase_screen_not_fullwave"
    n_fail += _ok("mps_seed_reproducible_normalized", ok)

    rz = make_path_mps_delay_kernel(_mode(), ConstN3D(), tau, zero, phase, cfg.FREQ_MHZ)
    ok = rz["mps_status"] == "sigma_phi_zero_delta"
    ok = ok and np.count_nonzero(rz["mps_delay_kernel"]) == 1
    ok = ok and abs(float(np.sum(rz["mps_delay_kernel"])) - 1.0) < 1e-12
    n_fail += _ok("mps_sigma_zero_delta", ok)

    mode = {"tau_ms": 3.0, "Pr_W": 2.0, "scatter_method": "MPS_required", "sigma_tau_ms": 0.3}
    _ta, pd = apply_scatter_kernel([mode], tau, tau_res_ms=0.01)
    ok = mode["pd_kernel_source"] == "mps_required_no_weak_formula"
    ok = ok and abs(float(np.sum(pd)) - 2.0) < 1e-9
    ok = ok and np.count_nonzero(pd) == 1
    n_fail += _ok("mps_required_no_weak_formula", ok)

    print("mps_scatter_kernel: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
