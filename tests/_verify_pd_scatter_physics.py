"""Verify P1 scatter-to-P-D physical kernel behavior."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.pd_scatter_kernel import apply_scatter_kernel, make_delay_kernel


def _ok(label, condition):
    print("{} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run() -> int:
    tau = np.arange(0.0, 10.0, 0.02)
    n_fail = 0
    delta = make_delay_kernel({"tau_ms": 4.0, "Pr_W": 2.0, "scatter_method": "Rytov", "S4": 0.0}, tau)
    n_fail += _ok("s4_zero_delta", np.count_nonzero(delta) == 1 and abs(float(np.sum(delta)) - 1.0) < 1e-12)

    modes = [{"tau_ms": 4.0, "Pr_W": 2.0, "scatter_method": "Born", "sigma_phi_sq": 0.25, "scatter_DeltaL_km": 120.0}]
    _, pd = apply_scatter_kernel(modes, tau)
    n_fail += _ok("born_power_conserved", abs(float(np.sum(pd)) - 2.0) < 1e-10 and float(np.max(pd)) < 2.0)

    real_kernel = np.zeros_like(tau)
    real_kernel[100:105] = 1.0
    _, pd_mps = apply_scatter_kernel([{"tau_ms": 4.0, "Pr_W": 3.0, "scatter_method": "MPS", "mps_delay_kernel": real_kernel}], tau)
    n_fail += _ok("mps_real_kernel_only", abs(float(np.sum(pd_mps)) - 3.0) < 1e-10 and np.count_nonzero(pd_mps) == 5)

    _, pd_req = apply_scatter_kernel([{"tau_ms": 4.0, "Pr_W": 3.0, "scatter_method": "MPS_required", "delta_tau_ms": 0.02}], tau)
    n_fail += _ok("mps_required_no_fake_kernel", abs(float(np.sum(pd_req)) - 3.0) < 1e-10 and np.count_nonzero(pd_req) > 5)
    print("{} fail(s)".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
