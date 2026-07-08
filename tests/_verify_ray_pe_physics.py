"""Verify P1 Ray-PE/SSF physical diagnostics."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.hybrid_ray_pe import apply_ray_pe_coupling, apply_ray_pe_coupling_3d
from models.pe_propagator import PEPropagator


def _grid():
    x = np.linspace(0.0, 40.0, 41)
    z = np.linspace(200.0, 260.0, 121)
    return np.zeros((len(x), len(z))), x, z, {"z0_km": 230.0, "Lz_km": 20.0, "Lx_km": 15.0}


def _params(dx, dz_m):
    p = dict(cfg.PE)
    p.update({"dx_km": dx, "dz_m": dz_m, "n_pml": 0, "store_history": False, "earth_flat": False, "w0_km": 8.0})
    return p


def _mode():
    states = []
    for x in np.linspace(5.0, 35.0, 9):
        states.append(np.array([x, 0.0, 230.0 + 0.01 * (x - 20.0) ** 2, 1.0, 0.0, 0.0]))
    return {"trajectory_3d": states, "wave_mode": "O"}


def _ok(label, condition):
    print("{} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run() -> int:
    Ne, x, z, bubble = _grid()
    Pr = 1e-9
    n_fail = 0
    pe = PEPropagator(cfg.FREQ_MHZ, _params(1.0, 500.0))
    Pr_diag, diag = apply_ray_pe_coupling(pe, Ne, x, z, Pr, pe.k0, 5.0, 230.0, 0.0, bubble, {"enable": True, "apply_to_power": False})
    n_fail += _ok("uniform_diag_power_unchanged", abs(Pr_diag - Pr) <= 1e-21 and 0.98 <= diag["pe_power_ratio"] <= 1.02)

    pe_fine = PEPropagator(cfg.FREQ_MHZ, _params(0.5, 250.0))
    _, d_fine = apply_ray_pe_coupling(pe_fine, Ne, x, z, Pr, pe_fine.k0, 5.0, 230.0, 0.0, bubble, {"enable": True, "apply_to_power": False})
    n_fail += _ok("dx_dz_halving_converges", abs(d_fine["pe_power_ratio"] - diag["pe_power_ratio"]) < 0.05)

    irregularity = {**cfg.BUBBLE_3D, "enable": True, "x0_km": 20.0, "z0_km": 230.0, "sigma_perp_km": 50.0, "sigma_para_km": 300.0}
    _, d3 = apply_ray_pe_coupling_3d(pe, _mode(), Ne, x, z, Pr, pe.k0, irregularity, {"enable": True, "apply_to_power": False})
    fields = ["pe_boundary_amplitude_sqrt_W", "pe_boundary_phase_rad", "pe_boundary_curvature_km_inv", "pe_boundary_kx", "pe_boundary_polarization"]
    n_fail += _ok("boundary_field_complete", d3.get("pe_frame_status") == "local_2p5d" and all(k in d3 for k in fields) and d3["pe_boundary_amplitude_sqrt_W"] > 0.0)
    print("{} fail(s)".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
