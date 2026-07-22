"""Verify update2.1 Ray-PE/SSF energy coupling.

Run:
    conda run -n pytorch_cpu python tests/_verify_ray_pe_energy.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.hybrid_ray_pe import apply_ray_pe_coupling
from legacy.models.pe_propagator import PEPropagator


def _local_grid():
    x = np.linspace(0.0, 40.0, 41)
    z = np.linspace(200.0, 260.0, 121)
    Ne = np.zeros((len(x), len(z)), dtype=float)
    bubble = {
        "z0_km": 230.0,
        "Lz_km": 20.0,
        "Lx_km": 15.0,
    }
    return Ne, x, z, bubble


def _pe_params(n_pml):
    params = dict(cfg.PE)
    params.update(
        {
            "dx_km": 1.0,
            "dz_m": 500.0,
            "n_pml": n_pml,
            "sigma_pml": 0.8,
            "w0_km": 8.0,
            "store_history": False,
            "earth_flat": False,
        }
    )
    return params


def _run_case(name, apply_to_power, n_pml, expect_loss=False):
    Ne, x, z, bubble = _local_grid()
    pe = PEPropagator(freq_MHz=cfg.FREQ_MHZ, pe_params=_pe_params(n_pml))
    Pr_in = 1.0e-9
    Pr_out, diag = apply_ray_pe_coupling(
        pe,
        Ne,
        x,
        z,
        Pr_in,
        pe.k0,
        x_entry_km=5.0,
        z_entry_km=230.0,
        beta_inc_deg=0.0,
        bubble_params=bubble,
        coupling_params={
            "enable": True,
            "apply_to_power": apply_to_power,
            "mode": "local_2p5d",
        },
    )
    ratio = diag["pe_power_ratio"]
    expected_Pr = Pr_in * ratio if apply_to_power else Pr_in
    ok_energy = abs(Pr_out - expected_Pr) <= max(1.0e-30, abs(expected_Pr) * 1.0e-12)
    ok_db = abs(diag["pe_Pr_correction_dB"] - 10.0 * np.log10(max(ratio, 1.0e-300))) < 1.0e-12
    ok_diag = diag["pe_coupled"] and diag["pe_status"] == "ok"
    ok_loss = (ratio < 0.999) if expect_loss else (0.98 <= ratio <= 1.02)
    ok = ok_energy and ok_db and ok_diag and ok_loss
    print(
        "{} ratio={:.6f} corr={:.3f}dB Pr_in={:.3e} Pr_out={:.3e} {}".format(
            name,
            ratio,
            diag["pe_Pr_correction_dB"],
            Pr_in,
            Pr_out,
            "OK" if ok else "FAIL",
        )
    )
    return 0 if ok else 1


def run() -> int:
    print("=" * 60)
    print("update2.1 Ray-PE/SSF energy coupling verification")
    print("=" * 60)
    n_fail = 0

    default_ok = (cfg.UPDATE2["PE_COUPLING"].get("enable") is False and
                  cfg.UPDATE2["PE_COUPLING"].get("apply_to_power") is False)
    print("default PE_COUPLING hybrid_off {}".format("OK" if default_ok else "FAIL"))
    n_fail += 0 if default_ok else 1

    n_fail += _run_case("diagnostic_only", False, n_pml=0, expect_loss=False)
    n_fail += _run_case("apply_to_power", True, n_pml=8, expect_loss=True)

    print("-" * 60)
    print("{} fail(s)".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())


