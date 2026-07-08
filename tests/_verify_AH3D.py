"""Verification for update2.1 AH3D ray tracing.

Run:
    conda run -n pytorch_cpu python tests/_verify_AH3D.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


def _synthetic_field():
    x = np.linspace(-50.0, 1500.0, 156)
    z = np.linspace(60.0, 500.0, 111)
    zz = z[None, :]
    peak = 1.35e12 * np.exp(-((zz - 260.0) / 85.0) ** 2)
    Ne = np.repeat(peak, x.size, axis=0)
    return x, z, Ne


def _rt_params():
    return {
        "ds_km": 1.0,
        "ds_min_km": 0.2,
        "ds_max_km": 3.0,
        "beta_min": 10.0,
        "beta_max": 70.0,
        "n_fan": 7,
        "z_stop_km": 500.0,
        "max_steps": 1800,
        "earth_flat": False,
    }


def _check_isotropic_limit() -> int:
    x, z, Ne = _synthetic_field()
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    iso = trace_single_ray_3d((0.0, 0.0), 28.0, 0.0, n3d, cfg.FREQ_MHZ, rt)
    failures = 0
    for mode in ("O", "X"):
        n_ah = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, mode, geomag={**cfg.GEOMAG, "fH_MHz": 0.0})
        ah = trace_single_ray_ah_3d((0.0, 0.0), 28.0, 0.0, n_ah, cfg.FREQ_MHZ, rt)
        dtau = abs(iso["tau_ms"] - ah["tau_ms"])
        dgp = abs(iso["group_path_km"] - ah["group_path_km"])
        ok = dtau < 0.05 and dgp < 10.0 and ah["branch_status"] in ("ok", "cutoff")
        fields_ok = all(k in ah for k in (
            "trajectory_3d",
            "geom_path_km",
            "group_path_km",
            "tau_ms",
            "h_reflect_km",
            "x_land_km",
            "y_land_km",
            "azimuth_deflect_deg",
            "beta_deg_tx",
            "az_deg_tx",
            "wave_mode",
            "branch_status",
            "group_index_model",
        ))
        ok = ok and fields_ok
        print(" isotropic_limit_{} dtau={:.6f} dgp={:.3f} {}".format(mode, dtau, dgp, "OK" if ok else "FAIL"))
        failures += 0 if ok else 1
    return failures


def _check_weak_magnetic_split() -> int:
    x, z, Ne = _synthetic_field()
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    geomag = {**cfg.GEOMAG, "fH_MHz": 0.05}
    n_o = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag=geomag)
    n_x = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "X", geomag=geomag)
    ray_o = trace_single_ray_ah_3d((0.0, 0.0), 28.0, 0.0, n_o, cfg.FREQ_MHZ, rt)
    ray_x = trace_single_ray_ah_3d((0.0, 0.0), 28.0, 0.0, n_x, cfg.FREQ_MHZ, rt)
    dtau = abs(ray_o["tau_ms"] - ray_x["tau_ms"])
    ok = np.isfinite(dtau) and 0.0 <= dtau < 0.2
    ok = ok and ray_o["wave_mode"] == "O" and ray_x["wave_mode"] == "X"
    print(" weak_magnetic_split dtau={:.6f} {}".format(dtau, "OK" if ok else "FAIL"))
    return 0 if ok else 1


def _check_near_cutoff() -> int:
    x, z, Ne = _synthetic_field()
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    n_x = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "X", geomag=cfg.GEOMAG)
    vals = [
        n_x.n2(500.0, 0.0, 220.0, k_vec=np.array([1.0, 0.0, 0.2])),
        n_x.n2(500.0, 0.0, 260.0, k_vec=np.array([1.0, 0.0, 0.2])),
        n_x.n2(500.0, 0.0, 300.0, k_vec=np.array([1.0, 0.0, 0.2])),
    ]
    ok = all(np.isfinite(v) and v >= n_x.n2_floor for v in vals)
    status = n_x.branch_status(500.0, 0.0, 260.0, k_vec=np.array([1.0, 0.0, 0.2]))
    ok = ok and status in ("ok", "cutoff")
    print(" near_cutoff status={} {}".format(status, "OK" if ok else "FAIL"))
    return 0 if ok else 1


def run() -> int:
    print("=" * 55)
    print(" update2.1 verification: AH3D")
    print("=" * 55)
    n_fail = 0
    n_fail += _check_isotropic_limit()
    n_fail += _check_weak_magnetic_split()
    n_fail += _check_near_cutoff()
    print("AH3D: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
