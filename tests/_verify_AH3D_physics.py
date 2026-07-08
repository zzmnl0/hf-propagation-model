"""Physics checks for stage-1 AH3D Hamiltonian tracing."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, group_index_ah3d, trace_single_ray_ah_3d


def _synthetic_field():
    x = np.linspace(-20.0, 800.0, 83)
    z = np.linspace(60.0, 420.0, 91)
    zz = z[None, :]
    peak = 1.05e12 * np.exp(-((zz - 250.0) / 80.0) ** 2)
    return x, z, np.repeat(peak, x.size, axis=0)


def _rt_params():
    return {"ds_km": 1.0, "ds_min_km": 0.2, "ds_max_km": 2.0, "beta_min": 10.0, "beta_max": 65.0, "n_fan": 5, "z_stop_km": 420.0, "max_steps": 1200, "earth_flat": False}


def _ok(label, cond):
    print(" {} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def _check_group_index_isotropic():
    x, z, ne = _synthetic_field()
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    n_ah = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0})
    k = np.array([1.0, 0.0, 0.2])
    n = np.sqrt(n_ah.n2(200.0, 0.0, 220.0, k_vec=k))
    ng = group_index_ah3d(200.0, 0.0, 220.0, k, n_ah)
    return _ok("group_index_isotropic", abs(ng - 1.0 / n) < 0.02)


def _check_trace_residual_and_limit():
    x, z, ne = _synthetic_field()
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    iso = trace_single_ray_3d((0.0, 0.0), 30.0, 0.0, n3d, cfg.FREQ_MHZ, rt)
    n_ah = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0})
    ray = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, n_ah, cfg.FREQ_MHZ, rt)
    ok = ray["hamiltonian_residual_max"] < 1e-4
    ok = ok and abs(ray["tau_ms"] - iso["tau_ms"]) < 0.08
    ok = ok and ray["group_index_model"] == "frequency_derivative_appleton_hartree"
    return _ok("trace_residual_isotropic_limit", ok)


def _check_ox_split_status():
    x, z, ne = _synthetic_field()
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    ray_o = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag=cfg.GEOMAG), cfg.FREQ_MHZ, rt)
    ray_x = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "X", geomag=cfg.GEOMAG), cfg.FREQ_MHZ, rt)
    dtau = abs(ray_o["tau_ms"] - ray_x["tau_ms"])
    ok = np.isfinite(dtau) and dtau > 1e-5
    ok = ok and ray_o["branch_status"] in ("ok", "cutoff") and ray_x["branch_status"] in ("ok", "cutoff")
    return _ok("ox_split_status", ok)


def run():
    print("=" * 55)
    print(" stage1 P0 verification: AH3D physics")
    print("=" * 55)
    n_fail = _check_group_index_isotropic() + _check_trace_residual_and_limit() + _check_ox_split_status()
    print("AH3D_physics: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
