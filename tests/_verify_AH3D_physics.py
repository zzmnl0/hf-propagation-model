"""Physics checks for stage-1 AH3D Hamiltonian tracing."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from utils import freq_to_k0
from legacy.models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D, group_index_ah3d, trace_single_ray_ah_3d


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
    direction = np.array([1.0, 0.0, 0.2])
    direction /= np.linalg.norm(direction)
    n = np.sqrt(n_ah.n2(200.0, 0.0, 180.0, k_vec=direction))
    k = freq_to_k0(cfg.FREQ_MHZ) * n * direction
    ng = group_index_ah3d(200.0, 0.0, 180.0, k, n_ah)
    return _ok("group_index_isotropic", abs(ng - 1.0 / n) < 0.02)


def _check_trace_residual_and_limit():
    x, z, ne = _synthetic_field()
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    iso = trace_single_ray_3d((0.0, 0.0), 30.0, 0.0, n3d, cfg.FREQ_MHZ, rt)
    n_ah = RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0}, integration_profile="smoke")
    ray = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, n_ah, cfg.FREQ_MHZ, rt)
    ok = ray["hamiltonian_residual_max"] < 1e-4
    ok = ok and abs(ray["tau_ms"] - iso["tau_ms"]) < 0.12
    ok = ok and ray["group_index_model"] == "dispersion_gradient_-Dk_Domega"
    ok = ok and ray["group_path_definition"] == "c_integral_dt"
    ok = ok and ray["branch_status"] == "ground" and ray["turning_point_count"] == 1
    return _ok("trace_residual_isotropic_limit", ok)


def _check_ox_split_status():
    x, z, ne = _synthetic_field()
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    rt = _rt_params()
    ray_o = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "O", geomag=cfg.GEOMAG, integration_profile="smoke"), cfg.FREQ_MHZ, rt)
    ray_x = trace_single_ray_ah_3d((0.0, 0.0), 30.0, 0.0, RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, "X", geomag=cfg.GEOMAG, integration_profile="smoke"), cfg.FREQ_MHZ, rt)
    dtau = abs(ray_o["tau_ms"] - ray_x["tau_ms"])
    ok = np.isfinite(dtau) and dtau > 1e-5
    ok = ok and ray_o["branch_status"] == "ground"
    ok = ok and ray_x["branch_status"] == "ground"
    return _ok("ox_split_status", ok)


def _check_branch_status_raw_n2():
    f = cfg.FREQ_MHZ
    ne_cut = (f * 1e6 / cfg.K_FP) ** 2
    x = np.linspace(0.0, 10.0, 3)
    z = np.linspace(100.0, 120.0, 3)
    near = np.full((x.size, z.size), 0.995 * ne_cut)
    evan = np.full((x.size, z.size), 1.005 * ne_cut)
    n_near = RefractiveIndexAH3D(
        RefractiveIndex3D(near, x, z, f, earth_flat=False),
        f, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0})
    n_evan = RefractiveIndexAH3D(
        RefractiveIndex3D(evan, x, z, f, earth_flat=False),
        f, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0})
    ok = n_near.branch_status(5.0, 0.0, 110.0, k_vec=[1.0, 0.0, 0.0]) == "near_cutoff"
    ok = ok and n_evan.branch_status(5.0, 0.0, 110.0, k_vec=[1.0, 0.0, 0.0]) == "evanescent"
    ok = ok and n_evan.n2_raw(5.0, 0.0, 110.0, k_vec=[1.0, 0.0, 0.0]) < 0.0
    ok = ok and n_evan.n2(5.0, 0.0, 110.0, k_vec=[1.0, 0.0, 0.0]) < 0.0
    return _ok("branch_status_raw_n2", ok)


def run():
    print("=" * 55)
    print(" stage1 P0 verification: AH3D physics")
    print("=" * 55)
    n_fail = (
        _check_group_index_isotropic()
        + _check_trace_residual_and_limit()
        + _check_ox_split_status()
        + _check_branch_status_raw_n2()
    )
    print("AH3D_physics: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
