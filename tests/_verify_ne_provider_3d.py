"""Verification for stage-1 3-D NeProvider."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.ne_provider import GridNeProvider, GridNeProvider3D, ParameterizedIonosphere3DProvider
from legacy.models.ray_tracer_3d import RefractiveIndex3D


def _field():
    x = np.linspace(0.0, 100.0, 6)
    y = np.linspace(-40.0, 40.0, 5)
    z = np.linspace(180.0, 320.0, 8)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    ne = 8.0e10 + 1.0e8 * xx + 2.0e8 * yy + 3.0e8 * zz
    return x, y, z, ne


def _ok(label, cond):
    print(" {} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def _check_grid3d():
    x, y, z, ne = _field()
    provider = GridNeProvider3D(ne, x, y, z)
    val = provider.get_Ne((50.0, 0.0, 250.0))
    grad = provider.get_grad_Ne((50.0, 0.0, 250.0), step_km=2.0)
    grid = provider.sample_grid3d([0.0, 100.0], [-40.0, 40.0], [180.0, 320.0])
    meta = provider.get_metadata()
    ok = val > 0.0 and grid.shape == (2, 2, 2) and meta["dimensions"] == "3d"
    ok = ok and np.allclose(grad, [1.0e8, 2.0e8, 3.0e8], rtol=0.02)
    return _ok("grid3d_interface", ok)


def _check_parameterized():
    x = np.linspace(0.0, 200.0, 9)
    z = np.linspace(180.0, 340.0, 9)
    bg = np.full((x.size, z.size), 1.0e11)
    base = GridNeProvider(bg, x, z)
    tid = {"enable": True, "n_components": 1, "amplitude_list": [0.1], "az_deg_list": [35.0], "period_s_list": [2400.0], "lambda_h_km_list": [300.0], "z0_km": 260.0, "sigma_z_km": 70.0}
    bubble = {"enable": True, "delta_max": 0.4, "x0_km": 100.0, "y0_km": 0.0, "z0_km": 260.0, "sigma_perp_km": 40.0, "sigma_para_km": 500.0}
    provider = ParameterizedIonosphere3DProvider(base, tid_params=tid, bubble_params=bubble, geomag={**cfg.GEOMAG, "link_bearing_deg": 0.0}, seed=7)
    off = ParameterizedIonosphere3DProvider(base, tid_params={"enable": False}, bubble_params={"enable": False}, geomag=cfg.GEOMAG)
    pts = np.array([[0.0, 0.0, 220.0], [100.0, 0.0, 260.0], [150.0, 30.0, 280.0]])
    vals = provider.sample_ray_path(pts)
    vals2 = ParameterizedIonosphere3DProvider(base, tid_params=tid, bubble_params=bubble, geomag={**cfg.GEOMAG, "link_bearing_deg": 0.0}, seed=7).sample_ray_path(pts)
    summary = provider.perturbation_summary([0.0, 100.0, 200.0], [-30.0, 0.0, 30.0], [220.0, 260.0, 300.0])
    ok = np.all(vals >= 0.0) and np.allclose(vals, vals2)
    ok = ok and abs(off.get_Ne((100.0, 20.0, 260.0)) - base.get_Ne((100.0, 20.0, 260.0))) < 1.0
    ok = ok and summary["max_relative_perturbation"] > 0.05
    return _ok("parameterized_3d", ok)


def _check_refractive_provider3d():
    x, y, z, ne = _field()
    provider = GridNeProvider3D(ne, x, y, z)
    n3d = RefractiveIndex3D(provider.get_Ne_grid(x, z), x, z, cfg.FREQ_MHZ, earth_flat=False, ne_provider=provider)
    a = n3d.n2(50.0, -20.0, 250.0)
    b = n3d.n2(50.0, 20.0, 250.0)
    ok = np.isfinite(a) and np.isfinite(b) and abs(a - b) > 1e-8
    return _ok("refractive_provider3d_y_dependence", ok)


def run():
    print("=" * 55)
    print(" stage1 P0 verification: 3-D NeProvider")
    print("=" * 55)
    n_fail = _check_grid3d() + _check_parameterized() + _check_refractive_provider3d()
    print("ne_provider_3d: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
