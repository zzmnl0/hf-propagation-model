"""Verification for update2.1 local ray-frame utilities."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.local_frame import local_to_global, make_ray_local_frame, sample_Ne_local_cross_section
from models.ray_tracer_3d import RefractiveIndex3D


def _straight_ray():
    states = []
    for s in np.linspace(0.0, 100.0, 11):
        states.append(np.array([s, 0.0, 100.0 + 0.2 * s, 1.0, 0.0, 0.2], dtype=float))
    return {"trajectory_3d": states}


def _linear_field():
    x = np.linspace(-100.0, 200.0, 61)
    z = np.linspace(60.0, 180.0, 49)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    Ne = 1.0e10 + 1.0e8 * xx + 2.0e8 * zz
    return x, z, Ne


def run() -> int:
    print("=" * 55)
    print(" update2.1 verification: local ray frame")
    print("=" * 55)
    n_fail = 0

    ray = _straight_ray()
    frame = make_ray_local_frame(ray, sample_index=5)
    B = np.asarray(frame["basis_matrix"])
    ident = B @ B.T
    ok_orth = np.allclose(ident, np.eye(3), atol=1e-12)
    ok_right = np.dot(np.cross(frame["s_hat"], frame["u_hat"]), frame["v_hat"]) > 0.999
    print(" frame_orthonormal {} right_handed {} {}".format(ok_orth, ok_right, "OK" if ok_orth and ok_right else "FAIL"))
    n_fail += 0 if ok_orth and ok_right else 1

    origin = np.asarray(frame["origin_km"])
    shifted = local_to_global(frame, s_km=10.0, u_km=0.0)
    ok_shift = np.allclose(shifted - origin, 10.0 * frame["s_hat"])
    print(" local_to_global_s_shift {}".format("OK" if ok_shift else "FAIL"))
    n_fail += 0 if ok_shift else 1

    x, z, Ne = _linear_field()
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    s_vals = np.array([-10.0, 0.0, 10.0])
    u_vals = np.array([-5.0, 0.0, 5.0])
    cs = sample_Ne_local_cross_section(n3d, ray, s_vals, u_vals, sample_index=5)
    ok_shape = cs["Ne_m3"].shape == (3, 3)
    x0 = cs["x_km"][1, 1]
    z0 = cs["z_km"][1, 1]
    expected = 1.0e10 + 1.0e8 * x0 + 2.0e8 * z0
    ok_center = abs(cs["Ne_m3"][1, 1] - expected) <= 1e-6 * expected
    ok_nonneg = np.all(cs["Ne_m3"] >= 0.0)
    ok_sample = ok_shape and ok_center and ok_nonneg
    print(" cross_section shape={} center={} nonneg={} {}".format(ok_shape, ok_center, ok_nonneg, "OK" if ok_sample else "FAIL"))
    n_fail += 0 if ok_sample else 1

    print("local_frame: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
