"""Verification update2.1/P3 local ray-frame utilities."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.local_frame import (
    local_to_global,
    make_boundary_field_descriptor,
    make_fdtd_adapter_descriptor,
    make_ray_local_frame,
    sample_Ne_local_cross_section,
)
from legacy.models.ray_tracer_3d import RefractiveIndex3D


def _straight_ray():
    states = []
    for s in np.linspace(0.0, 100.0, 11):
        states.append(np.array([s, 0.0, 100.0 + 0.2 * s, 1.0, 0.0, 0.2], dtype=float))
    return {"trajectory_3d": states, "Pr_W": 4.0, "wave_mode": "O"}


def _linear_field():
    x = np.linspace(-100.0, 200.0, 61)
    z = np.linspace(60.0, 180.0, 49)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    Ne = 1.0e10 + 1.0e8 * xx + 2.0e8 * zz
    return x, z, Ne


def _ok(label, cond):
    print("{} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def run() -> int:
    print("=" * 55)
    print(" update2.1/P3 verification: local ray frame")
    print("=" * 55)
    n_fail = 0

    ray = _straight_ray()
    frame = make_ray_local_frame(ray, sample_index=5)
    B = np.asarray(frame["basis_matrix"])
    ident = B @ B.T
    ok_orth = np.allclose(ident, np.eye(3), atol=1e-12)
    ok_right = np.dot(np.cross(frame["s_hat"], frame["u_hat"]), frame["v_hat"]) > 0.999
    n_fail += _ok("frame_orthonormal_right_handed", ok_orth and ok_right)

    origin = np.asarray(frame["origin_km"])
    shifted = local_to_global(frame, s_km=10.0, u_km=0.0)
    n_fail += _ok("local_to_global_s_shift", np.allclose(shifted - origin, frame["s_hat"] * 10.0))

    x, z, Ne = _linear_field()
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    cs = sample_Ne_local_cross_section(
        n3d, ray,
        s_offsets_km=np.array([-10.0, 0.0, 10.0]),
        u_offsets_km=np.array([-5.0, 0.0, 5.0]),
        sample_index=5)
    ok_sample = cs["Ne_m3"].shape == (3, 3)
    ok_sample = ok_sample and np.all(cs["Ne_m3"] >= 0.0)
    n_fail += _ok("cross_section_shape_nonnegative", ok_sample)

    bf = make_boundary_field_descriptor(ray, frame=frame, sample_index=5, Pr_W=4.0, k0_per_km=1.0, n3d=n3d)
    ok_bf = abs(bf["amplitude_sqrt_W"] - 2.0) < 1e-12
    ok_bf = ok_bf and bf["polarization"] == "O"
    ok_bf = ok_bf and "Ne_gradient_s_m3_per_km" in bf
    n_fail += _ok("boundary_field_contains_gradient", ok_bf)

    fdtd = make_fdtd_adapter_descriptor(ray, frame=frame)
    ok_fdtd = fdtd["fdtd_adapter_status"] == "descriptor_only"
    ok_fdtd = ok_fdtd and fdtd["field_to_pd"]["implemented"] is False
    ok_fdtd = ok_fdtd and fdtd["validity"] == "interface_only_not_executed"
    n_fail += _ok("fdtd_adapter_descriptor", ok_fdtd)

    print("local_frame: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
