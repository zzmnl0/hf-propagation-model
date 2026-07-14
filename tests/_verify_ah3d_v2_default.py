"""Acceptance checks for replacing the AH3D legacy entry with V2."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.point_to_point_ah_3d import _shoot_ah_3d
from models.ray_tracer_3d import RefractiveIndex3D
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


def run():
    x = np.linspace(-50.0, 1500.0, 156)
    z = np.linspace(0.0, 500.0, 126)
    Ne = np.repeat(1.35e12 * np.exp(-((z[None, :] - 260.0) / 85.0) ** 2), x.size, 0)
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    geomag = {**cfg.GEOMAG, "fH_MHz": 0.0}
    default_model = RefractiveIndexAH3D(
        n3d, cfg.FREQ_MHZ, "O", geomag=geomag, integration_profile="smoke")
    legacy_model = RefractiveIndexAH3D(
        n3d, cfg.FREQ_MHZ, "O", geomag=geomag, use_dispersion_v2=False)
    rt = {**cfg.RT, "max_steps": 3000, "z_stop_km": 500.0}
    seed_ray = trace_single_ray_ah_3d(
        (0.0, 0.0), 28.0, 0.0, default_model, cfg.FREQ_MHZ, rt)
    p2p = _shoot_ah_3d(
        (0.0, 0.0), (seed_ray["x_land_km"], 0.0), default_model,
        cfg.FREQ_MHZ, 28.0, max_iter=2, tol_km=2.0, t=0.0)
    checks = {
        "default_solver_v2": cfg.UPDATE2["AH3D_SOLVER"] == "v2" and default_model.use_dispersion_v2,
        "legacy_explicit_only": not legacy_model.use_dispersion_v2,
        "fan_profile_screening": seed_ray["integration_profile"] == "smoke"
                                 and seed_ray["convergence_status"] == "screening_only",
        "v2_ground_turning": seed_ray["branch_status"] == "ground"
                             and seed_ray["turning_point_count"] == 1,
        "diagnostic_semantics": "rejected_steps" not in seed_ray
                                and "rejected_steps_estimate" in seed_ray
                                and "nfev" in seed_ray,
        "p2p_formal_recheck": p2p is not None
                              and p2p["convergence_status"] == "passed"
                              and p2p["p2p_search_profile"] == "smoke"
                              and p2p["p2p_final_profile"] == "research",
    }
    for name, ok in checks.items():
        print(f" {name}: {'OK' if ok else 'FAIL'}")
    failures = sum(not ok for ok in checks.values())
    print(f"AH3D V2 default: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
