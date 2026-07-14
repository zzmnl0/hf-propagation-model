"""A3 checks: terminal sensitivity and O/X product semantics."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.ray_tracer_3d import RefractiveIndex3D
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d
from models.point_to_point_ah_3d import _candidate_seeds


def run():
    x = np.linspace(-20.0, 900.0, 93)
    z = np.linspace(0.0, 420.0, 85)
    ne = np.repeat(1.2e12 * np.exp(-((z[None, :] - 250.0) / 80.0) ** 2), x.size, 0)
    n3d = RefractiveIndex3D(ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    model = RefractiveIndexAH3D(
        n3d, cfg.FREQ_MHZ, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0},
        integration_profile="smoke")
    rt = {**cfg.RT, "max_steps": 1500, "z_stop_km": 420.0}
    beta, az = 32.0, 0.0
    ray = trace_single_ray_ah_3d(
        (0.0, 0.0), beta, az, model, cfg.FREQ_MHZ, rt,
        with_sensitivity=True)
    h_beta, h_az = 0.02, 0.04

    def landing(b, a):
        r = trace_single_ray_ah_3d((0.0, 0.0), b, a, model, cfg.FREQ_MHZ, rt)
        return np.array([r["x_land_km"], r["y_land_km"]])

    reference = np.column_stack([
        (landing(beta + h_beta, az) - landing(beta - h_beta, az)) / (2 * h_beta),
        (landing(beta, az + h_az) - landing(beta, az - h_az)) / (2 * h_az),
    ])
    tangent = np.asarray(ray["landing_jacobian_deg"])
    relative = np.linalg.norm(tangent - reference) / max(np.linalg.norm(reference), 1e-12)
    fake_fan = [
        {"wave_mode": "O", "branch_status": "ground", "beta_deg_tx": beta_i,
         "az_deg_tx": 0.0, "x_land_km": x_land}
        for beta_i, x_land in ((5.0, 800.0), (20.0, 1200.0),
                               (40.0, 900.0), (60.0, 1250.0), (80.0, 700.0))
    ]
    candidate_seeds, audit_count = _candidate_seeds(
        (0.0, 0.0), (1000.0, 0.0), model, cfg.FREQ_MHZ, [], fake_fan, 0.0)
    checks = {
        "ground_event": ray["branch_status"] == "ground",
        "event_sensitivity": ray["landing_event_sensitivity_corrected"],
        "finite_jacobian": tangent.shape == (2, 2) and np.all(np.isfinite(tangent)),
        "smoke_central_agreement": relative <= 3e-2,
        "fan_default_o_only": cfg.UPDATE2["RAY_FAN_INCLUDE_X"] is False,
        "p2p_default_both_contract": cfg.P2P_3D["jacobian"] == "variational",
        "research_fallback_threshold": cfg.P2P_3D["jacobian_rtol"] == 1e-3,
        "full_domain_candidates": audit_count == len(fake_fan)
                                  and len(candidate_seeds) >= 3,
    }
    for name, ok in checks.items():
        print(f" {name}: {'OK' if ok else 'FAIL'}")
    print(f" variational_relative_error: {relative:.3e}")
    failures = sum(not ok for ok in checks.values())
    print(f"A3 variational P2P: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
