"""A2 adaptive AH3D integration acceptance checks."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.ray_tracer_3d import RefractiveIndex3D
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


def run():
    x = np.linspace(-50.0, 1500.0, 156)
    z = np.linspace(0.0, 500.0, 126)
    Ne = np.repeat(1.35e12 * np.exp(-((z[None, :] - 260.0) / 85.0) ** 2), x.size, 0)
    n3d = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    n_ah = RefractiveIndexAH3D(
        n3d, cfg.FREQ_MHZ, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0},
        use_dispersion_v2=True,
    )
    rt = {**cfg.RT, "max_steps": 3000, "z_stop_km": 500.0}
    ray = trace_single_ray_ah_3d((0.0, 0.0), 28.0, 0.0, n_ah, cfg.FREQ_MHZ, rt)
    checks = {
        "ground_after_turning": ray["branch_status"] == "ground" and ray["turning_point_count"] == 1,
        "dispersion_constraint": ray["hamiltonian_residual_max"] <= 1e-6,
        "landing_convergence": ray["convergence_landing_delta_km"] <= 1.0,
        "height_convergence": ray["convergence_reflection_height_delta_km"] <= 0.5,
        "delay_convergence": ray["convergence_group_delay_relative"] <= 1e-3,
        "research_gate": ray["convergence_status"] == "passed",
        "no_floor": ray["n2_floor_propagation"] is False,
    }
    for name, ok in checks.items():
        print(f" {name}: {'OK' if ok else 'FAIL'}")
    failures = sum(not ok for ok in checks.values())
    print(f"AH3D A2 convergence: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
