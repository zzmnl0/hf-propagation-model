"""Smoke-profile V2 checks on every independent 3-D scenario."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d
from models.scenario_provider_factory import ScenarioProviderFactory


def run():
    factory = ScenarioProviderFactory(cfg)
    failures = []
    expected_under_resolved = {"spreadf", "es"}
    rt = {**cfg.RT_3D, "max_steps": 3000}
    for scenario in ("baseline", "tid", "spreadf", "bubble", "es"):
        provider = factory.create(scenario, seed=0, validation=True)
        ne = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z)
        n3d = RefractiveIndex3D(
            ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
            earth_flat=True, ne_provider=provider)
        model = RefractiveIndexAH3D(
            n3d, cfg.FREQ_MHZ, "O", geomag=cfg.GEOMAG,
            integration_profile="smoke")
        ray = trace_single_ray_ah_3d(
            cfg.TX_POS, 30.0, 0.0, model, cfg.FREQ_MHZ, rt)
        ok = (ray["ox_solver"] == "ah3d_v2_adaptive"
              and not ray["n2_floor_propagation"]
              and ray["integration_profile"] == "smoke")
        if scenario in expected_under_resolved:
            ok = ok and ray["under_resolved_Ne"]
        else:
            ok = ok and not ray["under_resolved_Ne"]
        print(f" {scenario}: {'OK' if ok else 'FAIL'} "
              f"status={ray['branch_status']} resolution={ray['resolution_check']}")
        if not ok:
            failures.append(scenario)
    print(f"AH3D V2 scenarios: {len(failures)} fail(s)")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(run())
