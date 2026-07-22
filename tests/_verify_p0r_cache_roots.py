"""Minimal P0R checks for exact caching and pre-formal root merging."""

import os
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import legacy.models.point_to_point_ah_3d as p2p


def _local(beta, az=0.0, x=1000.0, tau=5.0):
    return {
        "wave_mode": "O", "branch_status": "ground", "turning_point_count": 1,
        "beta_deg_tx": beta, "az_deg_tx": az, "x_land_km": x,
        "y_land_km": 0.0, "h_reflect_km": 250.0, "tau_ms": tau,
        "p2p_landing_error_km": abs(x - 1000.0),
    }


def run():
    calls = []
    original = p2p.trace_single_ray_ah_3d

    def fake_trace(*args, **kwargs):
        calls.append((args[1], args[2], kwargs.get("with_sensitivity", False)))
        return {"nfev": 4, "accepted_steps": 2}

    stats = {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {},
        "profile_accepted_steps": {}, "reference_ray_calls": 0,
        "reference_nfev": 0, "reference_accepted_steps": 0,
    }
    model = SimpleNamespace(integration_profile="smoke")
    cache = {}
    try:
        p2p.trace_single_ray_ah_3d = fake_trace
        p2p._trace_cached((0.0, 0.0), model, 10.0, 20.0, 0.0, 0.0,
                          False, "smoke", cache, stats)
        p2p._trace_cached((0.0, 0.0), model, 10.0, 20.0, 0.0, 0.0,
                          False, "smoke", cache, stats)
        p2p._trace_cached((0.0, 0.0), model, 10.0, 20.000001, 0.0, 0.0,
                          False, "smoke", cache, stats)
    finally:
        p2p.trace_single_ray_ah_3d = original

    duplicate = _local(20.0)
    duplicate_2 = _local(20.0 + 1e-5, x=1000.0001, tau=5.0000001)
    distinct = _local(20.001)
    groups = p2p._merge_local_roots([duplicate, duplicate_2, distinct], 0.05, 0.10)

    fan = [
        {"wave_mode": "O", "branch_status": "ground", "beta_deg_tx": beta,
         "az_deg_tx": 0.0, "x_land_km": x}
        for beta, x in ((5.0, 800.0), (20.0, 1200.0), (40.0, 900.0),
                        (60.0, 1250.0), (80.0, 700.0))
    ]
    candidate_diag = {}
    seeds, _ = p2p._candidate_seeds(
        (0.0, 0.0), (1000.0, 0.0), SimpleNamespace(wave_mode="O"),
        10.0, [], fan, 0.0, diagnostics=candidate_diag)
    checks = {
        "exact_cache_only": len(calls) == 2 and stats["cache_hits"] == 1,
        "duplicate_root_merged_once": len(groups) == 2 and len(groups[0]["members"]) == 2,
        "distinct_root_retained": len(groups[1]["members"]) == 1,
        "candidate_provenance": bool(seeds) and all("provenance" in seed for seed in seeds),
        "nearest_not_predeleted": candidate_diag["candidate_removed_as_covered"] == 0,
    }
    for name, ok in checks.items():
        print(f" {name}: {'OK' if ok else 'FAIL'}")
    failures = sum(not ok for ok in checks.values())
    print(f"P0R cache/root checks: {failures} fail(s)")
    return failures


if __name__ == "__main__":
    raise SystemExit(run())
