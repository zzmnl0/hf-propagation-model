"""Fast control-flow checks for the production multiscale P2P Jacobian."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models import point_to_point_ah_3d as p2p


class _Model:
    wave_mode = "O"
    integration_profile = "smoke"


def _stats():
    return {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {},
        "profile_accepted_steps": {}, "reference_ray_calls": 0,
        "reference_nfev": 0, "reference_accepted_steps": 0,
        "jacobian_rejected_local_roots": {},
    }


def _ray(beta, az, branch_status="ground"):
    return {
        "x_land_km": 10.0 * beta,
        "y_land_km": 10.0 * az,
        "branch_status": branch_status,
        "nfev": 1,
        "accepted_steps": 1,
    }


def run():
    original = p2p._trace_cached
    try:
        p2p._trace_cached = lambda tx, model, freq, beta, az, t, sensitivity, profile, cache, stats: _ray(beta, az)
        accepted = p2p._shoot_ah_3d(
            (0.0, 0.0), (400.0, 0.0), _Model(), 10.0, 40.0,
            defer_formal=True, trace_stats=_stats())

        rejected_stats = _stats()

        def invalid_perturbation(tx, model, freq, beta, az, t, sensitivity,
                                 profile, cache, stats):
            status = "ground" if beta == 40.0 and az == 0.0 else "evanescent"
            return _ray(beta, az, status)

        p2p._trace_cached = invalid_perturbation
        rejected = p2p._shoot_ah_3d(
            (0.0, 0.0), (400.0, 0.0), _Model(), 10.0, 40.0,
            defer_formal=True, trace_stats=rejected_stats)
    finally:
        p2p._trace_cached = original

    checks = {
        "production_default": cfg.P2P_3D["jacobian"] == "central_reference_multiscale",
        "stable_reference_accepted": accepted is not None and accepted.get(
            "p2p_jacobian_effective") == "central_reference_multiscale",
        "invalid_branch_rejected": rejected is None and rejected_stats[
            "jacobian_rejected_local_roots"].get(
                "jacobian_perturbation_not_ground") == 1,
    }
    for name, passed in checks.items():
        print("{}: {}".format(name, "OK" if passed else "FAIL"))
    failures = sum(not passed for passed in checks.values())
    print("P2P multiscale reference: {} fail(s)".format(failures))
    return failures


if __name__ == "__main__":
    raise SystemExit(run())
