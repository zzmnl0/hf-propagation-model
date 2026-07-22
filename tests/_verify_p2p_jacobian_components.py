"""Focused A3R-P1 Jacobian check on the frozen Coleman rays."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d
from tests._benchmark_p2p_ah3d import _ColemanChapmanProvider


def _landing(ray):
    return np.array([ray["x_land_km"], ray["y_land_km"]], dtype=float)


def _check(beta_deg):
    provider = _ColemanChapmanProvider()
    ne = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z)
    n3d = RefractiveIndex3D(ne, cfg.BG_X, cfg.BG_Z, 10.0,
                           earth_flat=True, ne_provider=provider)
    n_ah = RefractiveIndexAH3D(
        n3d, 10.0, wave_mode="O",
        geomag={"fH_MHz": 0.0, "dip_deg": 0.0, "decl_deg": 0.0},
        integration_profile="smoke")

    variational_ray = trace_single_ray_ah_3d(
        cfg.TX_POS, beta_deg, 0.0, n_ah, 10.0, with_sensitivity=True)
    variational = np.asarray(variational_ray["landing_jacobian_deg"], dtype=float)
    references = []
    scales = [0.5, 1.0, 2.0]
    for scale in scales:
        columns = []
        for column, step in ((0, 0.05 * scale), (1, 0.10 * scale)):
            plus = [beta_deg, 0.0]
            minus = [beta_deg, 0.0]
            plus[column] += step
            minus[column] -= step
            rp = trace_single_ray_ah_3d(
                cfg.TX_POS, plus[0], plus[1], n_ah, 10.0)
            rm = trace_single_ray_ah_3d(
                cfg.TX_POS, minus[0], minus[1], n_ah, 10.0)
            columns.append((_landing(rp) - _landing(rm)) / (2.0 * step))
        references.append(np.column_stack(columns))

    changes = [np.linalg.norm(references[i] - references[i - 1]) /
               max(np.linalg.norm(references[i]), 1e-12)
               for i in range(1, len(references))]
    while min(changes) > cfg.P2P_3D["jacobian_scale_rtol"]:
        finer = 0.5 * scales[0]
        if finer < cfg.P2P_3D["jacobian_min_scale"] - 1e-15:
            break
        columns = []
        for column, step in ((0, 0.05 * finer), (1, 0.10 * finer)):
            plus = [beta_deg, 0.0]
            minus = [beta_deg, 0.0]
            plus[column] += step
            minus[column] -= step
            rp = trace_single_ray_ah_3d(
                cfg.TX_POS, plus[0], plus[1], n_ah, 10.0)
            rm = trace_single_ray_ah_3d(
                cfg.TX_POS, minus[0], minus[1], n_ah, 10.0)
            columns.append((_landing(rp) - _landing(rm)) / (2.0 * step))
        scales.insert(0, finer)
        references.insert(0, np.column_stack(columns))
        changes.insert(0, np.linalg.norm(references[1] - references[0]) /
                       max(np.linalg.norm(references[1]), 1e-12))
    reference = references[int(np.argmin(changes))]
    column_errors = [np.linalg.norm(variational[:, j] - reference[:, j]) /
                     max(np.linalg.norm(reference[:, j]), 1e-12)
                     for j in range(2)]
    return column_errors, changes, variational_ray


def run():
    failures = 0
    for name, beta in (("low_ray", 35.339808), ("high_ray", 39.918777)):
        columns, changes, ray = _check(beta)
        well_conditioned = float(ray["landing_jacobian_condition"]) <= 1e10
        selected = int(np.argmin(changes))
        stable = changes[selected] <= cfg.P2P_3D["jacobian_scale_rtol"]
        passed = (ray["branch_status"] == "ground" and well_conditioned
                  and np.all(np.isfinite(columns)) and stable)
        print("{} columns={} step_changes={} {}".format(
            name, [float(v) for v in columns], [float(v) for v in changes],
            "OK" if passed else "FAIL"))
        failures += int(not passed)
    print("{} fail(s)".format(failures))
    return failures


if __name__ == "__main__":
    sys.exit(run())
