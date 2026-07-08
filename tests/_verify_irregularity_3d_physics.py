"""Verify P1 3-D irregularity physics diagnostics."""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.irregularity_3d import (
    bragg_mismatch,
    bubble_descriptor,
    es_descriptor,
    geomag_unit_vector,
    spreadf_descriptor,
    summarize_irregularities,
    tid_descriptor,
)


def _ok(label, condition):
    print("{} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run() -> int:
    n_fail = 0
    b = geomag_unit_vector()
    n_fail += _ok("geomag_unit", abs(float(np.linalg.norm(b)) - 1.0) < 1e-12)
    bm = bragg_mismatch(10.0, [1, 0, 0], [1, 0, 0], 0.015)
    n_fail += _ok("bragg_fields_finite", all(np.isfinite(v) for v in bm.values()))
    tid = tid_descriptor({"enable": True, "n_components": 2, "amplitude_list": [0.1, 0.05], "lambda_h_km_list": [300, 250], "az_deg_list": [0, 45]}, 10.0)
    n_fail += _ok("tid_components", tid["enabled"] and len(tid["components"]) == 2 and tid["scatter_DeltaL_km"] > 0)
    bub = bubble_descriptor({"enable": True, "delta_max": 0.5, "sigma_perp_km": 50, "sigma_para_km": 500}, 10.0)
    n_fail += _ok("bubble_anisotropy", bub["enabled"] and bub["anisotropy_ratio"] == 10.0 and len(bub["b_hat"]) == 3)
    sf = spreadf_descriptor({"enable": True, "Cs": 1e-3, "p": 3.0, "L0_km": 50.0})
    es = es_descriptor({"enable": True, "delta_h_m": 100.0}, {"L_perp_m": 300.0, "L_para_m": 9000.0, "p_es": 1.5})
    n_fail += _ok("spreadf_es_positive_scales", sf["outer_scale_km"] > 0 and es["anisotropy_ratio"] > 1.0)
    summary = summarize_irregularities({"TID": {"enable": False}, "BUBBLE_3D": {"enable": False}}, 10.0)
    n_fail += _ok("summary_keys", all(k in summary for k in ("TID", "bubble", "Spread-F", "Es")))
    print("{} fail(s)".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
