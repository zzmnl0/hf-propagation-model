"""Fast structural regression for the frozen P0C high/low reference."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a3r_benchmark"


def _fold_count(samples):
    x = np.asarray(samples, dtype=float)
    slopes = np.diff(x[:, 1]) / np.diff(x[:, 0])
    return int(np.count_nonzero(slopes[:-1] * slopes[1:] < 0.0))


def main() -> int:
    prep = json.loads((OUT / "a3r_p0c_reference_preparation.json").read_text(encoding="utf-8"))
    closure = json.loads((OUT / "a3r_p0c_closure.json").read_text(encoding="utf-8"))
    refined = prep["fold_candidates"][0]["refined_samples"]
    assert _fold_count([(v["beta_deg"], v["x_land_km"]) for v in refined]) == 1
    assert _fold_count([(2.0, 100.0), (3.0, 200.0), (4.0, 300.0)]) == 0
    assert prep["branch_semantics"] == "magnetoionic_degenerate_single_count"
    assert all(len(run["production"]) == 2 for run in closure["repeats"])
    acceptance = closure["acceptance"]
    assert acceptance["p0c_3_p0_reference_complete"] is True
    assert acceptance["p1_gate"] == "open"
    assert acceptance["p2_full_azimuth_audit_required"] is True
    assert acceptance["magnetoionic_caustic_validation_required"] is True
    print("P0C high/low reference: 0 fail(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
