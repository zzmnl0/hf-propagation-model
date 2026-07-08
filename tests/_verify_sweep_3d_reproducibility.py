"""Verification for stage-1 compact 3-D parameter sweep."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.parameter_sweep_3d import run_sweep


def run():
    print("=" * 55)
    print(" stage1 P0 verification: 3-D parameter sweep")
    print("=" * 55)
    out = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) / "output"
    rows1, path1 = run_sweep(output_dir=out)
    rows2, path2 = run_sweep(output_dir=out)
    n_fail = 0
    ok = len(rows1) == len(rows2) and path1.exists() and path2.exists()
    for a, b in zip(rows1, rows2):
        for key in ("x_land_km", "y_land_km", "tau_ms", "h_reflect_km", "max_relative_perturbation"):
            ok = ok and np.isclose(float(a[key]), float(b[key]), rtol=0.0, atol=1e-10)
    print(" sweep_reproducibility {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1
    print("sweep_3d_reproducibility: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
