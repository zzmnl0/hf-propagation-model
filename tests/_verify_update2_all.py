"""Update2.1 phase-1 verification."""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from legacy.models.ne_provider import GridNeProvider
from legacy.models.ray_tracer_3d import RefractiveIndex3D


def _check_update2_defaults() -> int:
    u2 = cfg.UPDATE2
    checks = [
        ("AH3D_ENABLE", u2["AH3D_ENABLE"] is False),
        ("AH3D_SOLVER", u2["AH3D_SOLVER"] == "v2"),
        ("OX_MODE", u2["OX_MODE"] == "ah3d"),
        ("EXTERNAL_NE_PROVIDER", u2["EXTERNAL_NE_PROVIDER"]["enable"] is False),
        ("PE_COUPLING.enable", u2["PE_COUPLING"]["enable"] is False),
        ("PE_COUPLING.apply_to_power", u2["PE_COUPLING"]["apply_to_power"] is False),
        ("SCATTER_PD.apply_to_pd", u2["SCATTER_PD"]["apply_to_pd"] is False),
        ("ABSORPTION_COMPLEX.enable", u2["ABSORPTION_COMPLEX"]["enable"] is False),
        ("RADAR_CLUTTER.enable", u2["RADAR_CLUTTER"]["enable"] is False),
        ("RADAR_CLUTTER.apply_to_power", u2["RADAR_CLUTTER"]["apply_to_power"] is False),
    ]
    n_fail = 0
    for name, ok in checks:
        print(" default {} {}".format(name, "OK" if ok else "FAIL"))
        n_fail += 0 if ok else 1
    return n_fail


def _check_provider_equivalence() -> int:
    x = np.linspace(0.0, 20.0, 5)
    z = np.linspace(80.0, 120.0, 5)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    Ne = 1.0e10 + 2.0e8 * xx + 1.0e7 * zz
    provider = GridNeProvider(Ne, x, z)

    n_grid = RefractiveIndex3D(Ne, x, z, cfg.FREQ_MHZ, earth_flat=False)
    n_provider = RefractiveIndex3D(
        np.zeros_like(Ne),
        x,
        z,
        cfg.FREQ_MHZ,
        earth_flat=False,
        ne_provider=provider,
    )

    samples = [(0.0, 0.0, 80.0), (8.0, 3.0, 96.0), (20.0, -2.0, 120.0)]
    n_fail = 0
    for pos in samples:
        ne_ref = n_grid.Ne_at(*pos)
        ne_new = n_provider.Ne_at(*pos)
        ne_ok = abs(ne_ref - ne_new) <= 1e-6 * max(ne_ref, 1.0)
        n2_ok = abs(n_grid.n2(*pos) - n_provider.n2(*pos)) <= 1e-12
        ok = ne_ok and n2_ok
        print(" provider_sample {} {}".format(pos, "OK" if ok else "FAIL"))
        n_fail += 0 if ok else 1
    return n_fail


def _check_traceability_script() -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    proc = subprocess.run(
        [sys.executable, "analysis/build_traceability_matrix.py"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    ok = proc.returncode == 0
    for rel in ("word/model_traceability_update2.md", "word/update2_validation.txt"):
        ok = ok and os.path.exists(os.path.join(root, rel))
    print(" traceability_script {}".format("OK" if ok else "FAIL"))
    if not ok:
        print(proc.stdout)
        print(proc.stderr)
    return 0 if ok else 1


def _check_local_frame_script() -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    proc = subprocess.run(
        [sys.executable, "tests/_verify_local_frame.py"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    ok = proc.returncode == 0
    print(" local_frame_script {}".format("OK" if ok else "FAIL"))
    if not ok:
        print(proc.stdout)
        print(proc.stderr)
    return 0 if ok else 1


def run() -> int:
    print("=" * 55)
    print(" update2.1 phase-1 verification")
    print("=" * 55)
    n_fail = 0
    n_fail += _check_update2_defaults()
    n_fail += _check_provider_equivalence()
    n_fail += _check_traceability_script()
    n_fail += _check_local_frame_script()
    print("update2_all: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())



