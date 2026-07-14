"""Verify that every main2 scene changes the propagated Ne field physically."""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import config as cfg
import main2


def run() -> int:
    providers = {s: main2._scenario_provider(s) for s in main2.SCENARIO_LABELS}
    base = providers["baseline"].get_Ne_grid(cfg.BG_X, cfg.BG_Z)
    mask = base > 0.01 * float(base.max())
    failures = []

    for scene in ("tid", "spreadf", "bubble", "es"):
        field = providers[scene].get_Ne_grid(cfg.BG_X, cfg.BG_Z)
        rel_max = float(np.max(np.abs(field[mask] - base[mask]) / base[mask]))
        if scene == "tid" and rel_max < 0.05:
            failures.append("TID F-region perturbation missing")
        elif scene == "spreadf" and rel_max < 0.01:
            failures.append("Spread-F density RMS missing")
        elif scene == "bubble" and rel_max < 0.20:
            failures.append("Bubble depletion missing")
        elif scene == "es" and not np.any(field > base):
            failures.append("Es enhancement missing")

    ne_es, z_es = main2._plot_ne_slice(providers["es"], "es")
    iz = int(np.argmin(np.abs(z_es - float(cfg.ES["h_Es_km"]))))
    if not np.any(ne_es[:, iz] > 0.0):
        failures.append("Es visualization grid misses thin layer")

    print("main2_scenario_fields: {} fail(s)".format(len(failures)))
    for failure in failures:
        print(" FAIL: {}".format(failure))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(run())
