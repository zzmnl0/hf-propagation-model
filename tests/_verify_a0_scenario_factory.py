"""Fast A0 checks for reproducible main/main2 ionosphere construction."""

from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, ROOT)

import config as cfg
import main
import main2
from models.scenario_provider_factory import ScenarioProviderFactory
from utils import spherical_destination


def check(condition, label, failures):
    if not condition:
        failures.append(label)
        print("FAIL: {}".format(label))


def run():
    failures = []
    points = np.array([
        [0.0, 0.0, 100.0],
        [0.25 * cfg.RX_RANGE, 25.0, 220.0],
        [0.50 * cfg.RX_RANGE, -40.0, 300.0],
        [cfg.RX_RANGE, 0.0, 400.0],
    ])
    for scene in main2.SCENARIO_LABELS:
        a = main._scenario_provider(scene)
        b = main2._scenario_provider(scene)
        ma, mb = a.get_metadata(), b.get_metadata()
        check(ma == mb, "{}_metadata_exact".format(scene), failures)
        check(ma["config_hash"] == mb["config_hash"],
              "{}_config_hash".format(scene), failures)
        check(ma["provider_hash"] == mb["provider_hash"],
              "{}_provider_hash".format(scene), failures)
        ga = a.sample_link_plane(
            cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, cfg.BG_Z)
        gb = b.sample_link_plane(
            cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, cfg.BG_Z)
        check(np.array_equal(ga, gb), "{}_slice_exact".format(scene), failures)
        geocentric = []
        for x, y, height in points:
            bearing = cfg.LINK_BEARING_DEG + np.degrees(np.arctan2(y, x))
            lat, lon = spherical_destination(
                cfg.TX_LAT, cfg.TX_LON, bearing, np.hypot(x, y), cfg.RE_KM)
            geocentric.append((lat, lon, height))
        va = np.array([a.get_Ne_geocentric(*p, t=0.0) for p in geocentric])
        vb = np.array([b.get_Ne_geocentric(*p, t=0.0) for p in geocentric])
        check(np.allclose(va, vb, rtol=1e-12, atol=0.0),
              "{}_samples".format(scene), failures)

    factory = ScenarioProviderFactory(cfg)
    explicit = factory.create("baseline", components=("tid",), seed=0, t_s=0.0)
    labelled = factory.create("tid", seed=0, t_s=0.0)
    check(explicit.get_metadata()["provider_hash"] ==
          labelled.get_metadata()["provider_hash"],
          "label_does_not_change_physical_provider_hash", failures)
    check(np.array_equal(
        explicit.sample_link_plane(
            cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, cfg.BG_Z),
        labelled.sample_link_plane(
            cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, cfg.BG_Z)),
          "label_does_not_change_provider_data", failures)
    try:
        factory.create("spreadf", seed=None, validation=True)
    except ValueError:
        pass
    else:
        failures.append("validation_rejects_seed_none")

    print("A0 scenario factory: {} fail(s)".format(len(failures)))
    return len(failures)


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
