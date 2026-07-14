"""Verify main.py output semantics without running expensive scenarios."""

from __future__ import annotations

import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
import main
from models.hybrid_model import build_pd_spectrum
from viz.plot_utils import save_modes_csv


def _mode():
    return {
        "label": "1F",
        "tau_ms": 3.0,
        "group_path_km": 900.0,
        "Pr_W": 1.0,
        "Pr_dBW": 0.0,
        "delta_tau_ms": 0.1,
        "h_reflect_km": 280.0,
        "A_tube_km2": 20.0,
        "tau_spread_ms": 0.2,
    }


def run() -> int:
    n_fail = 0
    print("=" * 55)
    print(" main output semantics verification")
    print("=" * 55)

    comm = main._tag_communication_modes([_mode()])[0]
    ok = comm["tau_oneway_ms"] == comm["tau_ms"] and comm["tau_2way_ms"] is None
    print(" comm_oneway_semantics {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    radar = main._tag_point_target_radar_modes([_mode()])[0]
    ok = (
        radar["tau_oneway_ms"] == radar["tau_ms"]
        and abs(radar["tau_2way_ms"] - 2.0 * radar["tau_ms"]) < 1e-12
        and radar["radar_observable"] == "two_way_point_target"
    )
    print(" point_target_two_way_semantics {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    clutter = main._tag_clutter_radar_modes([_mode()])[0]
    ok = (
        abs(clutter["tau_2way_ms"] - 2.0 * clutter["tau_oneway_ms"]) < 1e-12
        and clutter["radar_observable"] == "two_way_surface_clutter"
        and "Barrick_1972" in clutter["physics_basis"]
    )
    print(" clutter_two_way_semantics {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    tau1, pd1 = build_pd_spectrum([_mode()], method="tube")
    pd_enabled = main.RUN_PRODUCTS["pd"]
    main.RUN_PRODUCTS["pd"] = True
    try:
        tau2, pd2 = main._build_pd([clutter], method="tube", delay_key="tau_2way_ms")
    finally:
        main.RUN_PRODUCTS["pd"] = pd_enabled
    ok = abs(float(tau1[np.argmax(pd1)]) * 2.0 - float(tau2[np.argmax(pd2)])) < 0.05
    print(" radar_pd_uses_two_way_delay {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    td = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
    path = save_modes_csv(
        "semantics", [clutter], tau2, pd2, cfg.FREQ_MHZ, td,
        pd_delay_key="tau_2way_ms")
    with open(path, newline="", encoding="utf-8") as f:
        row = next(csv.DictReader(f))
    ok = (
        row["pd_delay_key"] == "tau_2way_ms"
        and abs(float(row["pd_delay_ms"]) - clutter["tau_2way_ms"]) < 1e-12
        and row["radar_observable"] == "two_way_surface_clutter"
    )
    print(" csv_records_delay_semantics {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    ok = (
        cfg.UPDATE2["AH3D_ENABLE"] is False
        and cfg.UPDATE2["OX_MODE"] == "ah3d"
        and cfg.UPDATE2["RADAR_CLUTTER"]["enable"] is False
        and cfg.UPDATE2["RADAR_CLUTTER"]["apply_to_power"] is False
    )
    print(" default_gates_still_off {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    with main._temporary_update2({
        "AH3D_ENABLE": True,
        "SCATTER_PD": {"apply_to_pd": True},
    }):
        inside_ok = (
            cfg.UPDATE2["AH3D_ENABLE"] is True
            and cfg.UPDATE2["SCATTER_PD"]["apply_to_pd"] is True
        )
    restored_ok = (
        cfg.UPDATE2["AH3D_ENABLE"] is False
        and cfg.UPDATE2["SCATTER_PD"]["apply_to_pd"] is False
    )
    ok = inside_ok and restored_ok
    print(" temporary_update2_restores_defaults {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    print("main_output_semantics: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
