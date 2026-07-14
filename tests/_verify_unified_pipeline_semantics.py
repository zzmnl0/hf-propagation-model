"""Verify P0 unified pipeline semantics without expensive scenario runs."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import config as cfg
import main2
from models import unified_hybrid_pipeline as up
from models.unified_hybrid_pipeline import UnifiedPropagationRequest, run_unified_pipeline
from viz.plot_utils import save_modes_csv


class FakeProvider:
    name = "fake_provider"

    def get_Ne_grid(self, x_km, z_km, t=0.0, y_km=0.0):
        return np.zeros((len(x_km), len(z_km)), dtype=float)

    def get_metadata(self):
        return {
            "name": self.name,
            "dimensions": "3d",
            "units": "m^-3",
            "tid_enabled": True,
            "bubble_enabled": False,
            "spreadf_diagnostic_only": False,
            "es_diagnostic_only": False,
        }


class DummyN3D:
    def Ne_at(self, x, y, z):
        return 0.0


def _mode():
    return {
        "label": "1F",
        "tau_ms": 3.0,
        "group_path_km": 900.0,
        "h_reflect_km": 280.0,
        "trajectory_3d": [(0.0, 0.0, 0.0), (450.0, 0.0, 280.0)],
    }


def _run_fake_pipeline(observable_type="communication_oneway", delay_key="tau_ms"):
    old_n3d = up.RefractiveIndex3D
    old_fan = up.shoot_rays_fan_3d
    old_p2p = up.find_all_rays_3d
    old_scatter = up._add_scatter_diagnostics
    old_annotate = up.annotate_mode_causes_3d
    try:
        up.RefractiveIndex3D = lambda *args, **kwargs: DummyN3D()
        up.shoot_rays_fan_3d = lambda *args, **kwargs: []
        up.find_all_rays_3d = lambda *args, **kwargs: [_mode()]
        up._add_scatter_diagnostics = lambda *args, **kwargs: None
        up.annotate_mode_causes_3d = lambda *args, **kwargs: None
        req = UnifiedPropagationRequest(
            scenario="verify",
            provider=FakeProvider(),
            x_km=np.array([0.0, 1.0]),
            z_km=np.array([0.0, 1.0]),
            make_ray_fan=False,
            observable_type=observable_type,
            pd_delay_key=delay_key,
        )
        return run_unified_pipeline(req)
    finally:
        up.RefractiveIndex3D = old_n3d
        up.shoot_rays_fan_3d = old_fan
        up.find_all_rays_3d = old_p2p
        up._add_scatter_diagnostics = old_scatter
        up.annotate_mode_causes_3d = old_annotate


def run() -> int:
    n_fail = 0
    print("=" * 55)
    print(" unified pipeline semantics verification")
    print("=" * 55)

    ok = os.path.basename(main2._OUT) == "output2"
    print(" main2_output2 {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    expected_components = {
        "baseline": [],
        "tid": ["TID3D"],
        "spreadf": ["SpreadF3D"],
        "bubble": ["field_aligned_bubble"],
        "es": ["Es3D"],
    }
    scene_fields = {}
    scene_ok = True
    for scene, expected in expected_components.items():
        provider = main2._scenario_provider(scene)
        meta = provider.get_metadata()
        enabled = [
            c["component_type"] for c in meta.get("components", [])
            if c.get("enabled", False)
        ]
        scene_fields[scene] = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z)
        scene_ok = scene_ok and enabled == expected
    baseline = scene_fields["baseline"]
    scene_ok = scene_ok and all(
        np.any(np.abs(scene_fields[s] - baseline) > 0.0)
        for s in ("tid", "spreadf", "bubble", "es")
    )
    print(" scenario_label_single_provider_field {}".format(
        "OK" if scene_ok else "FAIL"))
    n_fail += 0 if scene_ok else 1

    ok = (
        cfg.UPDATE2["AH3D_ENABLE"] is False
        and cfg.UPDATE2["OX_MODE"] == "ah3d"
        and cfg.UPDATE2["SCATTER_PD"]["apply_to_pd"] is False
        and cfg.UPDATE2["PE_COUPLING"]["enable"] is False
    )
    print(" defaults_still_off_on_import {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    src_main = Path(ROOT, "main.py").read_text(encoding="utf-8-sig")
    src_plot = Path(ROOT, "viz", "plot_utils.py").read_text(encoding="utf-8-sig")
    src_pipeline = Path(ROOT, "models", "unified_hybrid_pipeline.py").read_text(encoding="utf-8-sig")
    ok = (
        "cmap='jet'" in src_plot
        and "alpha=0.5" in src_plot
        and "cmap='Blues'" not in src_plot
        and "cmap='turbo'" not in src_plot
        and "cmap='Blues'" not in src_main
    )
    print(" ne_background_jet_alpha50 {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    ok = (
        "shoot_rays_fan_ah_3d" in src_pipeline
        and "RefractiveIndexAH3D" in src_pipeline
        and 'if request.solver == "ah3d":' in src_pipeline
    )
    print(" ah3d_tx_fan_source {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    comm = _run_fake_pipeline()
    radar = _run_fake_pipeline(
        observable_type="radar_point_target_twoway",
        delay_key="tau_2way_ms")
    cm = comm.modes[0]
    rm = radar.modes[0]
    ok = (
        cm["tau_oneway_ms"] == cm["tau_ms"]
        and cm["tau_2way_ms"] is None
        and rm["tau_2way_ms"] == 2.0 * rm["tau_oneway_ms"]
        and cm["group_path_km"] == rm["group_path_km"]
        and rm["delay_key"] == "tau_2way_ms"
    )
    print(" comm_radar_delay_semantics {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    ok = (
        cm["provider_name"] == "fake_provider"
        and cm["provider_components"] == "tid_enabled"
        and cm["solver_chain"].startswith("iso3d_legacy")
        and cm["solver_is_legacy"] is True
        and cm["ox_solver"] == "iso3d_legacy"
        and cm["absorption_applied_once"] is False
        and cm["observable_type"] == "communication_oneway"
        and cm["calibration_status"] == "reference_uncalibrated"
    )
    print(" mode_traceability_fields {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    out_dir = Path(ROOT, "output")
    path = save_modes_csv(
        "unified_semantics", comm.modes, comm.tau_axis, comm.pd_W,
        cfg.FREQ_MHZ, str(out_dir), pd_delay_key="tau_ms")
    with open(path, newline="", encoding="utf-8") as f:
        row = next(csv.DictReader(f))
    ok = (
        row["provider_name"] == "fake_provider"
        and row["solver_chain"].startswith("iso3d_legacy")
        and row["observable_type"] == "communication_oneway"
        and row["oneway_or_twoway"] == "one_way"
    )
    print(" csv_unified_fields {}".format("OK" if ok else "FAIL"))
    n_fail += 0 if ok else 1

    print("unified_pipeline_semantics: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
