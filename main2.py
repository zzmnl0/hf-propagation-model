"""Unified propagation entry point for stage-1 P0.

main.py remains the single-scenario validation/legacy comparison driver.
main2.py applies one propagation pipeline to multiple ionospheric scenes and
writes all products to output2/.
"""

from __future__ import annotations

import copy
import argparse
import os
import sys
from contextlib import contextmanager

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import config as cfg
from models.scenario_provider_factory import (
    SCENARIO_COMPONENTS,
    ScenarioProviderFactory,
    normalize_scenario_label,
)
from models.unified_hybrid_pipeline import (
    UnifiedPropagationRequest,
    run_unified_pipeline,
)
from viz.plot_utils import (
    build_ray_info_text,
    plot_pd_spectrum,
    plot_p2p_modes,
    plot_tx_ray_fan,
    print_mode_table,
    save_modes_csv,
)


_OUT = os.path.join(os.path.dirname(__file__), "output2")
os.makedirs(_OUT, exist_ok=True)

RUN_PRODUCTS = {
    "ray_fan": True,
    "ray_fan_x": False,
    "p2p": False,
    "pd": False,
    "csv": False,
}

SCENARIO_LABELS = tuple(SCENARIO_COMPONENTS)
_PROVIDER_FACTORY = ScenarioProviderFactory(cfg)


def _merge_dict(dst, src):
    for key, val in src.items():
        if isinstance(val, dict) and isinstance(dst.get(key), dict):
            _merge_dict(dst[key], val)
        else:
            dst[key] = val


@contextmanager
def _temporary_update2(overrides):
    old = copy.deepcopy(cfg.UPDATE2)
    try:
        _merge_dict(cfg.UPDATE2, overrides)
        yield
    finally:
        cfg.UPDATE2.clear()
        cfg.UPDATE2.update(old)


def _enable_reference_physics():
    _merge_dict(cfg.UPDATE2, {
        "AH3D_ENABLE": True,
        "OX_MODE": "ah3d",
        "RAY_FAN_INCLUDE_X": bool(RUN_PRODUCTS.get("ray_fan_x", False)),
        "ABSORPTION_COMPLEX": {
            "enable": True,
            "mode": "reference",
        },
        "PE_COUPLING": {
            "enable": False,
            "apply_to_power": False,
            "mode": "local_2p5d",
        },
        "SCATTER_PD": {
            "apply_to_pd": False,
            "kernel": "diagnostic_only",
        },
        "RADAR_CLUTTER": {
            "enable": True,
            "apply_to_power": True,
            "surface_type": "sea",
            "sigma0_dB": cfg.RADAR.get("sigma0_ground_dB", -20.0),
            "radial_velocity_ms": 0.0,
            "clutter_class": "surface_bragg",
        },
    })


def _save(fig, name):
    fig.savefig(os.path.join(_OUT, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: output2/{}".format(name))


def _make_link_info(traceability=None):
    trace = traceability or {}
    return {
        "TX_LAT": cfg.TX_LAT,
        "TX_LON": cfg.TX_LON,
        "freq_MHz": cfg.FREQ_MHZ,
        "Pt_W": cfg.PT_W,
        "RX_km": cfg.RX_RANGE,
        "IRI_DT": cfg.IRI_DT,
        "IRI_LAT": cfg.IRI_LAT,
        "IRI_LON": cfg.IRI_LON,
        "earth_model": cfg.EARTH_MODEL["name"],
        "earth_radius_km": cfg.EARTH_MODEL["radius_km"],
        "coordinate_model": "geocentric lat/lon/height -> ECEF state",
        "oneway_or_twoway": trace.get("oneway_or_twoway", "one_way"),
        "reflection_coupling": "on",
        "ray_pe_coupling": trace.get("ray_pe_coupling", "off"),
        "scatter_coupling": trace.get("scatter_power_coupling", "diagnostic_only"),
    }


def _make_pert_info(kind):
    info = {}
    if kind == "tid":
        info["TID"] = cfg.TID_3D
    if kind == "es":
        info["Es"] = cfg.ES
    if kind == "bubble":
        info["Bubble"] = cfg.BUBBLE_3D
    if kind == "spreadf":
        info["SpreadF"] = cfg.SPREAD_F
    return info


def _make_pd_info(result, delay_key="tau_ms"):
    modes = result.modes
    sample = next((m for m in modes if isinstance(m, dict)), {})
    obs = "two-way radar observable" if delay_key == "tau_2way_ms" else "one-way communication observable"
    pe_state = "on" if any(bool(m.get("pe_apply_to_power", False)) for m in modes) else "off"
    scatter_state = sample.get("scatter_status", "diagnostic/off")
    lines = [
        obs,
        "delay_key={}".format(delay_key),
        "earth={} R={}km".format(
            result.traceability.get("earth_model", "sphere"),
            result.traceability.get("earth_radius_km", "")),
        "coordinates={}".format(result.traceability.get("coordinate_model", "")),
        "provider={}".format(result.traceability.get("provider_components", "")),
        "solver={}".format(result.traceability.get("solver_chain", "")),
        "absorption={}".format(sample.get("abs_model", "none")),
        "Ray-PE={}".format(pe_state),
        "scatter={}".format(scatter_state),
        "pd_kernel={}".format(sample.get("pd_kernel", "")),
        "calibration={}".format(sample.get("calibration_status", "")),
    ]
    return "\n".join(lines)


def _plot_x_limits():
    return (float(cfg.BG_X[0]), float(cfg.BG_X[-1]))


def _plot_ne_slice(provider, kind):
    """Sample the same propagation provider on a visualization-aware grid."""
    z_plot = np.asarray(cfg.BG_Z, dtype=float)
    if kind == "es":
        h_es = float(cfg.ES["h_Es_km"])
        sigma_km = float(cfg.ES["delta_h_m"]) / 1000.0
        z_es = np.linspace(h_es - 5.0 * sigma_km,
                           h_es + 5.0 * sigma_km, 41)
        z_plot = np.unique(np.concatenate((z_plot, z_es, [h_es])))
    return provider.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_plot), z_plot


def _normalize_scenario_label(kind):
    return normalize_scenario_label(kind)


def _scenario_provider(kind):
    kind = _normalize_scenario_label(kind)
    configured_seed = cfg.SPREAD_F.get("seed", 0)
    seed = 0 if configured_seed is None else int(configured_seed)
    return _PROVIDER_FACTORY.create(kind, seed=seed, t_s=0.0, validation=True)


def _request_for(kind):
    kind = _normalize_scenario_label(kind)
    apply_scatter_pd = bool(cfg.UPDATE2.get("SCATTER_PD", {}).get("apply_to_pd", False))
    apply_pe = bool(cfg.UPDATE2.get("PE_COUPLING", {}).get("enable", False))
    pe_params = None
    if kind == "spreadf":
        pe_params = {
            **cfg.BUBBLE_3D,
            "enable": True,
            "x0_km": cfg.RX_RANGE / 2.0,
            "z0_km": cfg.SPREAD_F["h_screen_km"],
            "sigma_perp_km": cfg.SPREAD_F.get("L0_km", 50.0),
            "sigma_para_km": 300.0,
        }
    elif kind == "bubble":
        pe_params = {**cfg.BUBBLE_3D, "enable": True}
    return UnifiedPropagationRequest(
        scenario="unified_{}".format(kind),
        provider=_scenario_provider(kind),
        tx_geocentric=cfg.TX_GEOCENTRIC,
        rx_geocentric=cfg.RX_GEOCENTRIC,
        bearing_deg=cfg.LINK_BEARING_DEG,
        earth_radius_km=cfg.EARTH_MODEL["radius_km"],
        solver="ah3d_ecef",
        observable_type="communication_oneway",
        pd_delay_key="tau_ms",
        apply_scatter_to_pd=apply_scatter_pd,
        apply_ray_pe_reference=apply_pe,
        pe_irregularity_params=pe_params,
        environment=kind,
        metadata={"main2_output_dir": _OUT},
        make_ray_fan=bool(RUN_PRODUCTS.get("ray_fan", True)),
        make_p2p=bool(RUN_PRODUCTS.get("p2p", True)),
        make_pd=bool(RUN_PRODUCTS.get("pd", True)),
    )


def run_unified_scenario(kind):
    kind = _normalize_scenario_label(kind)
    print("=" * 55)
    print("  main2 unified scenario: {}".format(kind))
    print("=" * 55)
    request = _request_for(kind)
    result = run_unified_pipeline(request)
    Ne_plot, z_plot = _plot_ne_slice(request.provider, kind)
    print("  P2P modes: {}".format(len(result.modes)))
    if RUN_PRODUCTS.get("p2p", True):
        print_mode_table(result.modes)
    if result.main_mode:
        print("  Main mode: {} tau={:.3f} ms Pr={:.1f} dBW".format(
            result.main_mode.get("label", ""),
            float(result.main_mode.get("tau_ms", 0.0)),
            float(result.main_mode.get("Pr_dBW", -999.0))))
    for warning in result.warnings:
        print("  Warning: {}".format(warning))

    pert_info = _make_pert_info(kind)
    if RUN_PRODUCTS.get("ray_fan", True) and result.rays3d:
        fan_modes = "/".join(result.traceability.get("ray_fan_wave_modes", ["O"]))
        fig, _ = plot_tx_ray_fan(
            result.rays3d, Ne_plot, cfg.BG_X, z_plot,
            title="TX Ray Fan ({}) - unified {} Ne background".format(
                fan_modes, kind),
            link_info=_make_link_info(result.traceability), pert_info=pert_info,
            xlim=_plot_x_limits())
        _save(fig, "tx_fan_unified_{}.png".format(kind))

    if RUN_PRODUCTS.get("p2p", True) and result.modes:
        fig, _ = plot_p2p_modes(
            result.modes, Ne_plot, cfg.BG_X, z_plot,
            title="P2P Modes (O/X) - unified {} Ne background".format(kind),
            link_info=_make_link_info(result.traceability), pert_info=pert_info,
            xlim=_plot_x_limits())
        _save(fig, "p2p_modes_unified_{}.png".format(kind))

    if RUN_PRODUCTS.get("pd", True) and result.tau_axis.size and result.pd_W.size:
        info_text = _make_pd_info(result, delay_key=request.pd_delay_key)
        fig, _ = plot_pd_spectrum(
            result.tau_axis, result.pd_W, result.modes,
            title="P-D - unified {} (one-way)".format(kind),
            delay_key=request.pd_delay_key,
            xlabel="One-way Group Delay  (ms)",
            info_text=info_text)
        _save(fig, "pd_unified_{}.png".format(kind))

    if RUN_PRODUCTS.get("csv", True) and result.modes:
        csv_path = save_modes_csv(
            "unified_{}".format(kind), result.modes, result.tau_axis, result.pd_W,
            cfg.FREQ_MHZ, _OUT, pd_delay_key=request.pd_delay_key)
        print("  Saved: output2/{}".format(os.path.basename(csv_path)))
    return result


def main(scenarios=None):
    _enable_reference_physics()
    scenarios = SCENARIO_LABELS if scenarios is None else scenarios
    return [run_unified_scenario(s) for s in scenarios]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified AH3D scenario runner")
    parser.add_argument(
        "--scenario", action="append", dest="scenarios",
        help="Scenario label; repeat for multiple scenes: {}".format(
            ", ".join(SCENARIO_LABELS)))
    parser.add_argument(
        "--ray-fan-x", action="store_true",
        help="Include and distinguish X-mode rays in the ray fan; default is O only")
    args = parser.parse_args()
    RUN_PRODUCTS["ray_fan_x"] = bool(args.ray_fan_x)
    main(args.scenarios)
