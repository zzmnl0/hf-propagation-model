"""Unified propagation entry point for stage-1 P0.

main.py remains the single-scenario validation/legacy comparison driver.
main2.py applies one propagation pipeline to multiple ionospheric scenes and
writes all products to output2/.
"""

from __future__ import annotations

import copy
import argparse
import json
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
    "p2p_seeds": None,
    "p2p_profile": None,
    "p2p_final_profile": None,
    "p2p_max_nfev": None,
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
    from models.unified_hybrid_pipeline import reject_archived_physics
    reject_archived_physics()
    _merge_dict(cfg.UPDATE2, {
        "AH3D_ENABLE": True,
        "OX_MODE": "ah3d",
        "RAY_FAN_INCLUDE_X": bool(RUN_PRODUCTS.get("ray_fan_x", False)),
        "ABSORPTION_COMPLEX": {
            "enable": True,
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
            "enable": False,
            "apply_to_power": False,
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


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _ray_fan_is_ox(result):
    return set(result.traceability.get("ray_fan_wave_modes", [])) == {"O", "X"}


def _save_ray_fan_summary(kind, result):
    is_ox = _ray_fan_is_ox(result)
    stem = "ray_fan_ox_unified_{}" if is_ox else "ray_fan_unified_{}"
    path = os.path.join(_OUT, stem.format(kind) + ".json")
    disabled = ["P-D", "Ray-PE", "scatter_power_coupling"]
    if not is_ox:
        disabled.insert(0, "X_ray_fan")
    if not RUN_PRODUCTS.get("p2p", False):
        disabled.insert(1 if not is_ox else 0, "P2P")
    payload = {
        "schema_version": "formal_complex_scene_ray_fan_v2",
        "scenario": kind,
        "frequency_MHz": float(cfg.FREQ_MHZ),
        "oneway_or_twoway": "one_way",
        "wave_modes": result.traceability.get("ray_fan_wave_modes", []),
        "provider_hash": result.traceability.get("provider_hash", ""),
        "provider_config_hash": result.traceability.get(
            "provider_config_hash", ""),
        "provider_metadata": result.provider_metadata,
        "ray_fan_diagnostics": result.traceability.get(
            "ray_fan_diagnostics", {}),
        "termination_counts": result.traceability.get(
            "ray_fan_termination_counts", {}),
        "warnings": result.warnings,
        "disabled_products": disabled,
    }
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(_jsonable(payload), stream, ensure_ascii=False, indent=2)
    print("  Saved: output2/{}".format(os.path.basename(path)))


def _save_p2p_summary(kind, result):
    path = os.path.join(_OUT, "p2p_unified_{}.json".format(kind))
    modes = []
    fields = (
        "label", "wave_mode", "termination_status", "branch_status",
        "elevation_tx_deg", "bearing_tx_deg", "apex_height_km",
        "geom_path_km", "group_path_km", "tau_ms",
        "p2p_landing_error_km", "p2p_residual_east_north_km",
        "p2p_trace_profile", "p2p_optimization_trace_profile",
        "p2p_final_retrace_applied", "p2p_search_scope", "p2p_seed_id",
        "p2p_seed_semantics", "p2p_nfev", "p2p_njev",
    )
    for mode in result.modes:
        modes.append({key: _jsonable(mode.get(key)) for key in fields})
    payload = {
        "schema_version": "formal_complex_scene_seeded_p2p_v1",
        "scenario": kind,
        "frequency_MHz": float(cfg.FREQ_MHZ),
        "oneway_or_twoway": "one_way",
        "search_scope": "seeded_local_roots_not_full_domain_mode_audit",
        "candidate_nearest_per_branch": int(
            RUN_PRODUCTS.get("p2p_seeds")
            or cfg.P2P_3D.get("candidate_nearest", 3)),
        "solve_profile": str(
            RUN_PRODUCTS.get("p2p_profile") or "reference"),
        "final_retrace_profile": str(
            RUN_PRODUCTS.get("p2p_final_profile")
            or RUN_PRODUCTS.get("p2p_profile") or "reference"),
        "max_nfev": int(
            RUN_PRODUCTS.get("p2p_max_nfev")
            or cfg.P2P_3D.get("max_iter", 20)),
        "completeness_claim": False,
        "provider_hash": result.traceability.get("provider_hash", ""),
        "provider_config_hash": result.traceability.get(
            "provider_config_hash", ""),
        "wave_modes": result.traceability.get("p2p_wave_modes", []),
        "root_counts": {
            branch: sum(str(mode.get("wave_mode")) == branch
                        for mode in result.modes)
            for branch in ("O", "X")
        },
        "modes": modes,
        "warnings": result.warnings,
        "disabled_products": ["P-D", "Ray-PE", "scatter_power_coupling"],
    }
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(_jsonable(payload), stream, ensure_ascii=False, indent=2)
    print("  Saved: output2/{}".format(os.path.basename(path)))


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
    from models.unified_hybrid_pipeline import reject_archived_physics
    reject_archived_physics()
    kind = _normalize_scenario_label(kind)
    apply_scatter_pd = bool(cfg.UPDATE2.get("SCATTER_PD", {}).get("apply_to_pd", False))
    p2p_config = dict(cfg.P2P_3D)
    if RUN_PRODUCTS.get("p2p_seeds") is not None:
        p2p_config["candidate_nearest"] = int(RUN_PRODUCTS["p2p_seeds"])
    if RUN_PRODUCTS.get("p2p_profile") is not None:
        p2p_config["solve_profile"] = str(RUN_PRODUCTS["p2p_profile"])
    if RUN_PRODUCTS.get("p2p_final_profile") is not None:
        p2p_config["final_retrace_profile"] = str(
            RUN_PRODUCTS["p2p_final_profile"])
    if RUN_PRODUCTS.get("p2p_max_nfev") is not None:
        p2p_config["max_iter"] = int(RUN_PRODUCTS["p2p_max_nfev"])
    return UnifiedPropagationRequest(
        scenario="unified_{}".format(kind),
        provider=_scenario_provider(kind),
        tx_geocentric=cfg.TX_GEOCENTRIC,
        rx_geocentric=cfg.RX_GEOCENTRIC,
        bearing_deg=cfg.LINK_BEARING_DEG,
        earth_radius_km=cfg.EARTH_MODEL["radius_km"],
        p2p_config=p2p_config,
        solver="ah3d_ecef",
        observable_type="communication_oneway",
        pd_delay_key="tau_ms",
        apply_scatter_to_pd=apply_scatter_pd,
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
    if RUN_PRODUCTS.get("ray_fan", True):
        _save_ray_fan_summary(kind, result)
    if RUN_PRODUCTS.get("p2p", False):
        _save_p2p_summary(kind, result)

    pert_info = _make_pert_info(kind)
    if RUN_PRODUCTS.get("ray_fan", True) and result.rays3d:
        fan_modes = "/".join(result.traceability.get("ray_fan_wave_modes", ["O"]))
        fig, _ = plot_tx_ray_fan(
            result.rays3d, Ne_plot, cfg.BG_X, z_plot,
            title="TX Ray Fan ({}) - unified {} Ne background".format(
                fan_modes, kind),
            link_info=_make_link_info(result.traceability), pert_info=pert_info,
            xlim=_plot_x_limits())
        fan_stem = ("tx_fan_ox_unified_{}" if _ray_fan_is_ox(result)
                    else "tx_fan_unified_{}")
        _save(fig, fan_stem.format(kind) + ".png")

    if RUN_PRODUCTS.get("p2p", True):
        fig, _ = plot_p2p_modes(
            result.modes, Ne_plot, cfg.BG_X, z_plot,
            title=("Seeded receiver-hitting P2P roots (O/X) - "
                   "unified {} Ne background").format(kind),
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

    if RUN_PRODUCTS.get("csv", True) and RUN_PRODUCTS.get("p2p", False):
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
    parser.add_argument(
        "--p2p", action="store_true",
        help="Run seeded local O/X P2P roots; not a full-domain mode audit")
    parser.add_argument(
        "--p2p-seeds", type=int, choices=range(1, 4), default=None,
        help="Nearest ground-ray seeds per O/X branch; default uses config")
    parser.add_argument(
        "--p2p-profile", choices=("smoke", "research", "reference"),
        default=None,
        help="P2P optimization integration profile; default is reference")
    parser.add_argument(
        "--p2p-final-profile", choices=("research", "reference"),
        default=None,
        help="Forward-retrace converged roots at this stricter profile")
    parser.add_argument(
        "--p2p-max-nfev", type=int, choices=range(1, 21), default=None,
        help="Maximum bounded-TRF residual evaluations per local seed")
    parser.add_argument(
        "--csv", action="store_true",
        help="Save the P2P mode table, including a header when no root converges")
    parser.add_argument(
        "--no-ray-fan", action="store_true",
        help="Disable ray-fan image output; P2P may still trace internal seed fans")
    args = parser.parse_args()
    RUN_PRODUCTS["ray_fan_x"] = bool(args.ray_fan_x)
    RUN_PRODUCTS["p2p"] = bool(args.p2p)
    RUN_PRODUCTS["p2p_seeds"] = args.p2p_seeds
    RUN_PRODUCTS["p2p_profile"] = args.p2p_profile
    RUN_PRODUCTS["p2p_final_profile"] = args.p2p_final_profile
    RUN_PRODUCTS["p2p_max_nfev"] = args.p2p_max_nfev
    RUN_PRODUCTS["csv"] = bool(args.csv)
    RUN_PRODUCTS["ray_fan"] = not bool(args.no_ray_fan)
    main(args.scenarios)
