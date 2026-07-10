"""Unified propagation entry point for stage-1 P0.

main.py remains the single-scenario validation/legacy comparison driver.
main2.py applies one propagation pipeline to multiple ionospheric scenes and
writes all products to output2/.
"""

from __future__ import annotations

import copy
import os
import sys
from contextlib import contextmanager

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import config as cfg
from models.ionosphere_model import IonosphereModel
from models.ne_provider import GridNeProvider, ParameterizedIonosphere3DProvider
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


def _rays3d_to_fan(rays3d):
    out = []
    for r in rays3d:
        rr = dict(r)
        rr["trajectory"] = [(float(s[0]), float(s[2]))
                             for s in r.get("trajectory_3d", [])]
        out.append(rr)
    return out


def _make_link_info():
    return {
        "TX_LAT": cfg.TX_LAT,
        "TX_LON": cfg.TX_LON,
        "freq_MHz": cfg.FREQ_MHZ,
        "Pt_W": cfg.PT_W,
        "RX_km": cfg.RX_RANGE,
        "IRI_DT": cfg.IRI_DT,
        "IRI_LAT": cfg.IRI_LAT,
        "IRI_LON": cfg.IRI_LON,
    }


def _make_pert_info(kind):
    info = {}
    if kind in ("tid", "full"):
        info["TID"] = cfg.TID_3D
    if kind in ("es", "full"):
        info["Es"] = cfg.ES
    if kind in ("bubble", "full"):
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
        "provider={}".format(result.traceability.get("provider_components", "")),
        "solver={}".format(result.traceability.get("solver_chain", "")),
        "absorption={}".format(sample.get("abs_model", "none")),
        "Ray-PE={}".format(pe_state),
        "scatter={}".format(scatter_state),
        "pd_kernel={}".format(sample.get("pd_kernel", "")),
        "calibration={}".format(sample.get("calibration_status", "")),
    ]
    return "\n".join(lines)


def _build_base_provider():
    iono = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    return GridNeProvider(Ne_2d, cfg.BG_X, cfg.BG_Z)


def _scenario_provider(kind):
    base = _build_base_provider()
    if kind == "baseline":
        return ParameterizedIonosphere3DProvider(base, geomag=cfg.GEOMAG)
    if kind == "tid":
        return ParameterizedIonosphere3DProvider(
            base, tid_params={**cfg.TID_3D, "enable": True},
            geomag=cfg.GEOMAG)
    if kind == "bubble":
        return ParameterizedIonosphere3DProvider(
            base, bubble_params={**cfg.BUBBLE_3D, "enable": True},
            geomag=cfg.GEOMAG)
    if kind == "spreadf":
        iono = IonosphereModel(spread_f_params={**cfg.SPREAD_F, "enable": True})
        Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
        return ParameterizedIonosphere3DProvider(
            GridNeProvider(Ne_2d, cfg.BG_X, cfg.BG_Z),
            spreadf_params={**cfg.SPREAD_F, "enable": True},
            geomag=cfg.GEOMAG)
    if kind == "es":
        iono = IonosphereModel(es_params={**cfg.ES, "enable": True})
        Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
        return ParameterizedIonosphere3DProvider(
            GridNeProvider(Ne_2d, cfg.BG_X, cfg.BG_Z),
            es_params={**cfg.ES, "enable": True},
            geomag=cfg.GEOMAG)
    if kind == "full":
        iono = IonosphereModel(es_params={**cfg.ES, "enable": True})
        Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
        return ParameterizedIonosphere3DProvider(
            GridNeProvider(Ne_2d, cfg.BG_X, cfg.BG_Z),
            tid_params={**cfg.TID_3D, "enable": True},
            bubble_params={**cfg.BUBBLE_3D, "enable": True},
            es_params={**cfg.ES, "enable": True},
            geomag=cfg.GEOMAG)
    raise ValueError("unknown scenario: {}".format(kind))


def _request_for(kind):
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
    elif kind in ("bubble", "full"):
        pe_params = {**cfg.BUBBLE_3D, "enable": True}
    return UnifiedPropagationRequest(
        scenario="unified_{}".format(kind),
        provider=_scenario_provider(kind),
        solver="ah3d" if bool(cfg.UPDATE2.get("AH3D_ENABLE", False)) else "iso3d_legacy",
        observable_type="communication_oneway",
        pd_delay_key="tau_ms",
        apply_scatter_to_pd=apply_scatter_pd,
        apply_ray_pe_reference=apply_pe,
        pe_irregularity_params=pe_params,
        environment=kind,
        metadata={"main2_output_dir": _OUT},
    )


def run_unified_scenario(kind):
    print("=" * 55)
    print("  main2 unified scenario: {}".format(kind))
    print("=" * 55)
    request = _request_for(kind)
    result = run_unified_pipeline(request)
    print("  P2P modes: {}".format(len(result.modes)))
    print_mode_table(result.modes)
    if result.main_mode:
        print("  Main mode: {} tau={:.3f} ms Pr={:.1f} dBW".format(
            result.main_mode.get("label", ""),
            float(result.main_mode.get("tau_ms", 0.0)),
            float(result.main_mode.get("Pr_dBW", -999.0))))
    for warning in result.warnings:
        print("  Warning: {}".format(warning))

    pert_info = _make_pert_info(kind)
    if result.rays3d:
        fig, _ = plot_tx_ray_fan(
            _rays3d_to_fan(result.rays3d), result.Ne_2d, cfg.BG_X, cfg.BG_Z,
            title="TX Ray Fan - unified {}".format(kind),
            link_info=_make_link_info(), pert_info=pert_info)
        _save(fig, "tx_fan_unified_{}.png".format(kind))

    fig, _ = plot_p2p_modes(
        result.modes, result.Ne_2d, cfg.BG_X, cfg.BG_Z,
        title="P2P Modes - unified {}".format(kind),
        link_info=_make_link_info(), pert_info=pert_info)
    _save(fig, "p2p_modes_unified_{}.png".format(kind))

    info_text = _make_pd_info(result, delay_key=request.pd_delay_key)
    fig, _ = plot_pd_spectrum(
        result.tau_axis, result.pd_W, result.modes,
        title="P-D - unified {} (one-way)".format(kind),
        delay_key=request.pd_delay_key,
        xlabel="One-way Group Delay  (ms)",
        info_text=info_text)
    _save(fig, "pd_unified_{}.png".format(kind))

    csv_path = save_modes_csv(
        "unified_{}".format(kind), result.modes, result.tau_axis, result.pd_W,
        cfg.FREQ_MHZ, _OUT, pd_delay_key=request.pd_delay_key)
    print("  Saved: output2/{}".format(os.path.basename(csv_path)))
    return result


def main():
    _enable_reference_physics()
    scenarios = ["baseline", "tid", "spreadf", "bubble", "es", "full"]
    return [run_unified_scenario(s) for s in scenarios]


if __name__ == "__main__":
    main()
