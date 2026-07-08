"""3-D irregularity physics diagnostics for stage-1 hybrid propagation.

This module derives scattering descriptors from parameterized irregularities.
It does not modify Ne or propagation results; downstream Born/Rytov/MPS/PE
couplers decide how diagnostics are used behind explicit config switches.
"""
from __future__ import annotations

import numpy as np

from config import ES, ES_BG, GEOMAG, SCATTER, SPREAD_F


def geomag_unit_vector(geomag: dict | None = None) -> np.ndarray:
    g = GEOMAG if geomag is None else {**GEOMAG, **geomag}
    dip = np.radians(float(g.get("dip_deg", 48.7)))
    decl = np.radians(float(g.get("decl_deg", -5.5)))
    link = np.radians(float(g.get("link_bearing_deg", 0.0)))
    b = np.array([
        np.cos(dip) * np.cos(decl - link),
        np.cos(dip) * np.sin(decl - link),
        -np.sin(dip),
    ], dtype=float)
    return b / max(float(np.linalg.norm(b)), 1e-30)


def bragg_mismatch(freq_MHz: float,
                   k_hat,
                   irregularity_axis,
                   scale_km: float) -> dict:
    """Return Bragg/aspect diagnostics for a field-aligned irregularity scale."""
    lam_km = 0.299792458 / max(float(freq_MHz), 1e-12)
    k0 = 2.0 * np.pi / lam_km
    kh = np.asarray(k_hat, dtype=float)
    kh = kh / max(float(np.linalg.norm(kh)), 1e-30)
    ah = np.asarray(irregularity_axis, dtype=float)
    ah = ah / max(float(np.linalg.norm(ah)), 1e-30)
    aspect_cos = abs(float(np.dot(kh, ah)))
    q_bragg = 2.0 * k0
    q_scale = 2.0 * np.pi / max(float(scale_km), 1e-12)
    return {
        "bragg_q_km_inv": float(q_bragg),
        "irregularity_q_km_inv": float(q_scale),
        "bragg_mismatch": float(abs(q_scale - q_bragg) / max(q_bragg, 1e-30)),
        "aspect_cos": aspect_cos,
        "aspect_angle_deg": float(np.degrees(np.arccos(np.clip(aspect_cos, -1.0, 1.0)))),
    }


def tid_descriptor(tid_params: dict, freq_MHz: float, k_hat=(1.0, 0.0, 0.0)) -> dict:
    n_comp = int(tid_params.get("n_components", 1))
    amps = tid_params.get("amplitude_list", [tid_params.get("amplitude", 0.0)])
    lams = tid_params.get("lambda_h_km_list", [tid_params.get("lambda_h_km", 300.0)])
    azs = tid_params.get("az_deg_list", [0.0])
    comps = []
    for i in range(n_comp):
        lam = float(lams[i % len(lams)])
        az = np.radians(float(azs[i % len(azs)]))
        axis = np.array([np.cos(az), np.sin(az), 0.0], dtype=float)
        b = bragg_mismatch(freq_MHz, k_hat, axis, lam)
        comps.append({
            "amplitude": float(amps[i % len(amps)]),
            "lambda_h_km": lam,
            "az_deg": float(azs[i % len(azs)]),
            **b,
        })
    return {
        "type": "TID",
        "enabled": bool(tid_params.get("enable", False)),
        "components": comps,
        "scatter_DeltaL_km": float(tid_params.get("sigma_z_km", 90.0)) * 2.0,
    }


def bubble_descriptor(bubble_params: dict,
                      freq_MHz: float,
                      k_hat=(1.0, 0.0, 0.0),
                      geomag: dict | None = None) -> dict:
    b_hat = geomag_unit_vector(geomag)
    sigma_perp = float(bubble_params.get("sigma_perp_km", bubble_params.get("Lx_km", 50.0)))
    sigma_para = float(bubble_params.get("sigma_para_km", bubble_params.get("Lz_km", 150.0)))
    b = bragg_mismatch(freq_MHz, k_hat, b_hat, sigma_perp)
    return {
        "type": "field_aligned_bubble",
        "enabled": bool(bubble_params.get("enable", False)),
        "delta_max": float(bubble_params.get("delta_max", 0.0)),
        "sigma_perp_km": sigma_perp,
        "sigma_para_km": sigma_para,
        "anisotropy_ratio": float(sigma_para / max(sigma_perp, 1e-12)),
        "b_hat": b_hat.tolist(),
        "scatter_DeltaL_km": float(2.0 * sigma_perp),
        **b,
    }


def spreadf_descriptor(spreadf_params: dict | None = None,
                       scatter_params: dict | None = None) -> dict:
    sf = SPREAD_F if spreadf_params is None else {**SPREAD_F, **spreadf_params}
    sc = SCATTER if scatter_params is None else {**SCATTER, **scatter_params}
    return {
        "type": "Spread-F",
        "enabled": bool(sf.get("enable", False)),
        "Cs": float(sf.get("Cs", sc.get("Cs_rel", 0.0))),
        "spectral_index": float(sf.get("p", sc.get("p", 3.0))),
        "outer_scale_km": float(sf.get("L0_km", sc.get("L_outer_km", 10.0))),
        "inner_scale_m": float(sc.get("L_inner_m", 0.1)),
        "scatter_DeltaL_km": float(sc.get("DeltaL_km", 100.0)),
    }


def es_descriptor(es_params: dict | None = None,
                  es_bg_params: dict | None = None) -> dict:
    es = ES if es_params is None else {**ES, **es_params}
    bg = ES_BG if es_bg_params is None else {**ES_BG, **es_bg_params}
    l_perp_km = float(bg.get("L_perp_m", es.get("L1_m", 300.0))) / 1000.0
    l_para_km = float(bg.get("L_para_m", 10000.0)) / 1000.0
    return {
        "type": "Es",
        "enabled": bool(es.get("enable", False)),
        "foEs_MHz": float(es.get("foEs_MHz", 0.0)),
        "h_Es_km": float(es.get("h_Es_km", 110.0)),
        "delta_h_km": float(es.get("delta_h_m", 0.0)) / 1000.0,
        "L_perp_km": l_perp_km,
        "L_para_km": l_para_km,
        "anisotropy_ratio": float(l_para_km / max(l_perp_km, 1e-12)),
        "spectral_index": float(bg.get("p_es", 1.5)),
        "scatter_DeltaL_km": float(2.0 * es.get("delta_h_m", 0.0) / 1000.0),
    }


def summarize_irregularities(params: dict,
                             freq_MHz: float,
                             k_hat=(1.0, 0.0, 0.0)) -> dict:
    return {
        "TID": tid_descriptor(params.get("TID", params.get("tid", {})), freq_MHz, k_hat),
        "bubble": bubble_descriptor(params.get("BUBBLE_3D", params.get("bubble", {})), freq_MHz, k_hat),
        "Spread-F": spreadf_descriptor(params.get("SPREAD_F")),
        "Es": es_descriptor(params.get("ES"), params.get("ES_BG")),
    }
