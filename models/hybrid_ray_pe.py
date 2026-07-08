"""
2.5-D Ray-PE/SSF coupling utilities.

This module keeps the first update2.1 Ray-PE interface local and explicit:
it uses a 2-D PE cross-section around the ray/bubble entry, reports an
energy-based power correction, and only applies that correction when the
caller enables the config switch.
"""

from __future__ import annotations

import numpy as np

from config import BUBBLE, PE
from .pe_propagator import PEPropagator, construct_incident_field
from .local_frame import extract_irregularity_entry, make_boundary_field_descriptor


def _power_ratio_to_db(ratio: float) -> float:
    return float(10.0 * np.log10(max(float(ratio), 1e-300)))


def apply_ray_pe_coupling(
    pe: PEPropagator,
    Ne_2d: np.ndarray,
    x_array: np.ndarray,
    z_array: np.ndarray,
    Pr_W: float,
    k0_per_km: float,
    x_entry_km: float,
    z_entry_km: float,
    beta_inc_deg: float,
    bubble_params: dict,
    coupling_params: dict | None = None,
) -> tuple[float, dict]:
    """Run a local 2.5-D PE/SSF section and return optional power update.

    Parameters are in the same units used by the rest of the project:
    distances in km, electron density in m^-3, and power in W.

    Returns
    -------
    Pr_out_W:
        Input power, or input power multiplied by the SSF energy ratio when
        ``coupling_params['apply_to_power']`` is true.
    diagnostics:
        PE/SSF diagnostic fields suitable for attaching to a mode dict.
    """

    params = coupling_params or {}
    apply_to_power = bool(params.get("apply_to_power", False))

    z0 = float(bubble_params.get("z0_km", BUBBLE["z0_km"]))
    Lz = float(bubble_params.get("Lz_km", BUBBLE["Lz_km"]))
    Lx = float(bubble_params.get("Lx_km", BUBBLE["Lx_km"]))
    h_bot = z0 - Lz

    dx_km = float(pe.params.get("dx_km", PE["dx_km"]))
    dz_km = float(pe.params.get("dz_m", PE["dz_m"])) / 1000.0
    x_rng = (
        max(float(x_entry_km), float(x_array[0])),
        min(float(x_entry_km) + 2.0 * Lx, float(x_array[-1])),
    )
    z_rng = (
        max(h_bot, float(z_array[0])),
        min(z0 + Lz, float(z_array[-1])),
    )

    if x_rng[1] <= x_rng[0] or z_rng[1] <= z_rng[0]:
        diag = {
            "pe_coupled": False,
            "pe_apply_to_power": apply_to_power,
            "pe_status": "empty_domain",
            "pe_Pr_correction_dB": 0.0,
            "pe_power_ratio": 1.0,
            "pe_delta_tau_ms": 0.0,
        }
        return float(Pr_W), diag

    n_pe, x_pe, z_pe = pe.extract_domain(Ne_2d, x_array, z_array, x_rng, z_rng)
    A_inc = float(np.sqrt(max(float(Pr_W), 0.0)))
    u_init = construct_incident_field(
        A_inc, float(beta_inc_deg), float(z_entry_km), z_pe, k0_per_km
    )
    u_out, _history = pe.propagate(u_init, n_pe, dx_km, dz_km)

    dx_total_km = max(float(x_pe[-1] - x_pe[0]), 0.0)
    pe_res = pe.analyze(u_out, z_pe, dx_total_km)

    P_in = float(np.sum(np.abs(u_init) ** 2))
    P_out = float(np.sum(np.abs(u_out) ** 2))
    ratio = P_out / P_in if P_in > 0.0 else 1.0
    ratio = max(float(ratio), 0.0)

    Pr_out = float(Pr_W) * ratio if apply_to_power else float(Pr_W)
    diag = {
        "pe_coupled": True,
        "pe_apply_to_power": apply_to_power,
        "pe_status": "ok",
        "pe_Pr_correction_dB": _power_ratio_to_db(ratio),
        "pe_power_ratio": ratio,
        "pe_P_in": P_in,
        "pe_P_out": P_out,
        "pe_delta_tau_ms": float(pe_res.get("delta_tau_ms", 0.0)),
        "pe_mean_aoa_deg": float(pe_res.get("mean_aoa_deg", 0.0)),
        "pe_tau_extra_mean_ms": float(pe_res.get("tau_extra_mean_ms", 0.0)),
        "pe_x_entry_km": float(x_entry_km),
        "pe_z_entry_km": float(z_entry_km),
        "pe_dx_total_km": dx_total_km,
        "pe_boundary_amplitude_sqrt_W": A_inc,
        "pe_boundary_phase_rad": 0.0,
        "pe_boundary_curvature_km_inv": 0.0,
        "pe_boundary_polarization": "scalar",
    }
    return Pr_out, diag



def apply_ray_pe_coupling_3d(
    pe: PEPropagator,
    mode: dict,
    Ne_2d: np.ndarray,
    x_array: np.ndarray,
    z_array: np.ndarray,
    Pr_W: float,
    k0_per_km: float,
    irregularity_params: dict,
    coupling_params: dict | None = None,
) -> tuple[float, dict]:
    """3-D mode interface for local_2p5d Ray-PE/SSF diagnostics.

    The entry point and frame come from the 3-D ray and field-aligned
    irregularity geometry. The current PE solver remains the validated 2.5-D
    section; diagnostics explicitly report local_2p5d frame status.
    """
    params = coupling_params or {}
    apply_to_power = bool(params.get('apply_to_power', False))
    entry = extract_irregularity_entry(mode, irregularity_params)
    if entry is None:
        return float(Pr_W), {
            'pe_coupled': False,
            'pe_apply_to_power': apply_to_power,
            'pe_status': 'no_entry',
            'pe_frame_status': 'no_entry',
            'pe_Pr_correction_dB': 0.0,
            'pe_power_ratio': 1.0,
            'pe_delta_tau_ms': 0.0,
        }

    bubble_like = {
        'z0_km': float(irregularity_params.get('z0_km', entry['z_km'])),
        'Lz_km': float(irregularity_params.get('Lz_km', irregularity_params.get('sigma_perp_km', 30.0))),
        'Lx_km': float(irregularity_params.get('Lx_km', irregularity_params.get('sigma_perp_km', 50.0))),
    }
    mode_for_field = dict(mode)
    mode_for_field['Pr_W'] = Pr_W
    mode_for_field['k0_per_km'] = k0_per_km
    field = make_boundary_field_descriptor(
        mode_for_field, frame=entry['frame'], sample_index=entry['index'],
        Pr_W=Pr_W, k0_per_km=k0_per_km)

    Pr_out, diag = apply_ray_pe_coupling(
        pe, Ne_2d, x_array, z_array, Pr_W, k0_per_km,
        entry['x_km'], entry['z_km'], field.get('beta_inc_deg', entry['beta_inc_deg']),
        bubble_like, params)
    frame = entry['frame']
    diag.update({
        'pe_x_entry_km': float(entry['x_km']),
        'pe_y_entry_km': float(entry['y_km']),
        'pe_z_entry_km': float(entry['z_km']),
        'pe_entry_q': float(entry['q']),
        'pe_frame_status': 'local_2p5d',
        'pe_frame_sx': float(frame['s_hat'][0]),
        'pe_frame_sy': float(frame['s_hat'][1]),
        'pe_frame_sz': float(frame['s_hat'][2]),
        'pe_boundary_amplitude_sqrt_W': float(field.get('amplitude_sqrt_W', 0.0)),
        'pe_boundary_phase_rad': float(field.get('phase_rad', 0.0)),
        'pe_boundary_curvature_km_inv': float(field.get('curvature_km_inv', 0.0)),
        'pe_boundary_kx': float(np.asarray(field.get('k_hat', [1.0, 0.0, 0.0]), dtype=float)[0]),
        'pe_boundary_ky': float(np.asarray(field.get('k_hat', [1.0, 0.0, 0.0]), dtype=float)[1]),
        'pe_boundary_kz': float(np.asarray(field.get('k_hat', [1.0, 0.0, 0.0]), dtype=float)[2]),
        'pe_boundary_polarization': str(field.get('polarization', mode.get('wave_mode', 'scalar'))),
    })
    if 'pe_angle_spread_deg' not in diag:
        diag['pe_angle_spread_deg'] = float(diag.get('pe_mean_aoa_deg', 0.0))
    return Pr_out, diag
