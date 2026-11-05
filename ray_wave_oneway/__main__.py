"""Standalone controlled one-way reflection/scattering observables."""

from __future__ import annotations

import json
import sys
from dataclasses import asdict

import numpy as np

from .cases import run_controlled_case
from .config import ControlledCase, ONEWAY_AH_FD_SETTINGS, ONEWAY_L2_SETTINGS
from .observables import doppler_hz


def main():
    if len(sys.argv) != 1:
        print("This controlled entry uses parameters in ray_wave_oneway.config", file=sys.stderr)
        return 2
    case = ControlledCase()
    output = case.output_dir
    output.mkdir(parents=True, exist_ok=True)
    data = run_controlled_case(case)
    np.savez_compressed(output / "observables.npz", **data)
    center = case.sweep_count // 2
    slow_center = case.slow_time_count // 2
    frequency_spacing = float(data["frequency_grid_hz"][1] - data["frequency_grid_hz"][0])
    delay_bin_s = 1 / (case.sweep_count * frequency_spacing)
    doppler_bin_hz = 1 / (case.slow_time_count * case.slow_time_step_s)
    side_u_doppler = doppler_hz(data["slow_delta_u"], data["slow_time_s"])
    side_v_doppler = doppler_hz(data["slow_delta_v"], data["slow_time_s"])
    expected_u = data["slow_delta_u"][slow_center] * np.exp(-2j * np.pi * data["slow_time_s"] / 8)
    expected_v = data["slow_delta_v"][slow_center] * np.exp(-2j * np.pi * data["slow_time_s"] / 10)
    translation_u = float(np.max(np.abs(data["slow_delta_u"] - expected_u)) / abs(expected_u[slow_center]))
    translation_v = float(np.max(np.abs(data["slow_delta_v"] - expected_v)) / abs(expected_v[slow_center]))
    budget = (1 - data["frequency_reflected_power"] - data["frequency_transmitted_power"]
              - data["frequency_absorbed_power"])
    summary = {
        "route_id": "ray_wave_oneway_controlled_periodic_s4a",
        "status": "CONTROLLED_CASE_COMPUTED",
        "scope": "uniform magnetized exterior and deterministic drifting periodic Ne layer",
        "field_convention": "exp(+i*k.r-i*omega*t)",
        "frequency_phase_delay": "+d(arg H)/d(omega)",
        "physical_doppler": "f_out-f_in = -d(arg H)/dt/(2*pi)",
        "input_normalization": "source and receiver projected to current power ports and normalized to 1 W",
        "H_definition": "complex RX amplitude for 1 W projected input and unit-power RX mode",
        "physical_template_power": "Pin=source_projection_norm^2 W; Prx=abs(H)^2*Pin",
        "reference_frequency_hz": case.frequency_hz,
        "background_Ne_m3_constant_across_frequency": case.background_Ne_m3,
        "frequency_grid_hz": data["frequency_grid_hz"].tolist(),
        "frequency_spacing_hz": frequency_spacing,
        "frequency_span_hz": float(data["frequency_grid_hz"][-1] - data["frequency_grid_hz"][0]),
        "delay_fft_bin_s": delay_bin_s,
        "delay_aperture_resolution_s": 1 / float(data["frequency_grid_hz"][-1] - data["frequency_grid_hz"][0]),
        "delay_unambiguous_period_s": 1 / frequency_spacing,
        "delay_transform": "rectangular window, negative-exponent DFT, power sums to mean frequency power",
        "delay_broadening_status": "UNRESOLVED_AT_FINITE_BANDWIDTH",
        "group_delay_1khz_s": float(data["group_delay_1khz_s"]),
        "group_delay_500hz_s": float(data["group_delay_500hz_s"]),
        "reference_Pin_physical_template_W": float(data["frequency_Pin_physical_template_W"][center]),
        "reference_Prx_physical_source_W": float(data["frequency_Prx_physical_source_W"][center]),
        "reference_Prx_per_Pin": float(data["frequency_Prx_per_Pin"][center]),
        "slow_time_s": data["slow_time_s"].tolist(),
        "slow_time_step_s": case.slow_time_step_s,
        "doppler_fft_bin_hz": doppler_bin_hz,
        "doppler_transform": "rectangular window, positive-exponent DFT, power sums to mean slow-time power",
        "doppler_fft_resolution_status": "SUB_BIN_SHIFT_UNRESOLVED_BY_FFT",
        "receiver_doppler_phase_derivative_hz": float(data["slow_doppler_hz"][slow_center]),
        "u_sideband_doppler_phase_derivative_hz": float(side_u_doppler[slow_center]),
        "v_sideband_doppler_phase_derivative_hz": float(side_v_doppler[slow_center]),
        "u_translation_relative_error": translation_u,
        "v_translation_relative_error": translation_v,
        "max_frequency_power_budget_residual_W": float(np.max(np.abs(budget))),
    }
    manifest = {
        "route_id": summary["route_id"],
        "command": "python -m ray_wave_oneway",
        "parameters": asdict(case),
        "oneway_numerics": {
            "integration_L2": dict(ONEWAY_L2_SETTINGS),
            "ah_spatial_derivative": dict(ONEWAY_AH_FD_SETTINGS),
        },
        "periodic_density": "Ne_bg*(1+eta*(cos(2*pi*(u-vu*t)/Lu)+0.6*cos(2*pi*(v-vv*t)/Lv)))",
        "velocity_u_m_s": case.period_u_m / 8, "velocity_v_m_s": case.period_v_m / 10,
        "actual_medium_time_s": case.medium_time_s,
        "output": {"observables.npz": sorted(data), "summary.json": "physical definitions and headline observables"},
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False), encoding="utf-8")
    print(f"ray-wave results: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
