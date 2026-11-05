"""Check the standalone S4a output against physical and transform controls."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ray_wave_oneway.config import ControlledCase
from ray_wave_oneway.observables import doppler_hz, phase_delay_s, power_delay, power_doppler
from ray_wave_oneway.propagation import _power_inner


def compare(metrics, name, value, limit):
    number = float(value)
    metrics[name] = {"value": number, "limit": float(limit),
                     "pass": bool(np.isfinite(number) and number <= limit)}


def main():
    case = ControlledCase()
    output = case.output_dir
    metrics = {}
    try:
        if not (output / "observables.npz").exists() or not (output / "summary.json").exists():
            raise FileNotFoundError("run python -m ray_wave_oneway before this validation")
        data = np.load(output / "observables.npz", allow_pickle=False)
        summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        frequencies = data["frequency_grid_hz"]
        full_frequencies = data["all_frequencies_hz"]
        H, H0 = data["frequency_H"], data["frequency_H0"]
        compare(metrics, "frequency_Ne_fixed", np.max(np.abs(data["frequency_background_ne"]
                - case.background_Ne_m3)) / case.background_Ne_m3, 1e-14)
        compare(metrics, "frequency_count_and_grid", float(
            len(frequencies) != 9 or len(full_frequencies) != 13), 0.)
        compare(metrics, "frequency_scatter_nonzero", max(0., 1e-8 - np.max(np.abs(H - H0))), 0.)
        compare(metrics, "group_delay_step_convergence_s", abs(
            float(data["group_delay_1khz_s"]) - float(data["group_delay_500hz_s"])),
            max(1e-9, .001 * abs(float(data["group_delay_500hz_s"]))))
        compare(metrics, "source_normalization", np.max(np.abs(
            data["frequency_Pin_physical_template_W"]
            - data["frequency_source_projection_norm"]**2)), 1e-12)
        compare(metrics, "received_power_ratio", np.max(np.abs(
            data["frequency_Prx_per_Pin"] - abs(H)**2)), 1e-12)
        compare(metrics, "received_power_physical_source", np.max(np.abs(
            data["frequency_Prx_physical_source_W"]
            - data["frequency_Prx_per_Pin"] * data["frequency_Pin_physical_template_W"])), 1e-12)

        n = len(data["indices"])
        template_source_E, template_source_H = data["source_template_E"], data["source_template_H"]
        template_receiver_E, template_receiver_H = data["receiver_template_E"], data["receiver_template_H"]
        replay_error, baseline_error, template_error = 0., 0., 0.
        thickness = 2 * 299792458.0 / case.frequency_hz
        for i in range(len(full_frequencies)):
            source, receiver = data["frequency_source"][i], data["frequency_receiver"][i]
            before, after = data["frequency_phase_before"][i], data["frequency_phase_after"][i]
            at_window = (source * np.exp(1j * before)).reshape(2 * n)
            after_phase = np.exp(1j * after)
            transmitted = (data["frequency_S"][i, 2 * n:, :2 * n] @ at_window).reshape(n, 2) * after_phase
            transmitted0 = (data["frequency_S0"][i, 2 * n:, :2 * n] @ at_window).reshape(n, 2) * after_phase
            replay_error = max(replay_error, abs(np.vdot(receiver, transmitted) - H[i]),
                               abs(np.vdot(receiver, transmitted0) - H0[i]))
            k_s = data["frequency_port_k_ecef_rad_km"][i, :, :, 0] / 1000
            analytic0 = source * np.exp(1j * k_s * (case.before_m + thickness + case.after_m))
            baseline_error = max(baseline_error, abs(np.vdot(receiver, analytic0) - H0[i]))
            for side, target, electric, magnetic in (
                ("source", source, template_source_E, template_source_H),
                ("receiver", receiver, template_receiver_E, template_receiver_H),
            ):
                port_E = data["frequency_source_port_electric" if side == "source" else "frequency_port_electric"][i]
                port_H = data["frequency_source_port_magnetic" if side == "source" else "frequency_port_magnetic"][i]
                raw = np.array([[_power_inner(port_E[j, k], port_H[j, k], electric[j], magnetic[j])
                                 for k in range(2)] for j in range(n)])
                template_error = max(template_error, np.max(np.abs(raw / np.linalg.norm(raw) - target)))
        compare(metrics, "frequency_complete_S_replay", replay_error, 1e-8)
        compare(metrics, "zero_structure_background_AH", baseline_error, 1e-5)
        compare(metrics, "fixed_physical_templates", template_error, 1e-10)
        compare(metrics, "AH_launch_k", np.max(data["frequency_launch_k_max_relative_error"]), 1e-8)
        compare(metrics, "AH_fixed_plane", np.max(data["frequency_plane_error_max_km"]), 1e-6)
        compare(metrics, "frequency_power_budget", np.max(np.abs(
            1 - data["frequency_reflected_power"] - data["frequency_transmitted_power"]
            - data["frequency_absorbed_power"])), 1e-5)

        slow_time = data["slow_time_s"]
        midpoint = len(slow_time) // 2
        for axis, period in (("u", 8.), ("v", 10.)):
            side = data["slow_delta_" + axis]
            predicted = side[midpoint] * np.exp(-2j * np.pi * slow_time / period)
            compare(metrics, f"slow_{axis}_translation_relative", np.max(np.abs(side - predicted))
                    / abs(side[midpoint]), .01)
            measured = doppler_hz(side, slow_time)[midpoint]
            compare(metrics, f"slow_{axis}_Doppler_sign_and_value_hz", abs(measured - 1 / period), .005)
        compare(metrics, "slow_density_recomputed", float(np.array_equal(
            data["slow_density"][0], data["slow_density"][-1])), 0.)
        compare(metrics, "slow_background_stationary", np.max(np.abs(
            data["slow_H0"] - data["slow_H0"][midpoint])), 1e-8)
        compare(metrics, "delay_parseval", abs(np.sum(data["power_delay_W_per_W"])
                - np.mean(np.abs(H[:len(frequencies)])**2)), 1e-12)
        compare(metrics, "doppler_parseval", abs(np.sum(data["power_doppler_W_per_W"])
                - np.mean(np.abs(data["slow_H"])**2)), 1e-12)
        compare(metrics, "delay_resolution_recorded", abs(summary["delay_fft_bin_s"]
                - 1 / (len(frequencies) * (frequencies[1] - frequencies[0]))), 1e-15)
        compare(metrics, "delay_broadening_marked_unresolved", float(
            summary["delay_broadening_status"] != "UNRESOLVED_AT_FINITE_BANDWIDTH"), 0.)

        # The independent complex exponentials test the signs of both transforms.
        delay_axis, _ = power_delay(frequencies, np.ones(len(frequencies), complex))
        target_delay = 2 / (len(frequencies) * (frequencies[1] - frequencies[0]))
        pure_delay = np.exp(2j * np.pi * (frequencies - case.frequency_hz) * target_delay)
        _, delay_power = power_delay(frequencies, pure_delay)
        compare(metrics, "positive_delay_transform_sign_s", abs(delay_axis[np.argmax(delay_power)] - target_delay), 1e-14)
        compare(metrics, "positive_delay_phase_slope_s", abs(phase_delay_s(
            np.exp(2j * np.pi * 1000 * target_delay), np.exp(-2j * np.pi * 1000 * target_delay), 1000.)
            - target_delay), 1e-14)
        target_doppler = 2 / (len(slow_time) * case.slow_time_step_s)
        pure_doppler = np.exp(-2j * np.pi * target_doppler * slow_time)
        axis, power = power_doppler(slow_time, pure_doppler)
        compare(metrics, "positive_doppler_transform_sign_hz", abs(axis[np.argmax(power)] - target_doppler), 1e-12)
        compare(metrics, "positive_doppler_phase_slope_hz", abs(
            doppler_hz(pure_doppler, slow_time)[midpoint] - target_doppler), 1e-12)
        failed = [name for name, metric in metrics.items() if not metric["pass"]]
        status, detail = ("PASS", None) if not failed else ("FAIL", ", ".join(failed))
    except Exception as exc:
        status, detail = "STOPPED", str(exc)
    output.mkdir(parents=True, exist_ok=True)
    (output / "validation.json").write_text(json.dumps(
        {"status": status, "detail": detail, "metrics": metrics}, indent=2), encoding="utf-8")
    print(f"S4a {status}: {sum(item['pass'] for item in metrics.values())}/{len(metrics)} checks")
    if detail:
        print(detail)
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
