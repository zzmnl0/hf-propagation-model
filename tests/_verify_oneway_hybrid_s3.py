"""S3 periodic full-wave scatter, independent sideband and AH/RX checks."""

from __future__ import annotations

import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import psutil
import scipy
from scipy.constants import c

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _verify_oneway_hybrid_s2 import (
    AFTER, BASIS, BEFORE, COEFFICIENT, EARTH_RADIUS_KM, ENTRY, F0, FIELD,
    PERIODS, THICKNESS, TIME_S, UniformProvider, make_models,
)
from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from ray_wave_oneway.local_wave import MediumSnapshot
from ray_wave_oneway.config import ONEWAY_AH_FD_SETTINGS
from ray_wave_oneway.periodic import solve_periodic_slab, apply_periodic_window
from ray_wave_oneway.propagation import make_field_template, prepare_windows, propagate_oneway, reconstruct_fields
from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF, group_velocity_ah3d_ecef


OUT = Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "isolation" / "s3"
ETA = .05
START = time.monotonic()
PROCESS = psutil.Process()


def indices_for(size):
    modes = np.arange(-(size // 2), size // 2 + 1)
    return np.array([(u, v) for u in modes for v in modes])


def compare(metrics, name, value, limit):
    number = float(value)
    metrics[name] = {"value": number, "limit": float(limit),
                     "pass": bool(np.isfinite(number) and number <= limit)}


def prepare(size, field=FIELD):
    if np.array_equal(field, FIELD):
        models = make_models(F0)
    else:
        models = {branch: RefractiveIndexAH3DECEF(
            UniformProvider(), earth_radius_km=EARTH_RADIUS_KM, freq_MHz=F0 / 1e6,
            wave_mode=branch, geomagnetic_provider=ConstantGeomagneticField3D(field),
            plasma_coefficient=COEFFICIENT,
            integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
        ) for branch in ("O", "X")}
    medium = MediumSnapshot.from_ah_model(models["O"], ENTRY, TIME_S)
    indices = indices_for(size)
    window = prepare_windows(models["O"], medium, ENTRY, BASIS, THICKNESS,
                             PERIODS, indices, TIME_S)
    return models, medium, window


def density_grid(medium, size, eta, both=True):
    grid = np.arange(size)
    u_term = np.cos(2 * np.pi * grid[:, None] / size)
    v_term = .6 * np.cos(2 * np.pi * grid[None, :] / size) if both else 0.0
    return medium.Ne_m3 * (1 + eta * (u_term + v_term)) * np.ones((1, size))


def uniform_s1_matrix(window):
    n = len(window["slabs"])
    matrix = np.zeros((4 * n, 4 * n), complex)
    for i, slab in enumerate(window["slabs"]):
        block = [2 * i, 2 * i + 1, 2 * n + 2 * i, 2 * n + 2 * i + 1]
        matrix[np.ix_(block, block)] = slab.S
    return matrix


def s2_exterior(models, window):
    indices = window["indices"]
    support = np.all(np.abs(indices) <= 1, axis=1)
    spatial = (np.exp(-.45 * np.sum(indices**2, axis=1))
               * np.exp(1j * (.23 * indices[:, 0] - .17 * indices[:, 1])) * support)
    source = np.column_stack((spatial, .45j * spatial))
    source /= np.linalg.norm(source)
    k_s = np.array([[port.k_suv_rad_m[0] for port in slab.right_ports.forward]
                    for slab in window["slabs"]])
    receiver = source * np.exp(1j * k_s * (BEFORE + THICKNESS + AFTER))
    source_template = make_field_template(window, source, "source")
    receiver_template = make_field_template(window, receiver, "receiver")
    return propagate_oneway(window, models, source_template, receiver_template,
                            BEFORE, AFTER, TIME_S)


def born_sideband(window, periodic, eta, shift):
    indices = window["indices"]
    center = int(np.flatnonzero(np.all(indices == [0, 0], axis=1))[0])
    side = int(np.flatnonzero(np.all(indices == [shift, 0], axis=1))[0])
    source_modes = window["slabs"][center].left_ports.forward
    output_modes = window["slabs"][side].right_ports.forward
    source_index = int(np.argmax([abs(mode.electric_suv[2]) for mode in source_modes]))
    output_index = int(np.argmax([abs(mode.electric_suv[2]) for mode in output_modes]))
    source = source_modes[source_index]
    output = output_modes[output_index]
    q0, q1 = source.k_suv_rad_m[0], output.k_suv_rad_m[0]
    k0 = 2 * np.pi * F0 / c
    epsilon_fourier = -.2 * eta / 2
    expected = (1j * k0**2 * epsilon_fourier * THICKNESS / (2 * q1)
                * np.exp(1j * (q0 + q1) * THICKNESS / 2)
                * np.sinc((q0 - q1) * THICKNESS / (2 * np.pi))
                * source.electric_suv[2] / output.electric_suv[2])
    n = len(indices)
    row, column = 2 * n + 2 * side + output_index, 2 * center + source_index
    actual = periodic.S[row, column]
    noise = abs(periodic.S0[row, column])
    return actual, expected, noise


def write_result(metrics, summary, arrays, status, detail=None):
    OUT.mkdir(parents=True, exist_ok=True)
    memory = PROCESS.memory_info()
    payload = {"status": status, "detail": detail, "metrics": metrics, "summary": summary,
               "resource": {"available_bytes": psutil.virtual_memory().available,
                            "rss_bytes": memory.rss,
                            "peak_rss_bytes": getattr(memory, "peak_wset", memory.rss),
                            "elapsed_s": time.monotonic() - START}}
    (OUT / "metrics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if arrays:
        np.savez_compressed(OUT / "representative_fields.npz", **arrays)
    manifest = {
        "task": "AH3D-ONEWAY-HYBRID-20261105-S3", "status": status,
        "command": "python -B tests/_verify_oneway_hybrid_s3.py",
        "python_executable": sys.executable, "python_version": sys.version.split()[0],
        "numpy_version": np.__version__, "scipy_version": scipy.__version__,
        "inputs": {"frequency_hz": F0, "medium_time_s": TIME_S,
                   "coordinate_model": "geocentric_spherical_ecef", "earth_radius_km": EARTH_RADIUS_KM,
                   "entrance_ecef_km": ENTRY.tolist(), "basis_suv_ecef_columns": BASIS.tolist(),
                   "background_Ne_m3": .2 * F0**2 / COEFFICIENT, "plasma_coefficient": COEFFICIENT,
                   "B_ecef_T": FIELD.tolist(), "periods_m": PERIODS.tolist(),
                   "thickness_m": THICKNESS, "before_distance_m": BEFORE, "after_distance_m": AFTER,
                   "structure": "Ne_bg * (1 + eta*(cos(2*pi*u/Lu)+0.6*cos(2*pi*v/Lv)))",
                   "eta": ETA, "lossy_collision_s": 1e5, "grid_sizes": [5, 7],
                   "born_eta": [4e-5, 2e-5], "reference_input_power_W": 1.0},
        "outputs": {"metrics.json": "scattering, conservation, convergence and ray/RX checks",
                    "representative_fields.npz": sorted(arrays)},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main():
    metrics, summary, arrays = {}, {}, {}
    try:
        print("S3 constructing 5x5 periodic Maxwell operators", flush=True)
        models5, medium5, window5 = prepare(5)
        uniform5 = solve_periodic_slab(window5, medium5, np.full((5, 5), medium5.Ne_m3))
        expected = uniform_s1_matrix(window5)
        compare(metrics, "uniform_S_vs_S1", np.max(np.abs(uniform5.S - expected)), 1e-6)
        compare(metrics, "uniform_S0_vs_S1", np.max(np.abs(uniform5.S0 - expected)), 1e-6)
        compare(metrics, "uniform_S_vs_S0", np.max(np.abs(uniform5.S - uniform5.S0)), 1e-8)
        density5 = density_grid(medium5, 5, ETA)
        real5 = solve_periodic_slab(window5, medium5, density5)
        loss5 = solve_periodic_slab(window5, replace(medium5, collision_s=1e5), density5)
        identity = np.eye(4 * real5.count)
        compare(metrics, "lossless_S_unitarity", np.linalg.norm(real5.S.conj().T @ real5.S - identity, ord=2), 1e-5)
        compare(metrics, "lossy_S_passivity", max(0., np.linalg.eigvalsh(loss5.S.conj().T @ loss5.S - identity).max()), 1e-5)
        rng = np.random.default_rng(20261105)
        coherent = rng.normal(size=4 * real5.count) + 1j * rng.normal(size=4 * real5.count)
        coherent /= np.linalg.norm(coherent)
        for name, slab in (("lossless", real5), ("lossy", loss5)):
            outgoing = slab.S @ coherent
            absorbed = slab.absorbed_power_W(coherent, samples=65)
            coarse = slab.absorbed_power_W(coherent, samples=33)
            compare(metrics, name + "_coherent_RTA", abs(1 - np.vdot(outgoing, outgoing).real - absorbed), 1e-5)
            compare(metrics, name + "_absorption_convergence", abs(absorbed - coarse), 1e-5)

        print("S3 tracing 5x5 AH ports", flush=True)
        ah5 = s2_exterior(models5, window5)
        print("S3 tracing reversed-B AH context", flush=True)
        reverse_models, _, reverse_window = prepare(5, -FIELD)
        reverse_ah = s2_exterior(reverse_models, reverse_window)
        k_errors = []
        polarization_distances = []
        for regular_slab, reverse_slab in zip(window5["slabs"], reverse_window["slabs"]):
            for regular, reverse in zip(regular_slab.left_ports.forward, reverse_slab.left_ports.forward):
                k_errors.append(np.linalg.norm(regular.k_ecef_rad_km - reverse.k_ecef_rad_km)
                                / np.linalg.norm(regular.k_ecef_rad_km))
                overlap = abs(np.vdot(regular.electric_suv, reverse.electric_suv))
                overlap /= np.linalg.norm(regular.electric_suv) * np.linalg.norm(reverse.electric_suv)
                polarization_distances.append(max(0., 1 - overlap))
        compare(metrics, "AH_reverse_B_equal_wavevector", max(k_errors), 1e-8)
        compare(metrics, "AH_reverse_B_distinct_polarization", max(
            0., 1e-3 - max(polarization_distances)), 0.)
        try:
            apply_periodic_window(real5, window5, reverse_ah)
        except ValueError:
            context_rejected = True
        else:
            context_rejected = False
        compare(metrics, "AH_reverse_B_context_rejected", float(not context_rejected), 0.)
        zero = apply_periodic_window(uniform5, window5, ah5)
        actual5 = apply_periodic_window(real5, window5, ah5)
        lossy5 = apply_periodic_window(loss5, window5, ah5)
        try:
            apply_periodic_window(real5, dict(window5), ah5)
        except ValueError:
            mismatch_rejected = True
        else:
            mismatch_rejected = False
        compare(metrics, "AH_wrong_window_rejected", float(not mismatch_rejected), 0.)
        compare(metrics, "AH_zero_H_vs_S2", abs(zero["H_rx"] - ah5["H_rx_background"]), 1e-5)
        compare(metrics, "AH_zero_fields_vs_S2", np.max(np.abs(zero["transmitted"] - ah5["transmitted_background"])), 1e-5)
        for name, outcome in (("main_lossless", actual5), ("main_lossy", lossy5)):
            compare(metrics, name + "_RTA", abs(outcome["incident_power_W"] - outcome["reflected_power_W"]
                    - outcome["transmitted_power_W"] - outcome["absorbed_power_W"]), 1e-5)
            compare(metrics, name + "_receiver_budget", max(0., -outcome["uncaptured_power_W"]), 1e-10)
        electric, magnetic, flux = reconstruct_fields(window5, actual5["transmitted"], "receiver", 12)
        compare(metrics, "main_RX_field_flux", abs(flux - actual5["transmitted_power_W"]), 1e-5)
        indices5 = window5["indices"]
        new_bin = int(np.flatnonzero(np.all(indices5 == [2, 0], axis=1))[0])
        new_j = int(np.argmax(np.abs(actual5["transmitted"][new_bin])))
        new_value = actual5["transmitted"][new_bin, new_j]
        new_background = zero["transmitted"][new_bin, new_j]
        compare(metrics, "new_channel_above_background", max(
            0., 1e-8 - (abs(new_value) - abs(new_background))), 0.)
        new_port = window5["slabs"][new_bin].right_ports.forward[new_j]
        state = np.r_[window5["exit_ecef_km"], new_port.k_ecef_rad_km, 0.]
        velocity = group_velocity_ah3d_ecef(state, models5[new_port.label], t=TIME_S)
        group_direction = velocity / np.linalg.norm(velocity)
        central = int(np.flatnonzero(np.all(indices5 == [0, 0], axis=1))[0])
        converted_bin = int(np.flatnonzero(np.all(indices5 == [1, 0], axis=1))[0])
        central_o = next(i for i, port in enumerate(window5["slabs"][central].left_ports.forward)
                         if port.label == "O")
        side_x = next(i for i, port in enumerate(window5["slabs"][converted_bin].right_ports.forward)
                      if port.label == "X")
        row, col = 2 * real5.count + 2 * converted_bin + side_x, 2 * central + central_o
        converted = real5.S[row, col]
        converted_background = real5.S0[row, col]
        compare(metrics, "new_cross_bin_O_to_X", max(
            0., 1e-8 - (abs(converted) - abs(converted_background))), 0.)

        print("S3 constructing 7x7 and tracing AH ports", flush=True)
        models7, medium7, window7 = prepare(7)
        real7 = solve_periodic_slab(window7, medium7, density_grid(medium7, 7, ETA))
        ah7 = s2_exterior(models7, window7)
        actual7 = apply_periodic_window(real7, window7, ah7)
        common7 = [int(np.flatnonzero(np.all(window7["indices"] == index, axis=1))[0])
                   for index in window5["indices"]]
        compare(metrics, "grid_source_template_fixed", np.max(np.abs(
            ah7["source"][common7] - ah5["source"])), 1e-10)
        compare(metrics, "grid_receiver_template_fixed", np.max(np.abs(
            ah7["receiver"][common7] - ah5["receiver"])), 1e-10)
        compare(metrics, "grid_RX_amplitude_relative", abs(actual7["H_rx"] - actual5["H_rx"]) / abs(actual7["H_rx"]), .01)
        high_mask = np.any(np.abs(window7["indices"]) > 2, axis=1)
        high_power = float(np.vdot(actual7["transmitted"][high_mask], actual7["transmitted"][high_mask]).real)
        compare(metrics, "grid_new_high_order_power_W", high_power, .001)
        compare(metrics, "grid7_RTA", abs(actual7["incident_power_W"] - actual7["reflected_power_W"]
                - actual7["transmitted_power_W"] - actual7["absorbed_power_W"]), 1e-5)

        print("S3 independent B=0 Born sidebands", flush=True)
        _, unmagnetized, born_window = prepare(5, np.zeros(3))
        born_errors = []
        born_data = []
        for eta in (4e-5, 2e-5):
            periodic = solve_periodic_slab(born_window, unmagnetized, density_grid(unmagnetized, 5, eta, both=False))
            for shift in (-1, 1):
                value, prediction, noise = born_sideband(born_window, periodic, eta, shift)
                error = abs(value - prediction)
                compare(metrics, f"born_relative_{eta}_{shift}", error / abs(prediction), .01)
                born_errors.append(error)
                born_data.append((value, prediction, noise))
        compare(metrics, "born_half_perturbation_decreases_error", max(
            0., max(born_errors[2:]) - max(born_errors[:2])), 0.)
        noise_floor = max(20 * np.finfo(float).eps, 10 * max(item[2] for item in born_data))
        compare(metrics, "born_half_perturbation_noise_floor", max(born_errors[2:]), noise_floor)

        summary.update({
            "frequency_hz": F0, "medium_time_s": TIME_S, "grid_sizes": [5, 7],
            "eta": ETA, "AH_rays_5": len(ah5["rays"]), "AH_rays_7": len(ah7["rays"]),
            "AH_rays_reverse_B": len(reverse_ah["rays"]),
            "reverse_B_k_max_relative_error": max(k_errors),
            "reverse_B_polarization_max_distance": max(polarization_distances),
            "matching_condition_5": real5.matching_condition, "matching_condition_7": real7.matching_condition,
            "H_rx_5": [float(actual5["H_rx"].real), float(actual5["H_rx"].imag)],
            "H_rx_7": [float(actual7["H_rx"].real), float(actual7["H_rx"].imag)],
            "H_rx_background": [float(zero["H_rx"].real), float(zero["H_rx"].imag)],
            "lossless_power_W": {key: actual5[key] for key in ("incident_power_W", "reflected_power_W", "transmitted_power_W", "absorbed_power_W", "rx_power_W", "uncaptured_power_W")},
            "lossy_power_W": {key: lossy5[key] for key in ("incident_power_W", "reflected_power_W", "transmitted_power_W", "absorbed_power_W", "rx_power_W", "uncaptured_power_W")},
            "new_bin": [2, 0], "new_bin_mode": new_port.label,
            "new_bin_complex_amplitude": [float(new_value.real), float(new_value.imag)],
            "new_bin_background_abs": float(abs(new_background)),
            "new_bin_k_ecef_rad_km": new_port.k_ecef_rad_km.tolist(),
            "new_bin_group_direction_ecef": group_direction.tolist(),
            "cross_bin_O_to_X_amplitude": [float(converted.real), float(converted.imag)],
            "cross_bin_O_to_X_background_abs": float(abs(converted_background)),
            "grid7_new_high_order_power_W": high_power,
        })
        arrays.update({"indices_5": indices5, "indices_7": window7["indices"],
                       "S_5": real5.S, "S0_5": real5.S0, "S_lossy_5": loss5.S,
                       "S_7": real7.S, "S0_7": real7.S0,
                       "density_5_m3": density5, "density_7_m3": density_grid(medium7, 7, ETA),
                       "source_5": ah5["source"], "receiver_5": ah5["receiver"],
                       "phase_before_5": ah5["phase_before_rad"], "phase_after_5": ah5["phase_after_rad"],
                       "tau_before_5": ah5["tau_before_s"], "tau_after_5": ah5["tau_after_s"],
                       "transmitted_5": actual5["transmitted"], "transmitted_background_5": zero["transmitted"],
                       "reflected_5": actual5["reflected"], "transmitted_7": actual7["transmitted"],
                       "RX_E_suv": electric, "RX_H_suv": magnetic,
                       "born_values": np.array([x[0] for x in born_data]),
                       "born_predictions": np.array([x[1] for x in born_data]),
                       "born_background_leakage": np.array([x[2] for x in born_data])})
        failed = [name for name, item in metrics.items() if not item["pass"]]
        status = "PASS" if not failed else "FAIL"
        write_result(metrics, summary, arrays, status, None if not failed else ", ".join(failed))
        print(f"S3 {status}: {len(metrics) - len(failed)}/{len(metrics)} checks")
        for name in failed:
            print(f"FAIL {name}: {metrics[name]}")
        return 0 if not failed else 1
    except Exception as exc:
        write_result(metrics, summary, arrays, "STOPPED", str(exc))
        print(f"S3 STOPPED: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
