"""Three-group S2 validation on a finite periodic beam and real AH rays."""

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

from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from ray_wave_oneway.local_wave import MediumSnapshot
from ray_wave_oneway.config import ONEWAY_AH_FD_SETTINGS
from models.ne_provider_geocentric import GeocentricNeProvider3D, WGS84NeProviderAdapter
from ray_wave_oneway.propagation import (
    _power_inner,
    _trace_to_plane,
    make_field_template,
    prepare_windows,
    project_field_template,
    propagate_oneway,
    reconstruct_fields,
)
from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF, WGS84_COORDINATE_MODEL


OUT = Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "isolation" / "s2"
F0 = 10e6
DF = 1e3
COEFFICIENT = 80.61638587963631
EARTH_RADIUS_KM = 6371.0
TIME_S = 37.0
ENTRY = np.array([EARTH_RADIUS_KM + 250.0, 0., 0.])
BASIS = np.eye(3)
PERIODS = np.array([500., 600.])
THICKNESS = 2 * c / F0
BEFORE = 100.0
AFTER = 100.0
FIELD = 45e-6 * np.ones(3) / np.sqrt(3)
INDICES = np.array([(u, v) for u in range(-2, 3) for v in range(-2, 3)])
START = time.monotonic()
PROCESS = psutil.Process()


class UniformProvider(GeocentricNeProvider3D):
    earth_radius_km = EARTH_RADIUS_KM

    def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
        return .2 * F0**2 / COEFFICIENT

    def get_metadata(self):
        return {**super().get_metadata(), "earth_radius_km": self.earth_radius_km}


def make_models(frequency):
    provider = UniformProvider()
    field = ConstantGeomagneticField3D(FIELD)
    return {branch: RefractiveIndexAH3DECEF(
        provider, earth_radius_km=EARTH_RADIUS_KM, freq_MHz=frequency / 1e6,
        wave_mode=branch, geomagnetic_provider=field, plasma_coefficient=COEFFICIENT,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
    ) for branch in ("O", "X")}


def setup(frequency, *, background_layer=False):
    models = make_models(frequency)
    midpoint = ENTRY + BASIS[:, 0] * THICKNESS / 2000
    background = MediumSnapshot.from_ah_model(models["O"], midpoint, TIME_S)
    layer = (background if background_layer else replace(
        background, Ne_m3=.3 * F0**2 / COEFFICIENT, collision_s=1e5,
    ))
    return models, prepare_windows(models["O"], layer, ENTRY, BASIS, THICKNESS,
                                   PERIODS, INDICES, TIME_S)


def compare(metrics, name, value, limit):
    number = float(value)
    metrics[name] = {"value": number, "limit": float(limit),
                     "pass": bool(np.isfinite(number) and number <= limit)}


def unit_phase_slope(plus, minus):
    return float(np.angle(plus * np.conj(minus)) / (4 * np.pi * DF))


def geometric_expansion(window, source, use_background):
    output = []
    for i, slab in enumerate(window["slabs"]):
        left = np.array([mode.k_suv_rad_m[0] for mode in slab.left_ports.forward])
        right = np.array([mode.k_suv_rad_m[0] for mode in slab.right_ports.forward])
        matrix = slab.S0 if use_background else slab.S
        output.append(np.exp(1j * right * AFTER) *
                      (matrix[2:, :2] @ (source[i] * np.exp(1j * left * BEFORE))))
    return np.asarray(output)


def rephased_response(window, run, source_template, receiver_template, phases_in, phases_out):
    a_raw, b_raw = [], []
    transformed = []
    for i, slab in enumerate(window["slabs"]):
        left = slab.left_ports.forward
        right = slab.right_ports.forward
        a_raw.append([_power_inner(
            mode.electric_suv * np.exp(1j * phases_in[i, j]),
            mode.magnetic_suv * np.exp(1j * phases_in[i, j]),
            source_template["electric"][i], source_template["magnetic"][i],
        ) for j, mode in enumerate(left)])
        b_raw.append([_power_inner(
            mode.electric_suv * np.exp(1j * phases_out[i, j]),
            mode.magnetic_suv * np.exp(1j * phases_out[i, j]),
            receiver_template["electric"][i], receiver_template["magnetic"][i],
        ) for j, mode in enumerate(right)])
        transformed.append(np.exp(-1j * phases_out[i, :, None]) * slab.S[2:, :2]
                           * np.exp(1j * phases_in[i, None, :]))
    a_raw = np.asarray(a_raw)
    b_raw = np.asarray(b_raw)
    a = a_raw / np.linalg.norm(a_raw)
    b = b_raw / np.linalg.norm(b_raw)
    output = np.array([np.exp(1j * run["phase_after_rad"][i]) *
                       (transformed[i] @ (a[i] * np.exp(1j * run["phase_before_rad"][i])))
                       for i in range(len(transformed))])
    return np.vdot(b, output), a, b


def write_result(metrics, summary, arrays, status, detail=None):
    OUT.mkdir(parents=True, exist_ok=True)
    memory = PROCESS.memory_info()
    data = {"status": status, "detail": detail, "metrics": metrics, "summary": summary,
            "resource": {"available_bytes": psutil.virtual_memory().available,
                         "rss_bytes": memory.rss,
                         "peak_rss_bytes": getattr(memory, "peak_wset", memory.rss),
                         "elapsed_s": time.monotonic() - START}}
    (OUT / "metrics.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    if arrays:
        np.savez_compressed(OUT / "representative_fields.npz", **arrays)
    manifest = {
        "task": "AH3D-ONEWAY-HYBRID-20261105-S2",
        "status": status,
        "command": "python -B tests/_verify_oneway_hybrid_s2.py",
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input": {
            "frequencies_hz": [F0 - DF, F0, F0 + DF],
            "medium_time_s": TIME_S,
            "coordinate_model": "geocentric_spherical_ecef",
            "earth_radius_km": EARTH_RADIUS_KM,
            "plasma_coefficient": COEFFICIENT,
            "background_Ne_m3": .2 * F0**2 / COEFFICIENT,
            "slab_Ne_m3": .3 * F0**2 / COEFFICIENT,
            "B_ecef_T": FIELD.tolist(), "background_collision_s": 0.0,
            "slab_collision_s": 1e5, "entrance_ecef_km": ENTRY.tolist(),
            "basis_suv_ecef_columns": BASIS.tolist(), "periods_m": PERIODS.tolist(),
            "spectral_indices": INDICES.tolist(), "slab_thickness_m": THICKNESS,
            "before_distance_m": BEFORE, "after_distance_m": AFTER,
            "reference_input_power_W": 1.0,
            "reconstructed_field_grid": "12x12 on one periodic cell",
        },
        "outputs": {"metrics.json": "validation metrics, powers, delays and resource measurement",
                    "representative_fields.npz": sorted(arrays)},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main():
    metrics, summary, arrays = {}, {}, {}
    try:
        frequencies = [F0 - DF, F0, F0 + DF]
        setups = [setup(frequency) for frequency in frequencies]
        center_window = setups[1][1]
        spatial = np.exp(-0.45 * np.sum(INDICES**2, axis=1)) * np.exp(
            1j * (0.23 * INDICES[:, 0] - 0.17 * INDICES[:, 1]))
        source_coeff = np.column_stack((spatial, 0.45j * spatial))
        source_coeff /= np.linalg.norm(source_coeff)
        k_s = np.array([[port.k_suv_rad_m[0] for port in slab.right_ports.forward]
                        for slab in center_window["slabs"]])
        receiver_coeff = source_coeff * np.exp(1j * k_s * (BEFORE + THICKNESS + AFTER))
        source_template = make_field_template(center_window, source_coeff, "source")
        receiver_template = make_field_template(center_window, receiver_coeff, "receiver")
        runs = []
        for frequency, (models, window) in zip(frequencies, setups):
            print(f"S2 tracing {frequency / 1e6:.6f} MHz", flush=True)
            runs.append(propagate_oneway(window, models, source_template, receiver_template,
                                         BEFORE, AFTER, TIME_S))
        run = runs[1]
        window = center_window

        # A: actual L2 AH launches, target-plane intersection, and phase.
        compare(metrics, "A_launch_k_max_relative", max(r["launch_k_max_relative_error"] for r in runs), 1e-8)
        compare(metrics, "A_target_plane_max_km", max(r["target_plane_max_error_km"] for r in runs), 1e-6)
        compare(metrics, "A_uniform_plane_phase_max_rad", max(r["uniform_phase_max_error_rad"] for r in runs), 1e-5)
        wgs_provider = WGS84NeProviderAdapter(UniformProvider())
        wgs_model = RefractiveIndexAH3DECEF(
            wgs_provider, earth_radius_km=EARTH_RADIUS_KM, freq_MHz=F0 / 1e6,
            wave_mode="O", geomagnetic_provider=ConstantGeomagneticField3D(FIELD),
            plasma_coefficient=COEFFICIENT, coordinate_model=WGS84_COORDINATE_MODEL,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
        )
        wgs_background = MediumSnapshot.from_ah_model(wgs_model, ENTRY, TIME_S)
        wgs_window = prepare_windows(wgs_model, wgs_background, ENTRY, BASIS, THICKNESS,
                                     PERIODS, np.array([[0, 0]]), TIME_S)
        wgs_ray = _trace_to_plane(wgs_model, wgs_window["slabs"][0].left_ports.forward[0],
                                  ENTRY - BASIS[:, 0] * BEFORE / 1000, ENTRY, BASIS,
                                  TIME_S, wgs_window["left_medium"])
        compare(metrics, "A_WGS84_launch_k_relative", wgs_ray["launch_k_relative_error"], 1e-8)
        compare(metrics, "A_WGS84_plane_error_km", abs(wgs_ray["terminal_plane_error_km"]), 1e-6)
        compare(metrics, "A_WGS84_phase_error_rad", abs(wgs_ray["phase_rad"] - wgs_ray["uniform_phase_rad"]), 1e-5)

        # B: same-background replacement and independently reconstructed flux.
        _, uniform_window = setup(F0, background_layer=True)
        compare(metrics, "B_same_layer_S_equals_S0", max(
            np.max(np.abs(s.S - s.S0)) for s in uniform_window["slabs"]), 1e-8)
        expected_background = geometric_expansion(window, run["source"], True)
        compare(metrics, "B_background_field_vs_analytic", np.max(np.abs(
            run["transmitted_background"] - expected_background)), 1e-5)
        compare(metrics, "B_background_rx_vs_analytic", abs(
            run["H_rx_background"] - np.vdot(run["receiver"], expected_background)), 1e-5)
        for side, coeff in (("source", run["source"]), ("receiver", run["transmitted"])):
            electric, magnetic, p_coarse = reconstruct_fields(window, coeff, side, 12)
            _, _, p_fine = reconstruct_fields(window, coeff, side, 24)
            p_spectrum = float(np.vdot(coeff, coeff).real)
            compare(metrics, f"B_{side}_flux_vs_spectrum", abs(p_coarse - p_spectrum) / p_spectrum, 1e-5)
            compare(metrics, f"B_{side}_grid_refinement", abs(p_fine - p_coarse) / p_spectrum, 1e-5)
            arrays[side + "_E_suv"] = electric
            arrays[side + "_H_suv"] = magnetic
        grid = 12
        u = np.arange(grid) * PERIODS[0] / grid
        v = np.arange(grid) * PERIODS[1] / grid
        uu, vv = np.meshgrid(u, v, indexing="ij")
        harmonics = np.exp(1j * (window["tangent_k_rad_m"][:, 0, None, None] * uu
                                  + window["tangent_k_rad_m"][:, 1, None, None] * vv))
        area = np.prod(PERIODS)
        rx_electric = np.einsum("nuv,nd->uvd", harmonics, receiver_template["electric"]) / np.sqrt(area)
        rx_magnetic = np.einsum("nuv,nd->uvd", harmonics, receiver_template["magnetic"]) / np.sqrt(area)
        out_electric, out_magnetic = arrays["receiver_E_suv"], arrays["receiver_H_suv"]
        inner = (np.cross(np.conj(rx_electric), out_magnetic)
                 + np.cross(out_electric, np.conj(rx_magnetic)))[..., 0] / 4
        physical_rx = area * inner.mean() / run["receiver_projection_norm"]
        compare(metrics, "B_rx_physical_field_projection", abs(physical_rx - run["H_rx"]), 1e-5)
        try:
            reconstruct_fields(window, run["source"], "source", 4)
        except ValueError:
            alias_rejected = True
        else:
            alias_rejected = False
        compare(metrics, "B_aliasing_grid_rejected", float(not alias_rejected), 0.)
        compare(metrics, "B_RTA_closure", abs(
            run["incident_power_W"] - run["reflected_power_W"] - run["transmitted_power_W"]
            - run["absorbed_power_W"]), 1e-5)
        compare(metrics, "B_nonnegative_uncaptured", max(0., -run["uncaptured_power_W"]), 1e-10)
        class VaryingProvider(UniformProvider):
            def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
                return super().get_Ne_geocentric(lat_deg, lon_deg, height_km, t) * (1 + .001 * height_km)

        varying_model = RefractiveIndexAH3DECEF(
            VaryingProvider(), earth_radius_km=EARTH_RADIUS_KM, freq_MHz=F0 / 1e6,
            geomagnetic_provider=ConstantGeomagneticField3D(FIELD), plasma_coefficient=COEFFICIENT,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
        )
        varying_layer = MediumSnapshot.from_ah_model(varying_model, ENTRY, TIME_S)
        try:
            prepare_windows(varying_model, varying_layer, ENTRY, BASIS, THICKNESS,
                            PERIODS, np.array([[0, 0]]), TIME_S)
        except ValueError:
            nonuniform_rejected = True
        else:
            nonuniform_rejected = False
        compare(metrics, "B_nonuniform_medium_rejected", float(not nonuniform_rejected), 0.)

        # C: physical layer, fixed-template projection, and frequency slopes.
        expected = geometric_expansion(window, run["source"], False)
        expected_H = np.vdot(run["receiver"], expected)
        compare(metrics, "C_layer_field_vs_analytic", np.max(np.abs(run["transmitted"] - expected)), 1e-5)
        compare(metrics, "C_layer_rx_vs_analytic", abs(run["H_rx"] - expected_H), 1e-5)
        compare(metrics, "C_layer_transmitted_power_vs_analytic", abs(
            run["transmitted_power_W"] - np.vdot(expected, expected).real), 1e-5)
        compare(metrics, "C_layer_rx_power_vs_analytic", abs(run["rx_power_W"] - abs(expected_H)**2), 1e-5)
        center_bin = int(np.flatnonzero(np.all(INDICES == 0, axis=1))[0])
        center_port = window["slabs"][center_bin].left_ports.forward[0]
        local_values, full_values = [], []
        for (_, current), current_run in zip(setups, runs):
            slab = current["slabs"][center_bin]
            current_left = slab.left_ports.forward[0]
            current_right = slab.right_ports.forward[0]
            a = _power_inner(current_left.electric_suv, current_left.magnetic_suv,
                             center_port.electric_suv, center_port.magnetic_suv)
            b = _power_inner(current_right.electric_suv, current_right.magnetic_suv,
                             center_port.electric_suv, center_port.magnetic_suv)
            a /= abs(a)
            b /= abs(b)
            local = np.conj(b) * slab.S[2, 0] * a
            phase = (current_run["phase_before_rad"][center_bin, 0]
                     + current_run["phase_after_rad"][center_bin, 0])
            local_values.append(local)
            full_values.append(local * np.exp(1j * phase))
        if min(abs(value) for value in local_values) < 1e-8:
            raise ValueError("selected single mode has negligible local transmission")
        local_tau = unit_phase_slope(local_values[2], local_values[0])
        channel_tau = unit_phase_slope(full_values[2], full_values[0])
        ah_tau = float(run["tau_before_s"][center_bin, 0] + run["tau_after_s"][center_bin, 0])
        compare(metrics, "C_single_channel_delay_closure_s", abs(channel_tau - (local_tau + ah_tau)),
                max(1e-9, .001 * abs(channel_tau)))
        rx_tau = unit_phase_slope(runs[2]["H_rx"], runs[0]["H_rx"])
        rng = np.random.default_rng(20261105)
        rephased = []
        for (_, current), current_run in zip(setups, runs):
            phase_in = rng.uniform(-np.pi, np.pi, (len(INDICES), 2))
            phase_out = rng.uniform(-np.pi, np.pi, (len(INDICES), 2))
            h, a, b = rephased_response(current, current_run, source_template, receiver_template,
                                        phase_in, phase_out)
            compare(metrics, f"C_rephase_H_{len(rephased)}", abs(h - current_run["H_rx"]), 1e-8)
            compare(metrics, f"C_rephase_source_{len(rephased)}", np.max(np.abs(
                a * np.exp(1j * phase_in) - current_run["source"])), 1e-8)
            compare(metrics, f"C_rephase_rx_{len(rephased)}", np.max(np.abs(
                b * np.exp(1j * phase_out) - current_run["receiver"])), 1e-8)
            rephased.append(h)
        compare(metrics, "C_rephase_rx_delay_s", abs(
            unit_phase_slope(rephased[2], rephased[0]) - rx_tau), 1e-12)
        compare(metrics, "C_source_template_center", np.max(np.abs(run["source"] - source_coeff)), 1e-10)
        compare(metrics, "C_receiver_template_center", np.max(np.abs(run["receiver"] - receiver_coeff)), 1e-10)

        summary.update({
            "frequencies_hz": frequencies, "medium_time_s": TIME_S,
            "angular_bins": len(INDICES), "rays_per_frequency": len(run["rays"]),
            "periods_m": PERIODS.tolist(), "thickness_m": THICKNESS,
            "before_distance_m": BEFORE, "after_distance_m": AFTER,
            "incident_power_W": run["incident_power_W"], "reflected_power_W": run["reflected_power_W"],
            "transmitted_power_W": run["transmitted_power_W"], "absorbed_power_W": run["absorbed_power_W"],
            "rx_power_W": run["rx_power_W"], "uncaptured_power_W": run["uncaptured_power_W"],
            "H_rx": [float(run["H_rx"].real), float(run["H_rx"].imag)],
            "H_rx_background": [float(run["H_rx_background"].real), float(run["H_rx_background"].imag)],
            "rx_coherent_delay_s": rx_tau, "single_channel_delay_s": channel_tau,
            "single_channel_AH_delay_s": ah_tau, "single_channel_local_delay_s": local_tau,
            "source_projection_norms": [r["source_projection_norm"] for r in runs],
            "receiver_projection_norms": [r["receiver_projection_norm"] for r in runs],
        })
        arrays.update({"indices": INDICES, "source_coefficients": run["source"],
                       "receiver_coefficients": run["receiver"], "transmitted": run["transmitted"],
                       "transmitted_background": run["transmitted_background"], "reflected": run["reflected"],
                       "S": np.array([s.S for s in window["slabs"]]),
                       "S0": np.array([s.S0 for s in window["slabs"]]),
                       "S_all_frequencies": np.array([
                           [slab.S for slab in current["slabs"]] for _, current in setups]),
                       "S0_all_frequencies": np.array([
                           [slab.S0 for slab in current["slabs"]] for _, current in setups]),
                       "phase_before_rad": run["phase_before_rad"], "phase_after_rad": run["phase_after_rad"],
                       "tau_before_s": run["tau_before_s"], "tau_after_s": run["tau_after_s"],
                       "ray_end_ecef_km": np.array([r["endpoint_ecef_km"] for r in run["rays"]]),
                       "ray_end_k_ecef_rad_km": np.array([r["endpoint_k_ecef_rad_km"] for r in run["rays"]]),
                       "H_rx_all_frequencies": np.array([r["H_rx"] for r in runs]),
                       "frequencies_hz": np.array(frequencies),
                       "source_coefficients_all_frequencies": np.array([r["source"] for r in runs]),
                       "receiver_coefficients_all_frequencies": np.array([r["receiver"] for r in runs]),
                       "source_ports_E_suv_all_frequencies": np.array([
                           [[port.electric_suv for port in slab.left_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "source_ports_H_suv_all_frequencies": np.array([
                           [[port.magnetic_suv for port in slab.left_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "receiver_ports_E_suv_all_frequencies": np.array([
                           [[port.electric_suv for port in slab.right_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "receiver_ports_H_suv_all_frequencies": np.array([
                           [[port.magnetic_suv for port in slab.right_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "source_ports_labels_all_frequencies": np.array([
                           [[port.label for port in slab.left_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "receiver_ports_labels_all_frequencies": np.array([
                           [[port.label for port in slab.right_ports.forward] for slab in current["slabs"]]
                           for _, current in setups]),
                       "single_channel_local_all_frequencies": np.array(local_values),
                       "single_channel_end_to_end_all_frequencies": np.array(full_values),
                       "phase_before_all_frequencies": np.array([r["phase_before_rad"] for r in runs]),
                       "phase_after_all_frequencies": np.array([r["phase_after_rad"] for r in runs]),
                       "tau_before_all_frequencies": np.array([r["tau_before_s"] for r in runs]),
                       "tau_after_all_frequencies": np.array([r["tau_after_s"] for r in runs])})
        failed = [name for name, item in metrics.items() if not item["pass"]]
        status = "PASS" if not failed else "FAIL"
        write_result(metrics, summary, arrays, status, None if not failed else ", ".join(failed))
        print(f"S2 {status}: {len(metrics) - len(failed)}/{len(metrics)} checks")
        for name in failed:
            print(f"FAIL {name}: {metrics[name]}")
        return 0 if not failed else 1
    except Exception as exc:
        write_result(metrics, summary, arrays, "STOPPED", str(exc))
        print(f"S2 STOPPED: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
