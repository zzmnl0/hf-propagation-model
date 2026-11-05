"""Controlled magnetized periodic-scatter case for the standalone route."""

from __future__ import annotations

import numpy as np
from scipy.constants import c

from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from models.ne_provider_geocentric import GeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF
from .config import ControlledCase, ONEWAY_AH_FD_SETTINGS
from .local_wave import MediumSnapshot
from .periodic import apply_periodic_window, solve_periodic_slab
from .propagation import make_field_template, prepare_windows, propagate_oneway
from .observables import doppler_hz, phase_delay_s, power_delay, power_doppler


class _UniformProvider(GeocentricNeProvider3D):
    def __init__(self, radius_km, density_m3):
        self.earth_radius_km = radius_km
        self.density_m3 = density_m3

    def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
        return self.density_m3

    def get_metadata(self):
        return {**super().get_metadata(), "earth_radius_km": self.earth_radius_km}


def _context(case, frequency_hz, relative_time_s):
    provider = _UniformProvider(case.earth_radius_km, case.background_Ne_m3)
    magnetic = ConstantGeomagneticField3D(case.magnetic_field_T * np.ones(3) / np.sqrt(3))
    models = {branch: RefractiveIndexAH3DECEF(
        provider, earth_radius_km=case.earth_radius_km, freq_MHz=frequency_hz / 1e6,
        wave_mode=branch, geomagnetic_provider=magnetic, plasma_coefficient=case.plasma_coefficient,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
    ) for branch in ("O", "X")}
    entrance = np.array([case.earth_radius_km + case.slab_height_km, 0., 0.])
    time_s = case.medium_time_s + relative_time_s
    medium = MediumSnapshot.from_ah_model(models["O"], entrance, time_s)
    half = case.grid_size // 2
    indices = np.array([(u, v) for u in range(-half, half + 1) for v in range(-half, half + 1)])
    periods = np.array([case.period_u_m, case.period_v_m])
    thickness = 2 * c / case.frequency_hz
    window = prepare_windows(models["O"], medium, entrance, np.eye(3), thickness,
                             periods, indices, time_s)
    return models, medium, window


def _density(case, relative_time_s):
    u = np.arange(case.grid_size) * case.period_u_m / case.grid_size
    v = np.arange(case.grid_size) * case.period_v_m / case.grid_size
    velocity_u = case.period_u_m / 8.0
    velocity_v = case.period_v_m / 10.0
    u_phase = 2 * np.pi * (u[:, None] - velocity_u * relative_time_s) / case.period_u_m
    v_phase = 2 * np.pi * (v[None, :] - velocity_v * relative_time_s) / case.period_v_m
    return case.background_Ne_m3 * (1 + case.irregularity_fraction *
                                   (np.cos(u_phase) + .6 * np.cos(v_phase)))


def _templates(case, window):
    indices = window["indices"]
    support = np.all(np.abs(indices) <= 1, axis=1)
    spatial = (np.exp(-.45 * np.sum(indices**2, axis=1))
               * np.exp(1j * (.23 * indices[:, 0] - .17 * indices[:, 1])) * support)
    source = np.column_stack((spatial, .45j * spatial))
    source /= np.linalg.norm(source)
    k_s = np.array([[port.k_suv_rad_m[0] for port in slab.right_ports.forward]
                    for slab in window["slabs"]])
    receiver = source * np.exp(1j * k_s * (case.before_m + 2 * c / case.frequency_hz + case.after_m))
    return make_field_template(window, source, "source"), make_field_template(window, receiver, "receiver")


def _solve(case, context, relative_time_s, source_template, receiver_template):
    models, medium, window = context
    ah = propagate_oneway(window, models, source_template, receiver_template,
                          case.before_m, case.after_m, medium.time_s)
    density = _density(case, relative_time_s)
    periodic = solve_periodic_slab(window, medium, density)
    result = apply_periodic_window(periodic, window, ah)
    indices = window["indices"]
    center = int(np.flatnonzero(np.all(indices == [0, 0], axis=1))[0])
    u_side = int(np.flatnonzero(np.all(indices == [1, 0], axis=1))[0])
    v_side = int(np.flatnonzero(np.all(indices == [0, 1], axis=1))[0])
    center_o = next(i for i, port in enumerate(window["slabs"][center].left_ports.forward)
                    if port.label == "O")
    side_u_o = next(i for i, port in enumerate(window["slabs"][u_side].right_ports.forward)
                    if port.label == "O")
    side_v_o = next(i for i, port in enumerate(window["slabs"][v_side].right_ports.forward)
                    if port.label == "O")
    col = 2 * center + center_o
    row_u = 2 * periodic.count + 2 * u_side + side_u_o
    row_v = 2 * periodic.count + 2 * v_side + side_v_o
    delta = periodic.S - periodic.S0
    return {
        "H": result["H_rx"], "H0": result["H_rx_background"],
        "S": periodic.S, "S0": periodic.S0, "delta_u": delta[row_u, col], "delta_v": delta[row_v, col],
        "density": density, "background_ne": medium.Ne_m3,
        "source": ah["source"], "receiver": ah["receiver"],
        "source_projection_norm": ah["source_projection_norm"],
        "receiver_projection_norm": ah["receiver_projection_norm"],
        "phase_before": ah["phase_before_rad"], "phase_after": ah["phase_after_rad"],
        "tau_before": ah["tau_before_s"], "tau_after": ah["tau_after_s"],
        "reflected_power": result["reflected_power_W"],
        "transmitted_power": result["transmitted_power_W"],
        "absorbed_power": result["absorbed_power_W"],
        "port_electric": np.array([[p.electric_suv for p in slab.right_ports.forward]
                                   for slab in window["slabs"]]),
        "port_magnetic": np.array([[p.magnetic_suv for p in slab.right_ports.forward]
                                   for slab in window["slabs"]]),
        "source_port_electric": np.array([[p.electric_suv for p in slab.left_ports.forward]
                                          for slab in window["slabs"]]),
        "source_port_magnetic": np.array([[p.magnetic_suv for p in slab.left_ports.forward]
                                          for slab in window["slabs"]]),
        "port_labels": np.array([[p.label for p in slab.right_ports.forward] for slab in window["slabs"]]),
        "port_k_ecef_rad_km": np.array([[p.k_ecef_rad_km for p in slab.right_ports.forward]
                                        for slab in window["slabs"]]),
        "launch_k_max_relative_error": ah["launch_k_max_relative_error"],
        "plane_error_max_km": ah["target_plane_max_error_km"],
    }


def run_controlled_case(case=None):
    case = ControlledCase() if case is None else case
    if case.grid_size != 5 or case.sweep_count != 9 or case.slow_time_count != 9:
        raise ValueError("the controlled reference uses a 5x5 grid and nine samples per axis")
    center_context = _context(case, case.frequency_hz, 0.)
    source_template, receiver_template = _templates(case, center_context[2])
    frequency_grid = np.linspace(case.frequency_hz - case.sweep_half_span_hz,
                                 case.frequency_hz + case.sweep_half_span_hz, case.sweep_count)
    probe_frequencies = case.frequency_hz + np.array([-1000., 1000., -500., 500.])
    all_frequencies = np.r_[frequency_grid, probe_frequencies]
    frequency_results = []
    center_result = None
    for frequency in all_frequencies:
        print(f"ray-wave frequency {frequency / 1e6:.6f} MHz", flush=True)
        context = center_context if frequency == case.frequency_hz else _context(case, float(frequency), 0.)
        result = _solve(case, context, 0., source_template, receiver_template)
        frequency_results.append(result)
        if frequency == case.frequency_hz:
            center_result = result
    slow_time = (np.arange(case.slow_time_count) - case.slow_time_count // 2) * case.slow_time_step_s
    slow_results = []
    for time_offset in slow_time:
        print(f"ray-wave slow time {time_offset:+.2f} s", flush=True)
        if time_offset == 0:
            result = center_result
        else:
            context = _context(case, case.frequency_hz, float(time_offset))
            result = _solve(case, context, float(time_offset), source_template, receiver_template)
        slow_results.append(result)

    def collect(results, key):
        return np.array([item[key] for item in results])

    freq_H = collect(frequency_results, "H")
    freq_H0 = collect(frequency_results, "H0")
    slow_H = collect(slow_results, "H")
    slow_H0 = collect(slow_results, "H0")
    source_norm = collect(frequency_results, "source_projection_norm")
    receiver_norm = collect(frequency_results, "receiver_projection_norm")
    delay, pd = power_delay(frequency_grid, freq_H[:case.sweep_count])
    _, pd0 = power_delay(frequency_grid, freq_H0[:case.sweep_count])
    _, pd_scattered = power_delay(frequency_grid, (freq_H - freq_H0)[:case.sweep_count])
    doppler_axis, doppler_power = power_doppler(slow_time, slow_H)
    _, doppler_scatter_power = power_doppler(slow_time, slow_H - slow_H0)
    result = {
        "frequency_grid_hz": frequency_grid, "all_frequencies_hz": all_frequencies,
        "frequency_H": freq_H, "frequency_H0": freq_H0,
        "frequency_source_projection_norm": source_norm,
        "frequency_receiver_projection_norm": receiver_norm,
        "frequency_source_scale": 1 / source_norm,
        "frequency_receiver_scale": 1 / receiver_norm,
        "frequency_Pin_physical_template_W": source_norm**2,
        "frequency_Prx_physical_source_W": np.abs(freq_H)**2 * source_norm**2,
        "frequency_Prx_per_Pin": np.abs(freq_H)**2,
        "frequency_H_physical_template_overlap": freq_H * source_norm * receiver_norm,
        "delay_axis_s": delay, "power_delay_W_per_W": pd,
        "power_delay_background_W_per_W": pd0, "power_delay_scatter_W_per_W": pd_scattered,
        "slow_time_s": slow_time, "slow_H": slow_H, "slow_H0": slow_H0,
        "slow_doppler_hz": doppler_hz(slow_H, slow_time),
        "slow_scatter_doppler_hz": doppler_hz(slow_H - slow_H0, slow_time),
        "doppler_axis_hz": doppler_axis, "power_doppler_W_per_W": doppler_power,
        "power_doppler_scatter_W_per_W": doppler_scatter_power,
        "slow_delta_u": collect(slow_results, "delta_u"), "slow_delta_v": collect(slow_results, "delta_v"),
        "source_template_E": source_template["electric"], "source_template_H": source_template["magnetic"],
        "receiver_template_E": receiver_template["electric"], "receiver_template_H": receiver_template["magnetic"],
        "indices": center_context[2]["indices"],
    }
    for name in ("S", "S0", "source", "receiver", "phase_before", "phase_after", "tau_before", "tau_after",
                 "background_ne", "port_electric", "port_magnetic", "source_port_electric", "source_port_magnetic",
                 "port_labels", "port_k_ecef_rad_km",
                 "reflected_power", "transmitted_power", "absorbed_power", "launch_k_max_relative_error",
                 "plane_error_max_km"):
        result["frequency_" + name] = collect(frequency_results, name)
    for name in ("S", "S0", "source", "receiver", "phase_before", "phase_after", "density",
                 "background_ne", "reflected_power", "transmitted_power", "absorbed_power"):
        result["slow_" + name] = collect(slow_results, name)
    result["group_delay_1khz_s"] = np.array(phase_delay_s(freq_H[10], freq_H[9], 1000.))
    result["group_delay_500hz_s"] = np.array(phase_delay_s(freq_H[12], freq_H[11], 500.))
    return result
