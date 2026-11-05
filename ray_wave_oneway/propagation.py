"""One-way AH / planar full-wave coupling for a uniform exterior medium.

The transverse spectrum is a finite periodic cell.  Port coefficients carry
watts, while reconstructed E/H fields include the cell-area normalization.
"""

from __future__ import annotations

import numpy as np

from .local_wave import MediumSnapshot, solve_slab
from .ah_adapter import trace_to_plane
from models.ray_tracer_ah_3d_ecef import (
    WGS84_COORDINATE_MODEL,
    group_velocity_ah3d_ecef,
    launch_state_ah3d_ecef,
)
from utils import enu_basis, enu_basis_wgs84


def _same_medium(reference, sample):
    same_settings = (
        reference.frequency_hz == sample.frequency_hz
        and reference.time_s == sample.time_s
        and reference.plasma_coefficient == sample.plasma_coefficient
        and reference.earth_radius_km == sample.earth_radius_km
        and reference.coordinate_model == sample.coordinate_model
        and reference.near_cutoff_n2 == sample.near_cutoff_n2
        and reference.collision_s == sample.collision_s
    )
    density_close = np.isclose(reference.Ne_m3, sample.Ne_m3, rtol=1e-10, atol=1e-6)
    field_scale = max(float(np.linalg.norm(reference.B_ecef_T)), 1e-18)
    field_close = np.linalg.norm(reference.B_ecef_T - sample.B_ecef_T) <= 1e-10 * field_scale
    if not same_settings or not density_close or not field_close:
        raise ValueError("S2 requires a common, uniform, collisionless exterior medium")


def prepare_windows(model, layer, entrance_ecef_km, basis_suv_ecef,
                    thickness_m, periods_m, indices, medium_time_s):
    """Construct S/S0 for the fixed discrete transverse spectrum."""
    basis = np.asarray(basis_suv_ecef, dtype=float)
    entrance = np.asarray(entrance_ecef_km, dtype=float)
    periods = np.asarray(periods_m, dtype=float)
    bins = np.asarray(indices)
    if (entrance.shape != (3,) or not np.all(np.isfinite(entrance))
            or periods.shape != (2,) or not np.all(np.isfinite(periods))
            or np.any(periods <= 0) or bins.ndim != 2 or bins.shape[1] != 2
            or len(bins) == 0
            or not np.all(np.isfinite(bins)) or not np.all(bins == np.rint(bins))
            or len(np.unique(bins, axis=0)) != len(bins)):
        raise ValueError("invalid entrance, periods or transverse indices")
    bins = bins.astype(int)
    if not np.isfinite(thickness_m) or thickness_m <= 0:
        raise ValueError("thickness_m must be finite and positive")
    if (basis.shape != (3, 3) or not np.all(np.isfinite(basis))
            or not np.allclose(basis.T @ basis, np.eye(3), rtol=0, atol=1e-10)
            or abs(np.linalg.det(basis) - 1) > 1e-10):
        raise ValueError("basis must be right-handed and orthonormal")
    exit_position = entrance + basis[:, 0] * thickness_m / 1000.0
    midpoint = (entrance + exit_position) / 2
    left = MediumSnapshot.from_ah_model(model, entrance, medium_time_s)
    right = MediumSnapshot.from_ah_model(model, exit_position, medium_time_s)
    background = MediumSnapshot.from_ah_model(model, midpoint, medium_time_s)
    if left.collision_s != 0:
        raise ValueError("S2 exterior must be lossless")
    _same_medium(left, right)
    _same_medium(left, background)
    tangent_k = 2 * np.pi * bins / periods
    slabs = [solve_slab(left, layer, right, background, basis, ku, kv, thickness_m)
             for ku, kv in tangent_k]
    return {"basis": basis, "entrance_ecef_km": entrance,
            "exit_ecef_km": exit_position, "periods_m": periods,
            "indices": bins, "tangent_k_rad_m": tangent_k,
            "left_medium": left, "right_medium": right, "slabs": slabs}


def _ports(window, side):
    if side == "source":
        return [slab.left_ports.forward for slab in window["slabs"]]
    if side == "receiver":
        return [slab.right_ports.forward for slab in window["slabs"]]
    raise ValueError("side must be 'source' or 'receiver'")


def make_field_template(window, coefficients, side):
    """Freeze a physical E/H template at this frequency and reference plane."""
    values = np.asarray(coefficients, dtype=complex)
    ports = _ports(window, side)
    if values.shape != (len(ports), 2) or not np.all(np.isfinite(values)):
        raise ValueError("template coefficients must have shape (bins, 2)")
    norm = float(np.linalg.norm(values))
    if norm <= 0:
        raise ValueError("template has zero power")
    values = values / norm
    electric = np.array([sum(values[i, j] * port.electric_suv
                             for j, port in enumerate(group))
                         for i, group in enumerate(ports)])
    magnetic = np.array([sum(values[i, j] * port.magnetic_suv
                             for j, port in enumerate(group))
                         for i, group in enumerate(ports)])
    return {"indices": window["indices"].copy(), "periods_m": window["periods_m"].copy(),
            "basis": window["basis"].copy(), "electric": electric, "magnetic": magnetic}


def _power_inner(electric_a, magnetic_a, electric_b, magnetic_b):
    return 0.25 * (np.cross(np.conj(electric_a), magnetic_b)
                   + np.cross(electric_b, np.conj(magnetic_a)))[0]


def project_field_template(window, template, side):
    """Project the fixed physical template to power-normalized current ports."""
    if (not np.array_equal(window["indices"], template["indices"])
            or not np.array_equal(window["periods_m"], template["periods_m"])
            or not np.allclose(window["basis"], template["basis"], rtol=0, atol=1e-12)):
        raise ValueError("template and window spatial bases differ")
    ports = _ports(window, side)
    if (np.shape(template["electric"]) != (len(ports), 3)
            or np.shape(template["magnetic"]) != (len(ports), 3)
            or not np.all(np.isfinite(template["electric"]))
            or not np.all(np.isfinite(template["magnetic"]))):
        raise ValueError("invalid physical E/H template")
    projected = np.array([[_power_inner(port.electric_suv, port.magnetic_suv,
                                        template["electric"][i], template["magnetic"][i])
                           for port in group] for i, group in enumerate(ports)])
    raw_norm = float(np.linalg.norm(projected))
    if not np.isfinite(raw_norm) or raw_norm <= 1e-15:
        raise ValueError("template has zero overlap with propagating ports")
    return projected / raw_norm, raw_norm


def reconstruct_fields(window, coefficients, side, grid_size):
    """Reconstruct the periodic E/H field and integrated normal power."""
    ports = _ports(window, side)
    values = np.asarray(coefficients, dtype=complex)
    if (values.shape != (len(ports), 2) or not np.all(np.isfinite(values))
            or not isinstance(grid_size, int) or grid_size < 2):
        raise ValueError("invalid field coefficients or grid size")
    if grid_size <= int(np.max(np.ptp(window["indices"], axis=0))):
        raise ValueError("grid size aliases the selected transverse spectrum")
    periods = window["periods_m"]
    u = np.arange(grid_size) * periods[0] / grid_size
    v = np.arange(grid_size) * periods[1] / grid_size
    uu, vv = np.meshgrid(u, v, indexing="ij")
    angle = (window["tangent_k_rad_m"][:, 0, None, None] * uu
             + window["tangent_k_rad_m"][:, 1, None, None] * vv)
    harmonics = np.exp(1j * angle)
    electric_coeff = np.array([sum(values[i, j] * port.electric_suv
                                   for j, port in enumerate(group))
                               for i, group in enumerate(ports)])
    magnetic_coeff = np.array([sum(values[i, j] * port.magnetic_suv
                                   for j, port in enumerate(group))
                               for i, group in enumerate(ports)])
    area = float(np.prod(periods))
    electric = np.einsum("nuv,nd->uvd", harmonics, electric_coeff) / np.sqrt(area)
    magnetic = np.einsum("nuv,nd->uvd", harmonics, magnetic_coeff) / np.sqrt(area)
    flux = 0.5 * np.real(np.cross(electric, np.conj(magnetic))[..., 0])
    return electric, magnetic, float(area * flux.mean())


def _trace_to_plane(model, port, start_ecef_km, target_ecef_km, basis, time_s, reference):
    """Reconnect one outward port to a parallel plane through the actual AH solver."""
    _same_medium(reference, MediumSnapshot.from_ah_model(model, start_ecef_km, time_s))
    _same_medium(reference, MediumSnapshot.from_ah_model(model, target_ecef_km, time_s))
    lat, lon, height = model.coordinates_from_ecef(start_ecef_km)
    enu = (enu_basis_wgs84(lat, lon) if model.coordinate_model == WGS84_COORDINATE_MODEL
           else enu_basis(lat, lon, model.earth_radius_km))
    direction = port.k_ecef_rad_km / np.linalg.norm(port.k_ecef_rad_km)
    east, north, up = enu.T @ direction
    elevation = float(np.degrees(np.arctan2(up, np.hypot(east, north))))
    bearing = float(np.degrees(np.arctan2(east, north)) % 360)
    launch = launch_state_ah3d_ecef(
        model, tx_lat_deg=float(lat), tx_lon_deg=float(lon), tx_height_km=float(height),
        elevation_deg=elevation, bearing_deg=bearing, medium_time_s=time_s,
    )
    k_error = float(np.linalg.norm(launch[3:6] - port.k_ecef_rad_km)
                    / np.linalg.norm(port.k_ecef_rad_km))
    if k_error > 1e-8:
        raise ValueError("AH launch wavevector does not match the local port")
    velocity = group_velocity_ah3d_ecef(launch, model, t=time_s)
    normal_speed = float(np.dot(velocity, basis[:, 0]) / np.linalg.norm(velocity))
    normal_distance_km = float(np.dot(target_ecef_km - start_ecef_km, basis[:, 0]))
    if normal_speed <= 1e-8 or normal_distance_km <= 0:
        raise ValueError("AH energy flow must cross the target plane outward")
    length_km = normal_distance_km / normal_speed
    ray = trace_to_plane(
        model, start_ecef_km, launch[3:6], target_ecef_km, basis[:, 0],
        max_path_km=max(1.1 * length_km, length_km + 0.01), direction=1,
        medium_time_s=time_s, profile="L2",
    )
    if ray["termination_status"] != "receiver_plane_crossed":
        raise ValueError(f"AH ray failed: {ray['termination_status']}")
    points = ray["trajectory_ecef_km"]
    waves = ray["wavevector_ecef_rad_km"]
    for point in points:
        _same_medium(reference, MediumSnapshot.from_ah_model(model, point, time_s))
    terminal_error_km = float(np.dot(points[-1] - target_ecef_km, basis[:, 0]))
    phase_ray = float(ray["action_rad"][-1])
    plane_phase = phase_ray - float(np.dot(waves[-1], points[-1] - target_ecef_km))
    uniform_phase = float(port.k_suv_rad_m[0] * normal_distance_km * 1000)
    return {"phase_rad": plane_phase, "ray_phase_rad": phase_ray,
            "uniform_phase_rad": uniform_phase, "group_time_s": float(ray["group_time_s"][-1]),
            "endpoint_ecef_km": points[-1], "endpoint_k_ecef_rad_km": waves[-1],
            "terminal_plane_error_km": terminal_error_km, "launch_k_relative_error": k_error,
            "path_length_km": float(ray["path_km"][-1]), "group_direction_cosine": normal_speed,
            "accepted_steps": int(len(ray["path_km"]) - 1)}


def propagate_oneway(window, models_by_mode, source_template, receiver_template,
                     before_distance_m, after_distance_m, medium_time_s):
    """Run AH exterior segments, replace the window with S, and project RX."""
    if set(models_by_mode) != {"O", "X"}:
        raise ValueError("models_by_mode must contain O and X AH models")
    if any(model.wave_mode != branch for branch, model in models_by_mode.items()):
        raise ValueError("AH model keys and wave branches must agree")
    if (not np.isfinite(before_distance_m) or before_distance_m <= 0
            or not np.isfinite(after_distance_m) or after_distance_m <= 0):
        raise ValueError("exterior distances must be finite and positive")
    reference = window["left_medium"]
    if medium_time_s != reference.time_s:
        raise ValueError("window and propagation medium times differ")
    basis = window["basis"]
    source_position = window["entrance_ecef_km"] - basis[:, 0] * before_distance_m / 1000
    receiver_position = window["exit_ecef_km"] + basis[:, 0] * after_distance_m / 1000
    source, source_norm = project_field_template(window, source_template, "source")
    receiver, receiver_norm = project_field_template(window, receiver_template, "receiver")
    count = len(window["slabs"])
    before_phase = np.empty((count, 2))
    after_phase = np.empty((count, 2))
    before_tau = np.empty((count, 2))
    after_tau = np.empty((count, 2))
    launch_errors = []
    plane_errors = []
    phase_errors = []
    rays = []
    transmitted = np.empty((count, 2), dtype=complex)
    transmitted_background = np.empty_like(transmitted)
    reflected = np.empty_like(transmitted)
    incident_at_window = np.empty_like(transmitted)
    absorbed = 0.0
    for i, slab in enumerate(window["slabs"]):
        for j, port in enumerate(slab.left_ports.forward):
            if port.label not in models_by_mode:
                raise ValueError("S2 requires nondegenerate magnetized O/X ports")
            before = _trace_to_plane(models_by_mode[port.label], port, source_position,
                                     window["entrance_ecef_km"], basis, medium_time_s, reference)
            before_phase[i, j], before_tau[i, j] = before["phase_rad"], before["group_time_s"]
            rays.append(before)
        for j, port in enumerate(slab.right_ports.forward):
            if port.label not in models_by_mode:
                raise ValueError("S2 requires nondegenerate magnetized O/X ports")
            after = _trace_to_plane(models_by_mode[port.label], port, window["exit_ecef_km"],
                                    receiver_position, basis, medium_time_s, reference)
            after_phase[i, j], after_tau[i, j] = after["phase_rad"], after["group_time_s"]
            rays.append(after)
        at_window = source[i] * np.exp(1j * before_phase[i])
        incident_at_window[i] = at_window
        reflected[i] = slab.S[:2, :2] @ at_window
        transmitted[i] = (slab.S[2:, :2] @ at_window) * np.exp(1j * after_phase[i])
        transmitted_background[i] = (slab.S0[2:, :2] @ at_window) * np.exp(1j * after_phase[i])
        absorbed += slab.absorbed_power_W(np.r_[at_window, [0j, 0j]], quadrature_order=32)
    for ray in rays:
        launch_errors.append(ray["launch_k_relative_error"])
        plane_errors.append(abs(ray["terminal_plane_error_km"]))
        phase_errors.append(abs(ray["phase_rad"] - ray["uniform_phase_rad"]))
    return {
        "window_context": window,
        "source": source, "receiver": receiver,
        "source_projection_norm": source_norm, "receiver_projection_norm": receiver_norm,
        "incident_at_window": incident_at_window, "reflected": reflected,
        "transmitted": transmitted, "transmitted_background": transmitted_background,
        "H_rx": np.vdot(receiver, transmitted), "H_rx_background": np.vdot(receiver, transmitted_background),
        "incident_power_W": float(np.vdot(source, source).real),
        "reflected_power_W": float(np.vdot(reflected, reflected).real),
        "transmitted_power_W": float(np.vdot(transmitted, transmitted).real),
        "absorbed_power_W": float(absorbed),
        "rx_power_W": float(abs(np.vdot(receiver, transmitted))**2),
        "uncaptured_power_W": float(np.vdot(transmitted, transmitted).real - abs(np.vdot(receiver, transmitted))**2),
        "phase_before_rad": before_phase, "phase_after_rad": after_phase,
        "tau_before_s": before_tau, "tau_after_s": after_tau,
        "launch_k_max_relative_error": max(launch_errors),
        "target_plane_max_error_km": max(plane_errors),
        "uniform_phase_max_error_rad": max(phase_errors),
        "source_position_ecef_km": source_position, "receiver_position_ecef_km": receiver_position,
        "rays": rays,
    }
