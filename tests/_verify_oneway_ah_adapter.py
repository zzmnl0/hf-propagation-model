"""Independent controls for the separated AH-to-plane geometry adapter."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config as original_config
from models.dispersion_ah import C_KM_S
from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from models.ne_provider_geocentric import GeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    group_velocity_ah3d_ecef,
    integrate_fixed_path_ah3d_ecef,
    launch_state_ah3d_ecef,
)
from ray_wave_oneway.ah_adapter import trace_to_plane
from ray_wave_oneway.config import ONEWAY_AH_FD_SETTINGS


OUT = Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "ah_adapter"
RADIUS = 6371.0
HEIGHT = 200.0
FREQ_HZ = 10e6
COEFFICIENT = 80.61638587963631
X_START = RADIUS + HEIGHT
N2_START = 0.8
SLOPE = 0.02
TIME_S = 37.0


class FlatProvider(GeocentricNeProvider3D):
    earth_radius_km = RADIUS

    def __init__(self, slope, temporal_rate=0.):
        self.slope = float(slope)
        self.temporal_rate = float(temporal_rate)

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        del earth_radius_km
        n2 = (N2_START - self.slope * (float(ecef_km[0]) - X_START)
              - self.temporal_rate * (float(t) - TIME_S))
        return (1 - n2) * FREQ_HZ**2 / COEFFICIENT

    def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
        del lat_deg, lon_deg, height_km, t
        raise NotImplementedError("this control is expressed in ECEF")

    def get_metadata(self):
        return {**super().get_metadata(), "earth_radius_km": self.earth_radius_km,
                "n2_start": N2_START, "slope_per_km": self.slope,
                "temporal_rate_per_s": self.temporal_rate}


class BoundedProvider(FlatProvider):
    def __init__(self):
        super().__init__(0.)

    def domain_margin_ecef(self, ecef_km, earth_radius_km):
        del earth_radius_km
        return float(X_START + 0.16 - ecef_km[0])

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        if self.domain_margin_ecef(ecef_km, earth_radius_km) < 0:
            raise ValueError("outside valid density domain")
        return super().get_Ne_ecef(ecef_km, earth_radius_km, t=t)


def make_model(slope, field, temporal_rate=0.):
    return RefractiveIndexAH3DECEF(
        FlatProvider(slope, temporal_rate), earth_radius_km=RADIUS, freq_MHz=FREQ_HZ / 1e6,
        wave_mode="O", geomagnetic_provider=ConstantGeomagneticField3D(field),
        plasma_coefficient=COEFFICIENT,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
    )


def launch(model, elevation, time_s=TIME_S):
    return launch_state_ah3d_ecef(model, tx_lat_deg=0., tx_lon_deg=0.,
                                  tx_height_km=HEIGHT, elevation_deg=elevation,
                                  bearing_deg=90., medium_time_s=time_s)


def check(metrics, name, actual, expected, tolerance):
    error = float(abs(actual - expected))
    metrics[name] = {"actual": float(actual), "expected": float(expected),
                     "absolute_error": error, "tolerance": float(tolerance),
                     "pass": bool(np.isfinite(error) and error <= tolerance)}


def main():
    metrics = {}
    uniform = make_model(0., np.full(3, 45e-6 / np.sqrt(3)))
    initial = launch(uniform, 45.)
    velocity = group_velocity_ah3d_ecef(initial, uniform, t=TIME_S)
    unit_velocity = velocity / np.linalg.norm(velocity)
    distance = 0.4
    plane = initial[:3] + distance * unit_velocity
    result = trace_to_plane(uniform, initial[:3], initial[3:6], plane, unit_velocity,
                            max_path_km=0.6, direction=1, medium_time_s=TIME_S)
    assert result["termination_status"] == "receiver_plane_crossed"
    expected_action = float(np.dot(initial[3:6], unit_velocity) * distance)
    check(metrics, "uniform_endpoint_error_km", np.linalg.norm(result["endpoint_ecef_km"] - plane), 0., 2e-8)
    check(metrics, "uniform_wavevector_error_rad_km", np.linalg.norm(result["endpoint_k_ecef_rad_km"] - initial[3:6]), 0., 2e-8)
    check(metrics, "uniform_action_rad", result["action_rad"][-1], expected_action, 2e-6)
    check(metrics, "uniform_group_time_s", result["group_time_s"][-1], distance / np.linalg.norm(velocity), 2e-12)

    old = integrate_fixed_path_ah3d_ecef(
        uniform, tx_lat_deg=0., tx_lon_deg=0., tx_height_km=HEIGHT,
        elevation_deg=45., bearing_deg=90., path_length_km=distance,
        medium_time_s=TIME_S, profile="L2", stop_on_physical_events=False,
    )
    assert old["termination_status"] == "fixed_geometric_path_completed"
    check(metrics, "original_ah_endpoint_error_km", np.linalg.norm(result["endpoint_ecef_km"] - old["trajectory_ecef_km"][-1]), 0., 2e-8)
    check(metrics, "original_ah_wavevector_error_rad_km", np.linalg.norm(result["endpoint_k_ecef_rad_km"] - old["wavevector_ecef_km^-1"][-1]), 0., 2e-8)
    check(metrics, "original_ah_group_time_s", result["group_time_s"][-1], old["trajectory_group_time_s"][-1], 2e-12)

    old_settings = original_config.UPDATE2["AH3D_INTEGRATION"]["L2"]
    saved_settings = old_settings.copy()
    saved_fd = original_config.UPDATE2["AH3D_SPATIAL_FD_BASE_KM"]
    had_profile = "AH3D_INTEGRATION_PROFILE" in original_config.UPDATE2
    saved_profile = original_config.UPDATE2.get("AH3D_INTEGRATION_PROFILE")
    try:
        old_settings.update({"rtol": .1, "atol_r_km": .1, "atol_p": .1,
                             "ds_min_km": .2, "ds_max_km": .2})
        original_config.UPDATE2["AH3D_SPATIAL_FD_BASE_KM"] = 40.
        original_config.UPDATE2["AH3D_INTEGRATION_PROFILE"] = "research"
        isolated_model = make_model(0., np.full(3, 45e-6 / np.sqrt(3)))
        isolated_initial = launch(isolated_model, 45.)
        isolated = trace_to_plane(isolated_model, isolated_initial[:3], isolated_initial[3:6],
                                  plane, unit_velocity, max_path_km=.6, direction=1,
                                  medium_time_s=TIME_S, profile="L2")
    finally:
        old_settings.clear()
        old_settings.update(saved_settings)
        original_config.UPDATE2["AH3D_SPATIAL_FD_BASE_KM"] = saved_fd
        if had_profile:
            original_config.UPDATE2["AH3D_INTEGRATION_PROFILE"] = saved_profile
        else:
            del original_config.UPDATE2["AH3D_INTEGRATION_PROFILE"]
    assert all(getattr(isolated_model, key) == value for key, value in ONEWAY_AH_FD_SETTINGS.items())
    assert all(np.array_equal(isolated[key], result[key]) for key in (
        "path_km", "trajectory_ecef_km", "wavevector_ecef_rad_km", "group_time_s", "action_rad"))
    metrics["original_update2_settings_isolation"] = {"pass": True}
    try:
        trace_to_plane(uniform, initial[:3], initial[3:6], plane, unit_velocity,
                       max_path_km=.6, profile="research", medium_time_s=TIME_S)
    except ValueError as exc:
        assert "unknown one-way integration profile" in str(exc)
    else:
        raise AssertionError("unsupported integration profile was accepted")
    metrics["reject_unsupported_profile"] = {"pass": True}

    temporal_model = make_model(0., np.zeros(3), temporal_rate=0.001)
    temporal_initial = launch(temporal_model, 45., time_s=TIME_S + 5.)
    temporal_velocity = group_velocity_ah3d_ecef(temporal_initial, temporal_model, t=TIME_S + 5.)
    temporal_direction = temporal_velocity / np.linalg.norm(temporal_velocity)
    temporal_plane = temporal_initial[:3] + distance * temporal_direction
    temporal = trace_to_plane(temporal_model, temporal_initial[:3], temporal_initial[3:6],
                              temporal_plane, temporal_direction, max_path_km=.6,
                              direction=1, medium_time_s=TIME_S + 5.)
    assert temporal["termination_status"] == "receiver_plane_crossed"
    expected_n = np.sqrt(N2_START - 0.005)
    check(metrics, "frozen_medium_time_group_s", temporal["group_time_s"][-1],
          distance / (C_KM_S * expected_n), 2e-12)
    check(metrics, "frozen_medium_time_action_rad", temporal["action_rad"][-1],
          temporal_model.k0 * expected_n * distance, 2e-6)

    model = make_model(SLOPE, np.zeros(3))
    initial = launch(model, 45.)
    normal = np.array([1., 0., 0.])
    reflected = trace_to_plane(model, initial[:3], initial[3:6], initial[:3], normal,
                               max_path_km=200., direction=-1, medium_time_s=TIME_S,
                               max_step_km=2.)
    assert reflected["termination_status"] == "receiver_plane_crossed"
    assert reflected["path_km"][-1] > 0.
    k0 = model.k0
    tangential = np.linalg.norm(initial[4:6]) / k0
    a = N2_START - tangential**2
    turning = a / SLOPE
    upper = np.sqrt(turning)
    factor = 4 / np.sqrt(SLOPE)
    expected_y = factor * quad(lambda t: tangential, 0., upper)[0]
    expected_time = factor * quad(lambda t: 1 / C_KM_S, 0., upper)[0]
    expected_phase = factor * k0 * quad(lambda t: tangential**2 + SLOPE * t**2, 0., upper)[0]
    check(metrics, "reflection_receiver_plane_km", reflected["terminal_plane_distance_km"], 0., 2e-8)
    check(metrics, "reflection_tangential_displacement_km", reflected["endpoint_ecef_km"][1], expected_y, 3e-3)
    check(metrics, "reflection_group_time_s", reflected["group_time_s"][-1], expected_time, 2e-8)
    check(metrics, "reflection_action_rad", reflected["action_rad"][-1], expected_phase, .15)
    check(metrics, "reflection_tangential_k_drift_rad_km", np.max(np.linalg.norm(reflected["wavevector_ecef_rad_km"][:, 1:] - initial[4:6], axis=1)), 0., 2e-4)
    check(metrics, "reflection_max_shell_residual", reflected["max_normalized_dispersion_residual"], 0., 1e-8)
    assert np.max(reflected["trajectory_ecef_km"][:, 0] - X_START) > turning - 0.1
    assert np.any(reflected["wavevector_ecef_rad_km"][:, 0] > 0.)
    assert reflected["endpoint_k_normal_rad_km"] < 0.
    metrics["reflection_turning_height_km"] = {
        "sampled_max": float(np.max(reflected["trajectory_ecef_km"][:, 0] - X_START)),
        "expected": float(turning), "pass": True,
    }

    old_reflection = integrate_fixed_path_ah3d_ecef(
        model, tx_lat_deg=0., tx_lon_deg=0., tx_height_km=HEIGHT,
        elevation_deg=45., bearing_deg=90., path_length_km=reflected["path_km"][-1],
        medium_time_s=TIME_S, profile="L2", stop_on_physical_events=False,
    )
    assert old_reflection["termination_status"] == "fixed_geometric_path_completed"
    check(metrics, "nonuniform_original_ah_endpoint_error_km",
          np.linalg.norm(reflected["endpoint_ecef_km"] - old_reflection["trajectory_ecef_km"][-1]),
          0., 3e-3)
    check(metrics, "nonuniform_original_ah_group_time_s", reflected["group_time_s"][-1],
          old_reflection["trajectory_group_time_s"][-1], 2e-8)
    refined = trace_to_plane(model, initial[:3], initial[3:6], initial[:3], normal,
                             max_path_km=200., direction=-1, medium_time_s=TIME_S,
                             max_step_km=1.)
    check(metrics, "reflection_step_refinement_endpoint_km",
          np.linalg.norm(reflected["endpoint_ecef_km"] - refined["endpoint_ecef_km"]), 0., 1e-4)
    check(metrics, "reflection_step_refinement_action_rad", reflected["action_rad"][-1],
          refined["action_rad"][-1], 0.01)

    missed = trace_to_plane(uniform, initial[:3], launch(uniform, 45.)[3:6],
                            initial[:3] - 1.0 * normal, normal, max_path_km=0.05,
                            direction=-1, medium_time_s=TIME_S)
    assert missed["termination_status"] == "max_path_reached"
    metrics["missed_plane_status"] = {"pass": True}

    try:
        trace_to_plane(uniform, initial[:3], launch(uniform, 45.)[3:6] * 1.01,
                       plane, unit_velocity, max_path_km=1., medium_time_s=TIME_S)
    except ValueError as exc:
        assert "dispersion shell" in str(exc)
    else:
        raise AssertionError("off-shell port was accepted")
    metrics["reject_off_shell_port"] = {"pass": True}

    bounded_model = RefractiveIndexAH3DECEF(
        BoundedProvider(), earth_radius_km=RADIUS, freq_MHz=FREQ_HZ / 1e6,
        wave_mode="O", geomagnetic_provider=ConstantGeomagneticField3D(np.zeros(3)),
        plasma_coefficient=COEFFICIENT,
        integration_profile="L2", **ONEWAY_AH_FD_SETTINGS,
    )
    bounded_initial = launch(bounded_model, 90.)
    bounded_plane = bounded_initial[:3] + 0.05 * normal
    bounded = trace_to_plane(bounded_model, bounded_initial[:3], bounded_initial[3:6],
                             bounded_plane, normal, max_path_km=.4, direction=1,
                             medium_time_s=TIME_S, max_step_km=.5)
    fine = trace_to_plane(bounded_model, bounded_initial[:3], bounded_initial[3:6],
                          bounded_plane, normal, max_path_km=.4, direction=1,
                          medium_time_s=TIME_S, max_step_km=.01)
    assert bounded["termination_status"] == "receiver_plane_crossed"
    assert bounded["provider_boundary_retries"] > 0
    assert fine["termination_status"] == "receiver_plane_crossed"
    metrics["bounded_provider_retry"] = {"pass": True,
                                         "retries": bounded["provider_boundary_retries"]}
    check(metrics, "bounded_provider_endpoint_error_km",
          np.linalg.norm(bounded["endpoint_ecef_km"] - fine["endpoint_ecef_km"]), 0., 1e-8)
    check(metrics, "bounded_provider_group_time_s", bounded["group_time_s"][-1],
          fine["group_time_s"][-1], 1e-12)
    check(metrics, "bounded_provider_action_rad", bounded["action_rad"][-1],
          fine["action_rad"][-1], 1e-7)

    try:
        trace_to_plane(bounded_model, bounded_initial[:3], bounded_initial[3:6],
                       bounded_initial[:3] + .12 * normal, normal, max_path_km=.4,
                       direction=1, medium_time_s=TIME_S, max_step_km=.5)
    except RuntimeError as exc:
        assert "provider derivative domain" in str(exc)
        assert exc.__cause__ is not None and "domain" in str(exc.__cause__)
    else:
        raise AssertionError("plane outside the derivative domain was accepted")
    metrics["bounded_provider_unreachable_plane"] = {"pass": True}

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "trajectories.npz", uniform_position=result["trajectory_ecef_km"],
                        uniform_k=result["wavevector_ecef_rad_km"],
                        reflected_position=reflected["trajectory_ecef_km"],
                        reflected_k=reflected["wavevector_ecef_rad_km"],
                        reflected_path=reflected["path_km"], reflected_time=reflected["group_time_s"],
                        reflected_action=reflected["action_rad"])
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    failed = [name for name, value in metrics.items() if not value["pass"]]
    for name, value in metrics.items():
        print(f"{'PASS' if value['pass'] else 'FAIL'} {name}: {value}")
    print(f"{len(metrics) - len(failed)}/{len(metrics)} checks passed")
    if failed:
        raise AssertionError(f"failed checks: {failed}")


if __name__ == "__main__":
    main()
