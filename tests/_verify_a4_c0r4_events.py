"""Acceptance checks for C0R-4 spherical physical events and output semantics."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.ne_provider_geocentric import (
    ChapmanGeocentricNeProvider3D,
    GeocentricGridNeProvider3D,
    VacuumGeocentricNeProvider3D,
)
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)


R_KM = 6371.0
REPORT = ROOT / "output2" / "a4_pharlap" / "c0r4_event_validation.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _model(provider):
    return RefractiveIndexAH3DECEF(
        provider, earth_radius_km=R_KM, freq_MHz=10.0,
        wave_mode="O", gyro_MHz=0.0, b_ecef_hat=(0.0, 0.0, 1.0),
        integration_profile="smoke")


def _trace(model, *, height=0.0, elevation=35.0, h_stop=600.0,
           max_path=3000.0):
    return trace_ray_ah3d_ecef(
        model, tx_lat_deg=30.0, tx_lon_deg=120.0,
        tx_height_km=height, elevation_deg=elevation, bearing_deg=37.0,
        h_stop_km=h_stop, max_path_km=max_path, profile="smoke")


def run() -> int:
    vacuum_model = _model(VacuumGeocentricNeProvider3D(R_KM))

    upper_height = 100.0
    upper_elevation = 35.0
    upper = _trace(
        vacuum_model, elevation=upper_elevation,
        h_stop=upper_height, max_path=500.0)
    radial = np.sin(np.radians(upper_elevation))
    expected_upper_path = (
        -R_KM * radial
        + np.sqrt((R_KM * radial) ** 2
                  + 2.0 * R_KM * upper_height + upper_height**2)
    )

    tx_height = 10.0
    ground_elevation = -30.0
    ground = _trace(
        vacuum_model, height=tx_height, elevation=ground_elevation,
        h_stop=100.0, max_path=100.0)
    radius_tx = R_KM + tx_height
    radial = np.sin(np.radians(ground_elevation))
    expected_ground_path = (
        -radius_tx * radial
        - np.sqrt((radius_tx * radial) ** 2
                  - (radius_tx**2 - R_KM**2))
    )

    surface_downward_rejected = False
    try:
        _trace(vacuum_model, elevation=-1.0, h_stop=100.0, max_path=50.0)
    except ValueError as exc:
        surface_downward_rejected = "depart outward" in str(exc)

    latitude = np.array([29.0, 30.0, 31.0])
    longitude = np.array([119.0, 120.0, 121.0])
    height = np.array([0.0, 10.0, 20.0])
    Ne = np.zeros((latitude.size, longitude.size, height.size))
    grid_model = _model(GeocentricGridNeProvider3D(
        Ne, latitude, longitude, height, R_KM))
    grid_exit = _trace(
        grid_model, elevation=60.0, h_stop=100.0, max_path=100.0)

    chapman_model = _model(ChapmanGeocentricNeProvider3D(
        R_KM, foF2_MHz=7.0, hmF2_km=350.0, ymF2_km=100.0))
    reflected = _trace(
        chapman_model, elevation=35.0, h_stop=600.0, max_path=3000.0)

    upper_error = abs(upper["geom_path_km"] - expected_upper_path)
    ground_error = abs(ground["geom_path_km"] - expected_ground_path)
    grid_exit_height_error = abs(
        grid_exit["trajectory_physical_height_km"][-1] - 20.0)
    metrics = {
        "vacuum_upper_path_error_km": upper_error,
        "vacuum_upper_height_error_km": abs(
            upper["trajectory_physical_height_km"][-1] - upper_height),
        "vacuum_ground_path_error_km": ground_error,
        "vacuum_ground_height_error_km": abs(ground["landing_height_km"]),
        "grid_exit_height_error_km": grid_exit_height_error,
        "event_localization_error_max_km": max(
            upper["event_localization_error_max_km"],
            ground["event_localization_error_max_km"],
            grid_exit["event_localization_error_max_km"],
            reflected["event_localization_error_max_km"],
        ),
        "chapman_apex_height_km": reflected["apex_height_km"],
        "chapman_apex_radial_group_velocity_km_s": (
            reflected["apex_radial_group_velocity_km_s"]),
        "chapman_landing_height_km": reflected["landing_height_km"],
        "chapman_hamiltonian_residual_max": reflected[
            "hamiltonian_residual_max"],
    }
    checks = {
        "upper_boundary_status": upper["termination_status"] == "upper_boundary",
        "upper_boundary_dense_location": (
            upper_error <= 1.0e-5
            and metrics["vacuum_upper_height_error_km"] <= 1.0e-8),
        "ground_status_after_descent": ground["termination_status"] == "ground",
        "ground_dense_location": (
            ground_error <= 1.0e-5
            and metrics["vacuum_ground_height_error_km"] <= 1.0e-8),
        "surface_departure_guard": surface_downward_rejected,
        "grid_exit_distinct_status": grid_exit["termination_status"] == "grid_exit",
        "grid_exit_dense_location": grid_exit_height_error <= 1.0e-8,
        "chapman_reflected_ground": reflected["termination_status"] == "ground",
        "chapman_turning_apex": (
            reflected["turning_point_count"] >= 1
            and reflected["apex_height_km"] is not None
            and abs(reflected["apex_radial_group_velocity_km_s"]) <= 1.0e-3),
        "chapman_ground_height": (
            reflected["landing_height_km"] is not None
            and abs(reflected["landing_height_km"]) <= 1.0e-8),
        "event_localization_precision": (
            metrics["event_localization_error_max_km"] <= 1.0e-8),
        "explicit_trajectory_semantics": (
            "trajectory_3d" not in reflected
            and "trajectory_ecef_km" in reflected
            and "trajectory_geocentric_lat_lon_height" in reflected
            and "trajectory_physical_height_km" in reflected),
        "explicit_path_and_delay_semantics": all(
            key in reflected for key in (
                "geom_path_km", "group_path_km", "tau_ms")),
        "landing_semantics_only_for_ground": (
            upper["landing_lat_deg"] is None
            and grid_exit["landing_lat_deg"] is None
            and ground["landing_lat_deg"] is not None),
        "event_history_distinguishes_events": (
            upper["event_history"][-1]["event"] == "upper_boundary"
            and upper["event_history"][-1]["direction"] == "outward"
            and ground["event_history"][-1]["event"] == "ground"
            and ground["event_history"][-1]["direction"] == "descending"
            and grid_exit["event_history"][-1]["event"] == "grid_exit"
            and grid_exit["event_history"][-1]["direction"]
            == "provider_domain_exit"
            and reflected["event_history"][0]["event"] == "turning_apex"),
        "traceability_hashes_recorded": (
            len(reflected["provider_hash"]) == 64
            and len(reflected["coordinate_mapping_hash"]) == 64),
        "constraint_preserved": (
            reflected["hamiltonian_residual_max"] <= 1.0e-8),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    report = {
        "schema_version": "a4-c0r4-events-v1",
        "stage": "A4-C0R-4",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "spherical_ECEF_physical_events_and_outputs_only",
        "earth_radius_km": R_KM,
        "metrics": metrics,
        "checks": checks,
        "termination_statuses": {
            "vacuum_upper": upper["termination_status"],
            "vacuum_ground": ground["termination_status"],
            "grid_exit": grid_exit["termination_status"],
            "chapman_reflected": reflected["termination_status"],
        },
        "source_sha256": {
            "ecef_ah3d": _sha256(
                ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
            "provider": _sha256(
                ROOT / "models" / "ne_provider_geocentric.py"),
            "geometry": _sha256(ROOT / "utils.py"),
            "test": _sha256(Path(__file__)),
        },
        "claim_boundary": (
            "Spherical ground/upper/grid-exit/turning events and explicit ECEF/"
            "geocentric path outputs only; P2P and downstream production callers "
            "are not migrated until C0R-5"
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-4: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
