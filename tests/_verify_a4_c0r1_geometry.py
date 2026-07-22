"""Acceptance checks for the C0R-1 spherical/ECEF geometry kernel."""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from utils import (ecef_to_geocentric, elevation_bearing_from_ecef, enu_basis,
                   geocentric_to_ecef, launch_direction_ecef, rx_enu_residual,
                   spherical_destination)


REPORT = ROOT / "output2" / "a4_pharlap" / "c0r1_geometry_validation.json"
R_KM = 6371.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _angle_error_deg(actual: float, expected: float) -> float:
    return abs((actual - expected + 180.0) % 360.0 - 180.0)


def _invalid_inputs_rejected() -> bool:
    calls = (
        lambda: geocentric_to_ecef(91.0, 0.0, 0.0, R_KM),
        lambda: geocentric_to_ecef(0.0, 0.0, -R_KM, R_KM),
        lambda: ecef_to_geocentric([0.0, 0.0, 0.0], R_KM),
        lambda: enu_basis(0.0, np.nan, R_KM),
        lambda: launch_direction_ecef(0.0, 0.0, 91.0, 0.0, R_KM),
        lambda: elevation_bearing_from_ecef([0.0, 0.0, 0.0], 0.0, 0.0, R_KM),
        lambda: spherical_destination(0.0, 0.0, 0.0, np.nan, R_KM),
        lambda: rx_enu_residual([1.0, 2.0], 0.0, 0.0, 0.0, R_KM),
        lambda: geocentric_to_ecef(0.0, 0.0, 0.0, -R_KM),
    )
    rejected = 0
    for call in calls:
        try:
            call()
        except ValueError:
            rejected += 1
    return rejected == len(calls)


def run() -> int:
    lat = np.array([-89.999, -45.0, 0.0, 30.0, 89.999])
    lon = np.array([-179.9, -90.0, 0.0, 120.0, 179.9])
    height = np.array([0.0, 60.0, 300.0, 550.0, 1000.0])
    ecef = geocentric_to_ecef(lat, lon, height, R_KM)
    recovered = ecef_to_geocentric(ecef, R_KM)
    roundtrip_ecef = geocentric_to_ecef(
        recovered[:, 0], recovered[:, 1], recovered[:, 2], R_KM)
    roundtrip_error = np.linalg.norm(roundtrip_ecef - ecef, axis=1)

    orthogonality_error = 0.0
    minimum_determinant = 1.0
    elevation_error = 0.0
    bearing_error = 0.0
    launch_cases = (
        (30.0, 120.0, 35.0, 37.0),
        (-45.0, 179.9, 5.0, 359.0),
        (70.0, -170.0, -10.0, 181.0),
    )
    for lat0, lon0, elevation0, bearing0 in launch_cases:
        basis = enu_basis(lat0, lon0, R_KM)
        orthogonality_error = max(
            orthogonality_error,
            float(np.max(np.abs(basis.T @ basis - np.eye(3)))),
        )
        minimum_determinant = min(minimum_determinant, float(np.linalg.det(basis)))
        direction = launch_direction_ecef(
            lat0, lon0, elevation0, bearing0, R_KM)
        elevation1, bearing1 = elevation_bearing_from_ecef(
            direction, lat0, lon0, R_KM)
        elevation_error = max(elevation_error, abs(elevation1 - elevation0))
        bearing_error = max(bearing_error, _angle_error_deg(bearing1, bearing0))

    endpoint_lat, endpoint_lon = spherical_destination(
        0.0, 0.0, 90.0, 0.5 * np.pi * R_KM, R_KM)
    endpoint = geocentric_to_ecef(endpoint_lat, endpoint_lon, 0.0, R_KM)
    expected_endpoint = geocentric_to_ecef(0.0, 90.0, 0.0, R_KM)
    endpoint_error = float(np.linalg.norm(endpoint - expected_endpoint))
    periodic_error = float(np.linalg.norm(
        geocentric_to_ecef(30.0, 181.0, 100.0, R_KM)
        - geocentric_to_ecef(30.0, -179.0, 100.0, R_KM)))
    _, dateline_lon = spherical_destination(10.0, 179.5, 90.0, 200.0, R_KM)
    signed_a = spherical_destination(30.0, 120.0, 37.0, -50.0, R_KM)
    signed_b = spherical_destination(30.0, 120.0, 217.0, 50.0, R_KM)
    signed_distance_error = float(np.linalg.norm(
        geocentric_to_ecef(*signed_a, 0.0, R_KM)
        - geocentric_to_ecef(*signed_b, 0.0, R_KM)))

    rx_lat, rx_lon, rx_height = 32.5, 120.0, 0.0
    rx = geocentric_to_ecef(rx_lat, rx_lon, rx_height, R_KM)
    basis = enu_basis(rx_lat, rx_lon, R_KM)
    displaced = rx + basis @ np.array([1.25, -2.5, 0.3])
    residual_error = float(np.max(np.abs(
        rx_enu_residual(displaced, rx_lat, rx_lon, rx_height, R_KM)
        - np.array([1.25, -2.5]))))

    metrics = {
        "coordinate_roundtrip_max_error_km": float(np.max(roundtrip_error)),
        "enu_orthogonality_max_error": orthogonality_error,
        "enu_minimum_determinant": minimum_determinant,
        "launch_elevation_recovery_max_error_deg": elevation_error,
        "launch_bearing_recovery_max_error_deg": bearing_error,
        "known_endpoint_error_km": endpoint_error,
        "longitude_periodic_position_error_km": periodic_error,
        "dateline_endpoint_lon_deg": dateline_lon,
        "signed_distance_equivalence_error_km": signed_distance_error,
        "rx_enu_residual_max_error_km": residual_error,
    }
    checks = {
        "coordinate_roundtrip": metrics["coordinate_roundtrip_max_error_km"] <= 1.0e-9,
        "enu_orthogonality": metrics["enu_orthogonality_max_error"] <= 1.0e-12,
        "enu_right_handed": metrics["enu_minimum_determinant"] > 0.0,
        "launch_elevation_recovery": elevation_error <= 1.0e-10,
        "launch_bearing_recovery": bearing_error <= 1.0e-10,
        "known_spherical_endpoint": endpoint_error <= 1.0e-9,
        "longitude_periodicity": periodic_error <= 1.0e-9,
        "longitude_normalized": -180.0 <= dateline_lon < 180.0,
        "signed_surface_distance": signed_distance_error <= 1.0e-9,
        "rx_enu_residual": residual_error <= 1.0e-12,
        "invalid_inputs_rejected": _invalid_inputs_rejected(),
    }
    report = {
        "schema_version": "a4-c0r1-geometry-v1",
        "stage": "A4-C0R-1",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "spherical_geocentric_ecef_geometry_only",
        "earth_model": "sphere_user_radius",
        "earth_radius_km": R_KM,
        "units": {"position": "km", "height": "km", "angles": "degree"},
        "coordinate_semantics": {
            "latitude": "geocentric",
            "longitude": "east_positive_periodic",
            "height": "above_spherical_earth",
            "bearing": "clockwise_from_geocentric_north",
            "enu_basis_order": "east,north,up",
        },
        "metrics": metrics,
        "checks": checks,
        "source_sha256": {
            "utils.py": _sha256(ROOT / "utils.py"),
            "test": _sha256(Path(__file__)),
        },
        "claim_boundary": (
            "geometry kernel only; no Provider, AH3D, P2P, PE, scattering, "
            "power, WGS84, or PHaRLAP trajectory validation"
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-1: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
