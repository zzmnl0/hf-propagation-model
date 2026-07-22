"""Acceptance checks for the C0R-3 strict ECEF AH Hamiltonian core."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.dispersion_ah import C_KM_S
from models.ne_provider_geocentric import (
    ChapmanGeocentricNeProvider3D,
    GeocentricNeProvider3D,
    VacuumGeocentricNeProvider3D,
)
from models.ray_tracer_ah_3d_ecef import (
    COORDINATE_MODEL,
    RefractiveIndexAH3DECEF,
    SpatialDerivativeError,
    group_velocity_ah3d_ecef,
    hamiltonian_residual_ah3d_ecef,
    integrate_fixed_path_ah3d_ecef,
    launch_state_ah3d_ecef,
    project_dispersion_ah3d_ecef,
    spatial_gradient_ah3d_ecef,
)
from utils import geocentric_to_ecef, launch_direction_ecef


R_KM = 6371.0
REPORT = ROOT / "output2" / "a4_pharlap" / "c0r3_ecef_hamiltonian_validation.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class _UnresolvedProvider(GeocentricNeProvider3D):
    name = "test_unresolved_ecef_provider"

    def __init__(self, center_x: float):
        self.earth_radius_km = R_KM
        self.center_x = float(center_x)

    def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
        return 1.0e11

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        if float(earth_radius_km) != R_KM:
            raise ValueError("radius mismatch")
        offset = float(np.asarray(ecef_km)[0] - self.center_x)
        return 1.0e11 * (1.0 + 0.1 * np.sin(offset / 0.3))

    def get_metadata(self):
        return {**super().get_metadata(), "earth_radius_km": R_KM,
                "provider_hash": "unresolved-test"}


class _DomainRejectingProvider(_UnresolvedProvider):
    name = "test_domain_rejecting_ecef_provider"

    def get_Ne_ecef(self, ecef_km, earth_radius_km, t=0.0):
        offset = abs(float(np.asarray(ecef_km)[0] - self.center_x))
        if offset > 0.1:
            raise ValueError("outside frozen Provider domain")
        return 1.0e11


def _model(provider, wave_mode="O", gyro_MHz=0.0, **kwargs):
    return RefractiveIndexAH3DECEF(
        provider, earth_radius_km=R_KM, freq_MHz=10.0,
        wave_mode=wave_mode, gyro_MHz=gyro_MHz,
        b_ecef_hat=(0.2, -0.3, 0.9327379053), **kwargs)


def run() -> int:
    vacuum = VacuumGeocentricNeProvider3D(R_KM)
    vacuum_o = _model(vacuum, "O", gyro_MHz=1.2)
    vacuum_x = _model(vacuum, "X", gyro_MHz=1.2)
    position = geocentric_to_ecef(30.0, 120.0, 0.0, R_KM)
    direction = launch_direction_ecef(30.0, 120.0, 35.0, 37.0, R_KM)
    k_vacuum = vacuum_o.k0 * direction
    state_vacuum = np.r_[position, k_vacuum, 0.0]
    n2_o = vacuum_o.n2_raw(position, k_vacuum)
    n2_x = vacuum_x.n2_raw(position, k_vacuum)
    vacuum_residual = hamiltonian_residual_ah3d_ecef(
        state_vacuum, vacuum_o)
    vacuum_speed = float(np.linalg.norm(group_velocity_ah3d_ecef(
        state_vacuum, vacuum_o)))

    vacuum_result = integrate_fixed_path_ah3d_ecef(
        vacuum_o, tx_lat_deg=30.0, tx_lon_deg=120.0, tx_height_km=0.0,
        elevation_deg=35.0, bearing_deg=37.0, path_length_km=250.0,
        profile="smoke")
    trajectory = np.asarray(vacuum_result["trajectory_ecef_km"])
    path_samples = np.asarray(vacuum_result["trajectory_geom_path_km"])
    analytic = position[None, :] + path_samples[:, None] * direction[None, :]
    line_error = np.linalg.norm(trajectory - analytic, axis=1)
    path_relative = abs(
        vacuum_result["group_path_km"] - vacuum_result["geom_path_km"]
    ) / vacuum_result["geom_path_km"]

    degenerate_o = _model(vacuum, "O", gyro_MHz=0.0)
    degenerate_x = _model(vacuum, "X", gyro_MHz=0.0)
    b0_difference = abs(
        degenerate_o.n2_raw(position, k_vacuum)
        - degenerate_x.n2_raw(position, k_vacuum))

    chapman = ChapmanGeocentricNeProvider3D(
        R_KM, foF2_MHz=7.0, hmF2_km=350.0, ymF2_km=100.0)
    chapman_model = _model(chapman, "O", gyro_MHz=0.0)
    chapman_state = launch_state_ah3d_ecef(
        chapman_model, tx_lat_deg=30.0, tx_lon_deg=120.0,
        tx_height_km=250.0, elevation_deg=40.0, bearing_deg=0.0)
    gradient, gradient_diagnostics = spatial_gradient_ah3d_ecef(
        chapman_state, chapman_model)
    perturbed = chapman_state.copy()
    perturbed[3:6] *= 1.05
    before_projection = hamiltonian_residual_ah3d_ecef(
        perturbed, chapman_model)
    projected, _, projection_ok = project_dispersion_ah3d_ecef(
        perturbed, chapman_model, tolerance=1.0e-10)
    after_projection = hamiltonian_residual_ah3d_ecef(
        projected, chapman_model)
    chapman_result = integrate_fixed_path_ah3d_ecef(
        chapman_model, tx_lat_deg=30.0, tx_lon_deg=120.0,
        tx_height_km=250.0, elevation_deg=40.0, bearing_deg=0.0,
        path_length_km=2.0, profile="smoke")
    chapman_point = chapman_state[:3]
    chapman_geocentric = chapman.get_Ne_geocentric(30.0, 120.0, 250.0)
    chapman_ecef = chapman.get_Ne_ecef(chapman_point, R_KM)

    type_rejected = False
    try:
        _model(object())
    except TypeError:
        type_rejected = True
    radius_rejected = False
    try:
        RefractiveIndexAH3DECEF(
            vacuum, earth_radius_km=R_KM + 1.0, freq_MHz=10.0)
    except ValueError:
        radius_rejected = True

    unresolved_rejected = False
    unresolved = _model(_UnresolvedProvider(position[0]), gyro_MHz=0.0)
    unresolved_state = launch_state_ah3d_ecef(
        unresolved, tx_lat_deg=30.0, tx_lon_deg=120.0,
        tx_height_km=0.0, elevation_deg=35.0, bearing_deg=37.0)
    try:
        spatial_gradient_ah3d_ecef(unresolved_state, unresolved)
    except SpatialDerivativeError:
        unresolved_rejected = True

    domain_rejected = False
    domain_model = _model(_DomainRejectingProvider(position[0]), gyro_MHz=0.0)
    domain_state = launch_state_ah3d_ecef(
        domain_model, tx_lat_deg=30.0, tx_lon_deg=120.0,
        tx_height_km=0.0, elevation_deg=35.0, bearing_deg=37.0)
    try:
        spatial_gradient_ah3d_ecef(domain_state, domain_model)
    except SpatialDerivativeError:
        domain_rejected = True

    metrics = {
        "vacuum_n2_O_error": abs(n2_o - 1.0),
        "vacuum_n2_X_error": abs(n2_x - 1.0),
        "vacuum_initial_hamiltonian_residual": vacuum_residual,
        "vacuum_group_speed_relative_error": abs(vacuum_speed / C_KM_S - 1.0),
        "vacuum_straight_line_max_error_km": float(np.max(line_error)),
        "vacuum_group_geom_path_relative_difference": path_relative,
        "vacuum_hamiltonian_residual_max": vacuum_result[
            "hamiltonian_residual_max"],
        "B0_OX_n2_difference": b0_difference,
        "chapman_same_point_Ne_relative_difference": abs(
            chapman_geocentric - chapman_ecef) / chapman_geocentric,
        "chapman_spatial_gradient_norm": float(np.linalg.norm(gradient)),
        "chapman_fd_relative_spread": gradient_diagnostics["relative_spread"],
        "projection_residual_before": before_projection,
        "projection_residual_after": after_projection,
    }
    checks = {
        "strict_provider_type_rejected": type_rejected,
        "earth_radius_conflict_rejected": radius_rejected,
        "no_local_or_earth_flat_entry": (
            not hasattr(vacuum_o, "earth_flat")
            and not hasattr(vacuum_o, "n3d")
            and vacuum_result["coordinate_model"] == COORDINATE_MODEL),
        "vacuum_O_n2_unity": metrics["vacuum_n2_O_error"] <= 1.0e-12,
        "vacuum_X_n2_unity": metrics["vacuum_n2_X_error"] <= 1.0e-12,
        "vacuum_D_surface": vacuum_residual <= 1.0e-8,
        "vacuum_group_velocity": (
            metrics["vacuum_group_speed_relative_error"] <= 1.0e-7),
        "vacuum_fixed_path_completed": (
            vacuum_result["termination_status"]
            == "fixed_geometric_path_completed"),
        "vacuum_straight_line": (
            metrics["vacuum_straight_line_max_error_km"] <= 1.0e-2),
        "vacuum_group_equals_geom_path": path_relative <= 1.0e-4,
        "vacuum_constraint_preserved": (
            metrics["vacuum_hamiltonian_residual_max"] <= 1.0e-8),
        "B0_OX_degenerate": b0_difference <= 1.0e-12,
        "chapman_same_geocentric_ecef_point": (
            metrics["chapman_same_point_Ne_relative_difference"] <= 1.0e-14),
        "chapman_gradient_finite_and_stable": (
            np.all(np.isfinite(gradient))
            and gradient_diagnostics["relative_spread"]
            <= chapman_model.spatial_fd_stability_rtol),
        "dispersion_projection": (
            projection_ok and before_projection > 1.0e-4
            and after_projection <= 1.0e-8),
        "nonvacuum_traceability_hashes_recorded": (
            chapman_result["termination_status"]
            == "fixed_geometric_path_completed"
            and chapman_result["provider_hash"] == chapman_model.provider_hash
            and len(chapman_result["provider_hash"]) == 64
            and len(chapman_result["coordinate_mapping_hash"]) == 64),
        "unresolved_multiscale_gradient_rejected": unresolved_rejected,
        "provider_domain_crossing_rejected": domain_rejected,
        "traceability_hashes_recorded": (
            len(vacuum_result["provider_hash"]) == 64
            and len(vacuum_result["coordinate_mapping_hash"]) == 64),
        "formal_output_has_no_ambiguous_trajectory_3d": (
            "trajectory_3d" not in vacuum_result
            and "trajectory_ecef_km" in vacuum_result),
    }
    report = {
        "schema_version": "a4-c0r3-ecef-hamiltonian-v1",
        "stage": "A4-C0R-3",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "strict_ecef_AH_Hamiltonian_fixed_path_core_only",
        "coordinate_model": COORDINATE_MODEL,
        "earth_radius_km": R_KM,
        "metrics": metrics,
        "checks": checks,
        "source_sha256": {
            "ecef_ah3d": _sha256(
                ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
            "provider": _sha256(
                ROOT / "models" / "ne_provider_geocentric.py"),
            "dispersion": _sha256(ROOT / "models" / "dispersion_ah.py"),
            "geometry": _sha256(ROOT / "utils.py"),
            "test": _sha256(Path(__file__)),
        },
        "claim_boundary": (
            "ECEF Hamiltonian, multiscale spatial derivatives, DOP853 fixed-path "
            "integration, projection, and vacuum/group-velocity checks only; "
            "physical events and downstream callers are not migrated until C0R-4/5"
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-3: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
