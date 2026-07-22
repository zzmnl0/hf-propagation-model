"""Acceptance checks for C0R-5 P2P and explicit-ECEF downstream consumers."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.absorption_collision import complex_absorption_ecef_ray_dB
from legacy.models.local_frame import make_ray_local_frame, sample_Ne_local_cross_section
from legacy.models.hybrid_ray_pe import apply_ray_pe_coupling_3d
from models.ne_provider_geocentric import (
    ChapmanGeocentricNeProvider3D,
    GeocentricGridNeProvider3D,
)
from models.point_to_point_ah_3d_ecef import find_seeded_rays_ah3d_ecef
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)
from models.scatter_coupler import apply_scatter_regime_3d
from viz.plot_utils import _mode_xz


R_KM = 6371.0
REPORT = ROOT / "output2" / "a4_pharlap" / "c0r5_downstream_validation.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> int:
    chapman = ChapmanGeocentricNeProvider3D(
        R_KM, foF2_MHz=7.0, hmF2_km=350.0, ymF2_km=100.0)
    model = RefractiveIndexAH3DECEF(
        chapman, earth_radius_km=R_KM, freq_MHz=10.0,
        wave_mode="O", gyro_MHz=0.0, b_ecef_hat=(0.0, 0.0, 1.0),
        integration_profile="smoke")
    tx = np.array([30.0, 120.0, 0.0])
    reference = trace_ray_ah3d_ecef(
        model, tx_lat_deg=tx[0], tx_lon_deg=tx[1], tx_height_km=tx[2],
        elevation_deg=35.0, bearing_deg=37.0, h_stop_km=600.0,
        max_path_km=3000.0, profile="smoke")
    rx = np.asarray(reference["landing_geocentric_lat_lon_height"], dtype=float)
    roots = find_seeded_rays_ah3d_ecef(
        model, tx_geocentric=tx, rx_geocentric=rx,
        seeds=[{
            "seed_id": "legacy-35deg",
            "elevation_deg": 35.0,
            "bearing_deg": 37.0,
            "semantics": "legacy_coordinate_reference_seed_only",
        }],
        h_stop_km=600.0, max_path_km=3000.0,
        profile="smoke", tol_km=0.1, max_nfev=4,
        d_elevation_deg=0.05, d_bearing_deg=0.05)
    ray = roots[0] if roots else None

    if ray is None:
        metrics = {"p2p_root_count": 0}
        checks = {"p2p_root_found": False}
    else:
        frame = make_ray_local_frame(ray)
        basis = np.asarray(frame["basis_matrix"], dtype=float)
        frame_orthogonality = float(np.max(np.abs(
            basis @ basis.T - np.eye(3))))
        frame_handedness = float(np.linalg.det(basis))
        cross_section = sample_Ne_local_cross_section(
            chapman, ray, s_offsets_km=[-1.0, 0.0, 1.0],
            u_offsets_km=[-1.0, 0.0, 1.0])
        legacy_pe_grid_rejected = False
        try:
            apply_ray_pe_coupling_3d(
                None, ray, np.zeros((2, 2)), np.array([0.0, 1.0]),
                np.array([0.0, 1.0]), 1.0, model.k0,
                {"enable": True, "center_ecef_km": ray["apex_ecef_km"],
                 "b_ecef_hat": ray["b_ecef_hat"],
                 "sigma_perp_km": 1.0e6, "sigma_para_km": 1.0e6})
        except ValueError as exc:
            legacy_pe_grid_rejected = "cannot be interpreted as ECEF" in str(exc)

        latitude = np.array([-90.0, 0.0, 90.0])
        longitude = np.array([-180.0, -90.0, 0.0, 90.0])
        height = np.array([0.0, 300.0, 600.0])
        constant = GeocentricGridNeProvider3D(
            np.full((3, 4, 3), 1.0e11), latitude, longitude, height,
            R_KM, longitude_periodic_complete=True)
        absorption_one = complex_absorption_ecef_ray_dB(
            ray, constant,
            {"collision_model": "constant", "nu0_s": 1.0e5})
        absorption_two = complex_absorption_ecef_ray_dB(
            ray, constant,
            {"collision_model": "constant", "nu0_s": 1.0e5},
            two_way=True)

        scatter_modes = apply_scatter_regime_3d(
            [dict(ray)], chapman,
            {"Cs_rel": 1.0e-8, "p": 3.0, "L_outer_km": 10.0,
             "DeltaL_km": 100.0, "z_eff_km": 300.0,
             "Ne0_m3": 1.0e11},
            10.0, 1000.0, 1.0, 1.0, enable_mps=False)
        plot_x, plot_z = _mode_xz(ray)

        metrics = {
            "p2p_root_count": len(roots),
            "p2p_landing_error_km": ray["p2p_landing_error_km"],
            "p2p_jacobian_max_condition": max(
                item["condition"] for item in ray["p2p_jacobian_reference"]),
            "frame_orthogonality_error": frame_orthogonality,
            "frame_handedness": frame_handedness,
            "cross_section_min_Ne_m3": float(np.min(cross_section["Ne_m3"])),
            "absorption_one_way_dB": absorption_one["abs_total_dB"],
            "absorption_two_way_ratio": (
                absorption_two["abs_total_dB"]
                / max(absorption_one["abs_total_dB"], 1.0e-30)),
            "plot_final_cumulative_ground_range_km": float(plot_x[-1]),
            "plot_max_physical_height_km": float(max(plot_z)),
        }
        checks = {
            "p2p_root_found": len(roots) == 1,
            "p2p_RX_ENU_residual": ray["p2p_landing_error_km"] < 0.1,
            "p2p_multiscale_centered_jacobian": (
                ray["p2p_jacobian_method"] == "multiscale_centered_RX_ENU"
                and all(item["stable"]
                        for item in ray["p2p_jacobian_reference"])),
            "legacy_root_is_seed_only": (
                ray["p2p_seed_semantics"]
                == "legacy_coordinate_reference_seed_only"
                and ray["p2p_search_scope"]
                == "seeded_local_roots_not_full_domain_mode_audit"),
            "explicit_ECEF_dedup_semantics": (
                ray["p2p_dedup_semantics"]
                == "branch_status_turning_launch_ECEF_trajectory_apex_group_delay"),
            "local_frame_ECEF_radial_up": (
                frame["frame_coordinate_model"] == "ECEF_ray_tangent_radial_up"
                and "origin_ecef_km" in frame and "origin_km" not in frame),
            "local_frame_orthonormal_right_handed": (
                frame_orthogonality <= 1.0e-12
                and frame_handedness > 0.0),
            "local_cross_section_ECEF_provider": (
                cross_section["coordinate_model"]
                == "ECEF_ray_tangent_radial_up"
                and "ecef_x_km" in cross_section
                and "x_km" not in cross_section
                and np.all(np.isfinite(cross_section["Ne_m3"]))),
            "RayPE_rejects_legacy_grid_for_ECEF_ray": legacy_pe_grid_rejected,
            "absorption_uses_physical_height_ECEF": (
                absorption_one["abs_path_coordinate_model"]
                == "spherical_ECEF_physical_height"
                and absorption_one["abs_total_dB"] > 0.0),
            "absorption_two_way_semantics": abs(
                metrics["absorption_two_way_ratio"] - 2.0) <= 1.0e-12,
            "scatter_samples_ECEF_provider": (
                scatter_modes[0]["scatter_path_coordinate_model"]
                == "spherical_ECEF"),
            "plot_uses_geocentric_physical_height": (
                len(plot_x) == len(ray["trajectory_ecef_km"])
                and max(plot_z) == ray["apex_height_km"]
                and plot_x[-1] > 0.0),
            "formal_ray_has_no_legacy_aliases": all(
                key not in ray for key in (
                    "trajectory_3d", "x_land_km", "y_land_km",
                    "z_land_km", "h_reflect_km")),
            "traceability_hashes_recorded": (
                len(ray["provider_hash"]) == 64
                and len(ray["coordinate_mapping_hash"]) == 64),
        }
    checks = {name: bool(value) for name, value in checks.items()}
    report = {
        "schema_version": "a4-c0r5-downstream-v1",
        "stage": "A4-C0R-5",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "seeded_ECEF_P2P_and_explicit_downstream_consumers",
        "metrics": metrics,
        "checks": checks,
        "source_sha256": {
            "ecef_ray": _sha256(ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
            "ecef_p2p": _sha256(ROOT / "models" / "point_to_point_ah_3d_ecef.py"),
            "local_frame": _sha256(ROOT / "legacy" / "models" / "local_frame.py"),
            "ray_pe": _sha256(ROOT / "legacy" / "models" / "hybrid_ray_pe.py"),
            "absorption": _sha256(ROOT / "models" / "absorption_collision.py"),
            "scatter": _sha256(ROOT / "models" / "scatter_coupler.py"),
            "plot": _sha256(ROOT / "viz" / "plot_utils.py"),
            "test": _sha256(Path(__file__)),
        },
        "claim_boundary": (
            "Seeded local P2P roots and explicit ECEF downstream coordinate "
            "consumption only; not full-domain mode completeness, main entrypoint "
            "migration, full 3-D PE/full-wave, or calibrated power validation"
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-5: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
