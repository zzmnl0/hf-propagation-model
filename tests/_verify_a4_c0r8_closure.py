"""C0R-8 regression, convergence, and formal-entry closure gate."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from models.ne_provider_geocentric import ChapmanGeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    integrate_fixed_path_ah3d_ecef,
)


R_KM = 6371.0
OUT = ROOT / "output2" / "a4_pharlap"
REPORT = OUT / "a4_c0r8_validation.json"
LOG_DIR = OUT / "c0r8_logs"

REGRESSIONS = (
    "tests/_verify_a4_c0r0_freeze.py",
    "tests/_verify_a4_c0r1_geometry.py",
    "tests/_verify_a4_c0r2_provider.py",
    "tests/_verify_a4_c0r3_ecef_hamiltonian.py",
    "tests/_verify_a4_c0r4_events.py",
    "tests/_verify_a4_c0r5_downstream.py",
    "tests/_verify_a4_c0r6_entry_migration.py",
    "tests/_verify_a4_c0r7_pharlap_vacuum.py",
    "tests/_verify_main_output_semantics.py",
    "tests/_verify_main2_scenario_fields.py",
)

FORMAL_FILES = (
    "main.py", "main2.py",
    "models/unified_hybrid_pipeline.py",
    "models/scenario_provider_factory.py",
    "models/ne_provider_geocentric.py",
    "models/ray_tracer_ah_3d_ecef.py",
    "models/point_to_point_ah_3d_ecef.py",
    "models/absorption_collision.py",
    "models/scatter_coupler.py", "models/scatter_phase_screen.py",
    "viz/plot_utils.py",
)

HISTORICAL_REFERENCE_ONLY = (
    "legacy/models/ray_tracer.py", "legacy/models/point_to_point.py",
    "legacy/models/ne_provider.py", "legacy/models/parameterized_ionosphere_3d.py",
    "legacy/models/ray_tracer_3d.py", "legacy/models/point_to_point_3d.py",
    "legacy/models/ray_tracer_ah_3d.py", "legacy/models/ray_tracer_ah_3d_v2.py",
    "legacy/models/point_to_point_ah_3d.py", "legacy/models/p2p_mode_audit_3d.py",
    "legacy/models/hybrid_model.py", "legacy/models/local_frame.py",
    "legacy/models/hybrid_ray_pe.py", "legacy/models/pe_propagator.py",
    "legacy/models/tube_tracer.py", "legacy/models/tube_tracer_3d.py",
    "legacy/models/radar_clutter.py", "legacy/analysis/parameter_sweep_3d.py",
    "legacy/viz/plot_ne_background.py", "legacy/_bench.py",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def _run_regressions() -> dict:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for relative in REGRESSIONS:
        started = time.perf_counter()
        process = subprocess.run(
            [sys.executable, str(ROOT / relative)], cwd=ROOT,
            text=True, capture_output=True, timeout=1800)
        elapsed = time.perf_counter() - started
        log = LOG_DIR / (Path(relative).stem + ".log")
        log.write_text(process.stdout + "\n" + process.stderr, encoding="utf-8")
        results[relative] = {
            "passed": process.returncode == 0,
            "returncode": process.returncode,
            "elapsed_s": elapsed,
            "log": str(log), "log_sha256": _sha256(log),
        }
        print(f" {relative}: {'OK' if process.returncode == 0 else 'FAIL'} "
              f"({elapsed:.1f}s)")
    return results


def _resample(result: dict, sample: np.ndarray) -> np.ndarray:
    path = np.asarray(result["trajectory_geom_path_km"], dtype=float)
    points = np.asarray(result["trajectory_ecef_km"], dtype=float)
    return np.column_stack([
        np.interp(sample, path, points[:, axis]) for axis in range(3)])


def _step_convergence() -> tuple[dict, dict]:
    provider = ChapmanGeocentricNeProvider3D(
        R_KM, foF2_MHz=7.0, hmF2_km=350.0, ymF2_km=100.0)
    model = RefractiveIndexAH3DECEF(
        provider, earth_radius_km=R_KM, freq_MHz=10.0,
        wave_mode="O", gyro_MHz=0.0, b_ecef_hat=(0.0, 0.0, 1.0),
        integration_profile="reference")
    half_name = "c0r8_reference_half_step"
    base = dict(cfg.UPDATE2["AH3D_INTEGRATION"]["reference"])
    half = dict(base)
    half["ds_min_km"] = 0.5 * float(base["ds_min_km"])
    half["ds_max_km"] = 0.5 * float(base["ds_max_km"])
    cfg.UPDATE2["AH3D_INTEGRATION"][half_name] = half
    try:
        common = dict(
            tx_lat_deg=30.0, tx_lon_deg=120.0, tx_height_km=0.0,
            elevation_deg=35.0, bearing_deg=37.0, path_length_km=250.0)
        reference = integrate_fixed_path_ah3d_ecef(
            model, profile="reference", **common)
        refined = integrate_fixed_path_ah3d_ecef(
            model, profile=half_name, **common)
    finally:
        cfg.UPDATE2["AH3D_INTEGRATION"].pop(half_name, None)
    sample = np.linspace(0.0, 250.0, 201)
    error = np.linalg.norm(
        _resample(reference, sample) - _resample(refined, sample), axis=1)
    group_relative = abs(
        float(reference["group_path_km"]) - float(refined["group_path_km"])
    ) / max(float(refined["group_path_km"]), 1.0e-12)
    metrics = {
        "reference_ds_max_km": float(base["ds_max_km"]),
        "refined_ds_max_km": float(half["ds_max_km"]),
        "trajectory_max_change_km": float(np.max(error)),
        "trajectory_rms_change_km": float(np.sqrt(np.mean(error * error))),
        "group_path_relative_change": group_relative,
        "reference_hamiltonian_residual_max": float(
            reference["hamiltonian_residual_max"]),
        "refined_hamiltonian_residual_max": float(
            refined["hamiltonian_residual_max"]),
    }
    checks = {
        "reference_fixed_path_completed": (
            reference["termination_status"] == "fixed_geometric_path_completed"),
        "refined_fixed_path_completed": (
            refined["termination_status"] == "fixed_geometric_path_completed"),
        "step_halved": metrics["refined_ds_max_km"] == 0.5 * metrics["reference_ds_max_km"],
        "trajectory_step_convergence": metrics["trajectory_max_change_km"] <= 1.0e-3,
        "group_path_step_convergence": group_relative <= 1.0e-5,
        "reference_constraint": metrics["reference_hamiltonian_residual_max"] <= 1.0e-8,
        "refined_constraint": metrics["refined_hamiltonian_residual_max"] <= 1.0e-8,
    }
    return metrics, checks


def _static_closure() -> tuple[dict, dict]:
    forbidden = (
        "trajectory_3d", "earth_flat", "legacy_parabolic",
        "cfg.TX_POS", "cfg.RX_POS",
        "from .point_to_point_3d", "from legacy.models.point_to_point_3d",
        "from .ray_tracer_3d", "from legacy.models.ray_tracer_3d",
    )
    violations = []
    hashes = {}
    for relative in FORMAL_FILES:
        path = ROOT / relative
        text = path.read_text(encoding="utf-8")
        hashes[relative] = _sha256(path)
        violations.extend(
            {"path": relative, "token": token}
            for token in forbidden if token in text)
    main_text = (ROOT / "main.py").read_text(encoding="utf-8")
    config_text = (ROOT / "config.py").read_text(encoding="utf-8")
    checks = {
        "formal_files_exist": all((ROOT / path).is_file() for path in FORMAL_FILES),
        "formal_files_have_no_forbidden_legacy_tokens": not violations,
        "formal_geocentric_config_present": all(
            token in config_text for token in (
                "TX_GEOCENTRIC", "RX_GEOCENTRIC", "EARTH_MODEL")),
        "unreachable_legacy_tube_implementations_removed": all(
            token not in main_text for token in (
                "_run_radar_clutter_legacy_coordinate_reference_only",
                "_run_tube_radar_legacy_coordinate_reference_only")),
        "legacy_radar_products_remain_explicitly_rejected": all(
            token in main_text for token in (
                "def run_radar_clutter():", "def run_tube_radar():",
                "formal entry rejects the unported link-local")),
        "historical_reference_files_recorded": all(
            (ROOT / path).is_file() for path in HISTORICAL_REFERENCE_ONLY),
    }
    return {
        "formal_file_sha256": hashes,
        "forbidden_violations": violations,
        "historical_reference_only_not_production": list(HISTORICAL_REFERENCE_ONLY),
    }, checks


def run() -> int:
    started = time.perf_counter()
    static, static_checks = _static_closure()
    convergence_metrics, convergence_checks = _step_convergence()
    regressions = _run_regressions()
    checks = {
        **static_checks, **convergence_checks,
        "all_stage_regressions_passed": all(
            item["passed"] for item in regressions.values()),
    }
    report = {
        "schema_version": "a4-c0r8-coordinate-closure-v1",
        "stage": "A4-C0R-8",
        "status": "passed_pharlap_coordinate_contract_closed"
        if all(checks.values()) else "failed",
        "case_semantics": "formal_ECEF_entry_regression_and_legacy_isolation_gate",
        "checks": checks,
        "step_convergence_metrics": convergence_metrics,
        "regressions": regressions,
        "static_closure": static,
        "elapsed_s": time.perf_counter() - started,
        "source_sha256": {"test": _sha256(Path(__file__))},
        "claim_boundary": (
            "C0R-8 closes the spherical-geocentric/ECEF production coordinate "
            "contract and isolates historical local-coordinate modules. It does "
            "not validate ionospheric reflection against PHaRLAP, O/X magnetic "
            "component conventions, full-domain P2P completeness, full 3-D PE, "
            "scattering power, or absolute received power."
        ),
    }
    _write_json(REPORT, report)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-8: {failures} fail(s); status={report['status']}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
