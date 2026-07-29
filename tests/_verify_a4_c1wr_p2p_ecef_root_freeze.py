"""Read-only acceptance for the A4-C1WR-P2PF frozen result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr_p2pf"
REPORT = OUT / "a4_c1wr_p2pf_result.json"
MANIFEST = OUT / "a4_c1wr_p2pf_manifest.json"
EXPECTED_STATUS = "passed_current_ecef_chapman_p2p_root_freeze"
EXPECTED_SEMANTICS = "current_strict_spherical_ECEF_cubic_Chapman_P2P_reference"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    checks = {}
    try:
        report = json.loads(REPORT.read_text(encoding="utf-8"))
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except Exception:
        print(" report_and_manifest_readable: FAIL")
        print("A4-C1WR-P2PF: 1 fail(s)")
        return 1

    checks["status_passed"] = (
        report.get("status") == EXPECTED_STATUS
        and manifest.get("status") == EXPECTED_STATUS)
    checks["semantics_exact"] = (
        report.get("semantics") == EXPECTED_SEMANTICS
        and manifest.get("semantics") == EXPECTED_SEMANTICS)
    configuration = report.get("configuration", {})
    checks["fixed_configuration"] = (
        configuration.get("earth_radius_km") == 6371.0
        and configuration.get("frequency_MHz") == 10.0
        and configuration.get("B") == 0
        and configuration.get("collision") == 0
        and configuration.get("one_way") is True
        and configuration.get("bearing_deg") == 0.0
        and configuration.get("provider_interpolation") == "cubic"
        and configuration.get("spatial_fd_steps_km") == [0.125, 0.25, 0.5]
        and configuration.get("source_trace_grid_shape") == [149, 17, 464]
        and configuration.get("factorized_trace_grid_shape") == [4, 4, 464]
        and configuration.get("horizontal_factorization")
        == "exact_after_bitwise_horizontal_uniformity_check")
    checks["frozen_brackets"] = configuration.get("brackets_deg") == {
        "low_ray": [21.9, 22.0], "high_ray": [39.5, 39.6]}
    checks["legacy_semantics_preserved"] = report.get(
        "legacy_coordinate_reference_only_deg") == {
            "low_ray": 35.33980751304579,
            "high_ray": 39.918776604143765}
    rows = report.get("rows", [])
    checks["expected_root_rows"] = (
        len(rows) == 6
        and {(row["branch"], row["profile"], row["solver"]) for row in rows}
        == {
            (branch, profile, "bounded_Brent")
            for branch in ("low_ray", "high_ray")
            for profile in ("research", "reference")
        } | {
            (branch, "reference", "bounded_TRF_2D")
            for branch in ("low_ray", "high_ray")
        })
    checks["all_metrics_finite"] = bool(rows) and all(
        np.all(np.isfinite([
            row["elevation_deg"], row["rx_ecef_error_km"],
            row["rx_enu_error_km"], row["apex_km"], row["geom_path_km"],
            row["group_path_km"], row["tau_ms"],
            row["hamiltonian_residual_max"],
            row["spatial_fd_stability_ratio_max"], row["elapsed_s"],
        ])) for row in rows)
    gates = report.get("gates", {})
    checks["all_acceptance_gates_passed"] = bool(gates) and all(gates.values())
    checks["provider_hash_versioned"] = (
        len(configuration.get("provider_hash", "")) == 64
        and configuration.get("provider_hash")
        != configuration.get("prior_linear_provider_hash"))

    source_hashes = manifest.get("source_sha256", {})
    source_ok = bool(source_hashes)
    for name, expected in source_hashes.items():
        path = ROOT / name
        source_ok = source_ok and path.is_file() and _sha256(path) == expected
    checks["source_hashes_match"] = source_ok
    artifact_hashes = manifest.get("artifact_sha256", {})
    artifact_ok = bool(artifact_hashes)
    for name, expected in artifact_hashes.items():
        path = ROOT / name
        artifact_ok = artifact_ok and path.is_file() and _sha256(path) == expected
    checks["artifact_hashes_match"] = artifact_ok

    figure_ok = len(report.get("figures", [])) == 4
    for item in report.get("figures", []):
        path = ROOT / item["path"]
        try:
            image = mpimg.imread(path)
            figure_ok = (
                figure_ok and path.is_file()
                and _sha256(path) == item["sha256"]
                and image.shape[0] >= 900 and image.shape[1] >= 1400
                and np.all(np.isfinite(image))
                and float(np.std(image[..., :3])) >= 1.0e-3)
        except Exception:
            figure_ok = False
    checks["figures_valid"] = figure_ok

    npz_path = OUT / "a4_c1wr_p2pf_trajectories.npz"
    try:
        with np.load(npz_path, allow_pickle=False) as data:
            checks["trajectory_npz_complete"] = (
                len(data.files) == 16
                and all(np.asarray(data[name]).size > 1
                        and np.all(np.isfinite(data[name]))
                        for name in data.files))
    except Exception:
        checks["trajectory_npz_complete"] = False

    failures = [name for name, passed in checks.items() if not passed]
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    print(
        f"A4-C1WR-P2PF: {len(failures)} fail(s); "
        f"status={report.get('status')}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
