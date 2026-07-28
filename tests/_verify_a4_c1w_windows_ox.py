"""Read-only acceptance check for A4-C1W Windows O/X validation artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1w"
RESULT = OUT / "a4_c1w_validation.json"
MANIFEST = OUT / "a4_c1w_manifest.json"
CSV_PATH = OUT / "a4_c1w_ox_metrics.csv"
FIGURES = (
    "a4_c1w_ah3d_capability_matrix.png",
    "a4_c1w_igrf3d_field.png",
    "a4_c1w_wgs84_ox_equation_audit.png",
    "a4_c1w_ah3d_pharlap_ox_xyz.png",
    "a4_c1w_ox_acceptance_summary.png",
    "a4_c1w_pharlap_traceability.png",
)
STATUSES = {
    "complete_collisionless_ox_hamiltonian_validated",
    "implementation_defect_identified",
    "pharlap_behavior_mismatch_unresolved",
    "blocked_by_model_validity_boundary",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run() -> int:
    checks = {}
    checks["core_artifacts_exist"] = all(
        path.is_file() for path in (RESULT, MANIFEST, CSV_PATH))
    if not checks["core_artifacts_exist"]:
        for name, passed in checks.items():
            print(f" {name}: {'OK' if passed else 'FAIL'}")
        return 1
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks["status_enumerated"] = result.get("status") in STATUSES
    checks["scope_explicit"] = result.get("scope") == (
        "cold_collisionless_geometric_optics_WGS84_position_dependent_IGRF_OX")
    checks["prior_c1_failure_preserved"] = result.get("prior_c1_status") == (
        "failed_c1_chapman_high_low_external_comparison")
    checks["a4_c2_not_auto_admitted"] = result.get("a4_c2_admission") is False
    trace = result.get("pharlap_traceability", {})
    checks["pharlap_classification"] = trace.get("classification") == (
        "behavior_traceable_core_source_unavailable")
    checks["pharlap_core_not_claimed_open"] = (
        trace.get("fortran_source_files_present") is False)
    checks["local_igrf_absorption"] = (
        result.get("absorption", {}).get("geomagnetic_field_semantics")
        == "position_dependent_3d_provider")
    metrics = result.get("pharlap_behavior_metrics", {})
    checks["four_ox_cases"] = set(metrics) == {
        "low_O", "low_X", "high_O", "high_X"}
    checks["finite_metrics"] = all(
        np.all(np.isfinite(list(row.values()))) for row in metrics.values())
    comparable = {
        case_id for case_id, row in metrics.items()
        if row.get("comparison_valid") is True}
    checks["ground_ground_comparisons_explicit"] = comparable == {
        "low_O", "low_X", "high_X"}
    checks["high_o_not_misreported_as_landing_error"] = (
        metrics.get("high_O", {}).get("comparison_valid") is False
        and metrics.get("high_O", {}).get("pharlap_ray_label") == -2)
    checks["defect_has_numerical_rejection_evidence"] = (
        result.get("status") != "implementation_defect_identified"
        or result.get("numerical_domain_rejections") == ["high_O"])
    checks["defect_not_misclassified_as_model_boundary"] = (
        result.get("model_validity_boundary_encountered") is False)
    with CSV_PATH.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    checks["csv_four_cases"] = len(rows) == 4
    checks["figures_valid"] = all(
        (OUT / name).is_file()
        and (lambda image: image.ndim in (2, 3)
             and image.shape[0] >= 100 and image.shape[1] >= 100
             and np.nanstd(image) > 0.01)(mpimg.imread(OUT / name))
        for name in FIGURES)
    hashes = manifest.get("sha256", {})
    checks["manifest_hashes_match"] = all(
        (ROOT / relative).is_file()
        and _sha256(ROOT / relative) == expected
        for relative, expected in hashes.items())
    checks["manifest_status_matches"] = (
        manifest.get("status") == result.get("status"))
    failures = 0
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
        failures += int(not passed)
    print(f"A4-C1W: {failures} fail(s); status={result.get('status')}")
    return failures


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
