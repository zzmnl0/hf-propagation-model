"""Read-only acceptance check for A4-C1WR artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr"
RESULT = OUT / "a4_c1wr_attribution.json"
MANIFEST = OUT / "a4_c1wr_manifest.json"
CSV_PATH = OUT / "a4_c1wr_anchor_metrics.csv"
P2P_RESULT = OUT / "p2p_regression" / "a4_c1wr_p2p_high_low_regression.json"
C1_SHADOW = OUT / "c1_shadow" / "a4_c1wr_c1_shadow_summary.json"
FIGURES = (
    "a4_c1wr_boundary_crossing.png",
    "a4_c1wr_dr_scale_platform.png",
    "a4_c1wr_factor_attribution.png",
    "a4_c1wr_high_o_variant_trajectories.png",
    "a4_c1wr_regression_gate_matrix.png",
)
STATUSES = {
    "identified_and_repaired_grid_derivative",
    "identified_but_not_safely_repairable",
    "multiple_or_unresolved_grid_derivative_contributors",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run() -> int:
    checks = {
        "core_artifacts_exist": all(
            path.is_file() for path in (RESULT, MANIFEST, CSV_PATH)),
    }
    if not checks["core_artifacts_exist"]:
        print(" core_artifacts_exist: FAIL")
        return 1
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks.update({
        "status_enumerated": result.get("status") in STATUSES,
        "seven_variants": len(result.get("variants", [])) == 7,
        "four_cases": set(result.get("case_errors", {})) == {
            "low_O", "low_X", "high_O", "high_X"},
        "thresholds_unchanged": result.get("gates", {}).get(
            "thresholds_unchanged") is True,
        "prior_c1_failure_preserved": result.get("prior_c1_status") == (
            "failed_c1_chapman_high_low_external_comparison"),
        "prior_c1w_defect_preserved": result.get("prior_c1w_status") == (
            "implementation_defect_identified"),
        "a4_c2_not_admitted": result.get("a4_c2_admission") is False,
        "formula_units_explicit": (
            result.get("formula", {}).get("D_units") == "km^-2"
            and result.get("formula", {}).get("D_r_units") == "km^-3"),
        "production_steps_frozen": result.get("production_steps_km") == [
            0.125, 0.25, 0.5],
    })
    with CSV_PATH.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    checks["anchor_csv_complete"] = bool(rows) and {
        row["variant"] for row in rows} == set(result["variants"])
    numeric_keys = (
        "height_km", "reference_norm_km^-3",
        "reference_platform_spread_km^-3", "candidate_norm_km^-3",
        "candidate_spread_km^-3", "stability_limit_km^-3",
        "stability_ratio", "relative_reference_error",
    )
    checks["anchor_metrics_finite"] = all(
        all(np.isfinite(float(row[key])) for key in numeric_keys)
        for row in rows)
    checks["figures_valid"] = all(
        (OUT / name).is_file()
        and (lambda image: image.ndim in (2, 3)
             and image.shape[0] >= 100 and image.shape[1] >= 100
             and np.nanstd(image) > 0.01)(mpimg.imread(OUT / name))
        for name in FIGURES)
    checks["source_hashes_match"] = all(
        (ROOT / relative).is_file()
        and _sha256(ROOT / relative) == expected
        for relative, expected in manifest.get("source_sha256", {}).items())
    checks["artifact_hashes_match"] = all(
        (ROOT / relative).is_file()
        and _sha256(ROOT / relative) == expected
        for relative, expected in manifest.get("artifact_sha256", {}).items())
    checks["manifest_status_matches"] = (
        manifest.get("status") == result.get("status"))
    p2p = json.loads(P2P_RESULT.read_text(encoding="utf-8"))
    shadow = json.loads(C1_SHADOW.read_text(encoding="utf-8"))
    checks.update({
        "p2p_frozen_root_result_recorded": (
            p2p.get("status") == "failed_current_ecef_p2p_high_low"
            and p2p["gates"]["two_unique_roots"]
            and p2p["gates"]["both_ground_single_turning"]
            and p2p["gates"]["hamiltonian_residual"]
            and not p2p["gates"]["landing_error_below_5_km"]),
        "c1_shadow_expected_failure_preserved": (
            shadow.get("status") == "passed_restricted_c1_failure_preserved"
            and shadow.get("case_pass") == {
                "low_pharlap_anchor": True,
                "low_ah3d_anchor": True,
                "high_pharlap_anchor": False,
                "high_ah3d_anchor": False,
            }
            and shadow.get("original_c1_report_preserved") is True),
    })
    if result["status"] == "identified_and_repaired_grid_derivative":
        checks["repair_has_selected_variant"] = (
            result.get("selected_variant") in {
                "cubic_ne_linear_b", "linear_ne_cubic_b", "cubic_ne_cubic_b"})
        checks["repair_gates_all_pass"] = all(result.get("gates", {}).values())
        checks["production_provider_modified"] = (
            result.get("production_provider_modified") is True)
    failures = 0
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
        failures += int(not passed)
    print(f"A4-C1WR: {failures} fail(s); status={result.get('status')}")
    return failures


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
