"""Versioned rerun of the original B=0 C1-0/C1-1 restricted comparison."""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests import _run_a4_c1_validation as c1  # noqa: E402


def main() -> int:
    original_report = (
        ROOT / "output2" / "a4_pharlap" / "c1" / "same_angle"
        / "c1_1_validation.json")
    original_hash = c1._sha256(original_report)
    base = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr" / "c1_shadow"
    c1.BASE = base
    c1.ADMISSION = base / "admission"
    c1.SAME = base / "same_angle"
    c1.RAW = c1.SAME / "raw"
    c1.ADMISSION_REPORT = c1.ADMISSION / "c1_0_admission.json"
    c1.ADMISSION_FIGURE = c1.ADMISSION / "c1_0_source_status.png"
    c1.ADMISSION_MANIFEST = c1.ADMISSION / "c1_0_freeze_manifest.json"
    c1.PH_SUMMARY = c1.RAW / "c1_same_angle_pharlap_summary.json"
    c1.PH_LOG = c1.RAW / "c1_same_angle_pharlap.log"
    c1.METRICS_CSV = c1.SAME / "c1_same_angle_metrics.csv"
    c1.METRICS_JSON = c1.SAME / "c1_same_angle_metrics.json"
    c1.TRAJECTORY_FIGURE = c1.SAME / "c1_same_angle_trajectories.png"
    c1.RESIDUAL_FIGURE = c1.SAME / "c1_same_angle_residuals.png"
    c1.METRICS_FIGURE = c1.SAME / "c1_same_angle_metrics.png"
    c1.SAME_REPORT = c1.SAME / "c1_1_validation.json"
    c1.SAME_MANIFEST = c1.SAME / "c1_1_freeze_manifest.json"

    input_path, _ = c1._refined_input()
    sources = c1._current_sources(input_path)
    checks = {
        "refined_input_hash_unchanged": (
            c1._sha256(input_path) == c1.EXPECTED_INPUT_HASH),
        "pharlap_mex_available": c1.MEX.is_file(),
        "original_c1_report_frozen": original_report.is_file(),
        "versioned_output_separate": c1.SAME != original_report.parent,
    }
    source_tag = hashlib.sha256(json.dumps(
        sources, sort_keys=True).encode("ascii")).hexdigest()[:12]
    c1.ADMISSION.mkdir(parents=True, exist_ok=True)
    visualization = c1._plot_admission(checks, source_tag)
    admission_status = (
        "passed_c1_0_frozen_admission"
        if all(checks.values()) and visualization["passed"]
        else "blocked_c1_0")
    c1._atomic_json(c1.ADMISSION_REPORT, {
        "schema_version": "a4-c1wr-c1-shadow-admission-v1",
        "status": admission_status,
        "checks": checks,
        "source_sha256": sources,
        "visualization": visualization,
        "semantics": "read_only_frozen_input_admission_without_prep_rerun",
    })
    c1._atomic_json(c1.ADMISSION_MANIFEST, {
        "schema_version": "a4-c1wr-c1-shadow-admission-manifest-v1",
        "status": admission_status,
        "source_sha256": sources,
        "frozen_artifacts": [
            {"path": str(path), "sha256": c1._sha256(path)}
            for path in (c1.ADMISSION_REPORT, c1.ADMISSION_FIGURE)
        ],
    })
    if not admission_status.startswith("passed"):
        return 2

    result = c1.run_same_angle(workers=2)
    report = json.loads(c1.SAME_REPORT.read_text(encoding="utf-8"))
    original_preserved = c1._sha256(original_report) == original_hash
    expected_failure = report.get("status") == "failed_c1_1_numeric_thresholds"
    summary = {
        "schema_version": "a4-c1wr-c1-shadow-summary-v1",
        "status": (
            "passed_restricted_c1_failure_preserved"
            if result == 1 and expected_failure and original_preserved
            else "failed_restricted_c1_shadow"),
        "same_angle_status": report.get("status"),
        "case_pass": {
            item["case_id"]: bool(item["passed"])
            for item in report.get("metrics", [])
        },
        "original_c1_report_sha256": original_hash,
        "original_c1_report_preserved": original_preserved,
        "a4_c2_admission": False,
    }
    c1._atomic_json(base / "a4_c1wr_c1_shadow_summary.json", summary)
    return 0 if summary["status"].startswith("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
