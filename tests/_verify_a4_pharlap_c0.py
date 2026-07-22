"""Minimal machine-readable acceptance check for A4-C0."""

import json
from pathlib import Path


REPORT = Path(__file__).resolve().parents[1] / "output2" / "a4_pharlap" / "a4_c0_validation.json"


def run() -> int:
    if not REPORT.exists():
        print("A4-C0: FAIL missing validation report")
        return 1
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    checks = {
        "status_passed": report.get("status") == "passed",
        "c0_semantics": report.get("case_semantics") == "interface_and_coordinate_precheck_not_physics_validation",
        "no_external_validation_claim": report.get("external_validation_claim") == "none_C0_only",
        "all_checks_passed": all(report.get("checks", {}).values()),
        "hashes_present": all(report.get("software", {}).get(key) for key in (
            "pharlap_mex_sha256", "matlab_script_sha256")),
        "three_figures_verified": (
            set(report.get("figures", {})) == {"pharlap", "ah3d", "comparison"}
            and all(all(item.get("checks", {}).values())
                    for item in report.get("figures", {}).values())),
    }
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
