"""Fast structural verification of the bounded A3Rnew-P5new decision."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "output2/a3rnew_p5/a3rnew_p5_acceptance.json"
REPORT = ROOT / "word/a3rnew_p5_accuracy_acceptance_report.txt"


def main() -> int:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    checks = {
        "limited_gate_ready": result["status"] == "a4_limited_scope_ready",
        "all_acceptance_checks": all(result["acceptance_checks"].values()),
        "all_evidence_integrity": all(result["evidence_integrity"].values()),
        "all_focused_checks": all(
            item["status"] == "passed" for item in result["focused_checks"].values()),
        "bounded_completeness": result["statuses"]["completeness_status"]
        == "reference_case_limited_pass",
        "external_not_claimed": result["statuses"]["external_validation_status"]
        == "not_run_until_A4",
        "p2_6h_deferred": result["not_run"]["p2_6h"]
        == "deferred_full_completeness_audit",
        "equations_unchanged": result["equation_policy"]
        == "production equations and numerical tolerances unchanged",
        "report_exists": REPORT.exists(),
    }
    for name, passed in checks.items():
        print("{}: {}".format(name, "OK" if passed else "FAIL"))
    failures = sum(not passed for passed in checks.values())
    print("A3Rnew-P5new acceptance: {} fail(s)".format(failures))
    return failures


if __name__ == "__main__":
    raise SystemExit(main())
