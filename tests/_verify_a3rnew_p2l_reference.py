"""Verify the bounded semantics and pass state of the P2-6L report."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "output2" / "a3rnew_p2l" / "p2_6l_reference_acceptance.json"


def main() -> int:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    checks = {
        "limited_pass": report["status"] == "reference_case_limited_pass",
        "bounded_semantics": "not complex-scene mode completeness" in report["semantics"],
        "full_domain_lightweight": report["configuration"]["label"] == "lightweight_full_domain_preacceptance",
        "p0c_recalled": all(item["matched"] for item in report["p0c_matches"]),
        "strict_formal_modes": len(report["p0c_modes"]) == 2 and all(
            mode["p2p_jacobian_effective"] == "central_reference_multiscale"
            and mode["branch_status"] == "ground"
            and mode["convergence_status"] == "passed"
            for mode in report["p0c_modes"]),
        "complex_fields_deferred": report["not_run"]["tid_mode_completeness"] == "not_run"
        and report["not_run"]["high_density_frozen_reference_recall"] == "deferred_as_P2-6H",
    }
    for name, passed in checks.items():
        print("{}: {}".format(name, "OK" if passed else "FAIL"))
    failures = sum(not passed for passed in checks.values())
    print("P2-6L reference acceptance: {} fail(s)".format(failures))
    return failures


if __name__ == "__main__":
    raise SystemExit(main())
