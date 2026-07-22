"""Fast structural verification of the frozen A3Rnew-P3 result."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "output2" / "a3rnew_p3" / "a3rnew_p3_validation.json"


def main() -> int:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    failures = []

    def require(name, condition):
        if not condition:
            failures.append(name)

    require("stage_passed", result["status"] == "passed")
    require("all_checks", all(result["checks"].values()))
    require("coleman_not_strict",
            result["coleman_2011"]["reproduction_semantics"] ==
            "literature_guided_parameterized_reproduction"
            and result["coleman_2011"]["strict_reproduction"] is False)
    require("coleman_two_roots", len(result["coleman_2011"]["reference_roots_deg"]) == 2)
    require("nosikov_partial",
            result["nosikov_2020"]["reproduction_semantics"] ==
            "partial_reproduction_with_assumptions"
            and result["nosikov_2020"]["comparison_status"] == "partial_not_matched")
    require("analytical_limits", all(result["analytical_limits"]["checks"].values()))
    require("frequency_points",
            result["ray_frequency_audit"]["frequencies_MHz"] == [8.0, 9.5, 11.0, 12.5, 14.0]
            and len(result["ray_frequency_audit"]["records"]) == 10)
    require("a4_not_claimed", result["not_external_validation"] == "PHaRLAP comparison remains A4")
    require("p2_6h_deferred", any("P2-6H" in risk for risk in result["residual_risks"]))
    require("plots", all((RESULT.parent / name).exists() for name in (
        "coleman_high_low_trajectories.png", "ah3d_ox_frequency_trend.png")))

    for failure in failures:
        print(f" FAIL {failure}")
    print(f"A3Rnew-P3 validation: {len(failures)} fail(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
