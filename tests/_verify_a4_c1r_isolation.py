"""Strict no-recompute verification for A4-C1R artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import _image_check, _sha256


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1r"
ALLOWED = {
    0: {"passed_c1r0_failure_signature_frozen"},
    1: {"passed_c1r1_integration_not_primary"},
    2: {"passed_c1r2_frozen_field_not_primary",
        "completed_c1r2_field_effect_unresolved"},
    3: {"passed_c1r3_spatial_derivative_not_primary",
        "completed_c1r3_spatial_derivative_candidate",
        "completed_c1r3_spatial_derivative_unresolved"},
    4: {"completed_c1r4_remaining_implementation_attribution"},
    5: {"completed_c1r5_ah3d_defect_candidate_requires_repair_plan",
        "completed_c1r5_software_difference_or_unresolved"},
}


def main() -> int:
    failures = []
    statuses = {}
    for stage in range(6):
        manifest_path = OUT / f"c1r_{stage}" / f"c1r_{stage}_freeze_manifest.json"
        if not manifest_path.is_file():
            failures.append(f"C1R-{stage}: manifest missing")
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        statuses[stage] = manifest.get("status")
        if statuses[stage] not in ALLOWED[stage]:
            failures.append(f"C1R-{stage}: unexpected status {statuses[stage]}")
        if not all(manifest.get("checks", {}).values()):
            failures.append(f"C1R-{stage}: machine check failed")
        png_count = 0
        for item in manifest.get("frozen_artifacts", []):
            path = Path(item["path"])
            if not path.is_file() or _sha256(path) != item["sha256"]:
                failures.append(f"C1R-{stage}: artifact drift {path}")
                continue
            if path.suffix.lower() == ".png":
                png_count += 1
                if not _image_check(path)["passed"]:
                    failures.append(f"C1R-{stage}: invalid PNG {path}")
        if png_count < 1:
            failures.append(f"C1R-{stage}: no frozen PNG")

    c1 = json.loads(
        (ROOT / "output2" / "a4_pharlap" / "c1" / "same_angle"
         / "c1_1_validation.json").read_text(encoding="utf-8"))
    if c1.get("status") != "failed_c1_1_numeric_thresholds":
        failures.append("formal C1-1 status was overwritten")
    print(json.dumps({"statuses": statuses, "failures": failures},
                     indent=2, ensure_ascii=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
