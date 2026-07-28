"""Strict no-recompute verification for A4-C1P artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import _image_check, _sha256


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1p"
ALLOWED = {
    0: {"passed_c1p0_frozen_admission"},
    1: {"completed_c1p1_pharlap_convergence_scan"},
    2: {"completed_c1p2_state_semantics_audit"},
    3: {"completed_c1p3_local_rhs_audit"},
    4: {
        "completed_c1p4_pharlap_integration_primary_candidate",
        "completed_c1p4_pharlap_integration_contributor",
        "completed_c1p4_state_or_rhs_semantics_candidate",
        "completed_c1p4_software_difference_or_unresolved",
    },
}
FINAL_C1_STATUS = "failed_c1_chapman_high_low_external_comparison"


def main() -> int:
    failures = []
    statuses = {}
    for stage in range(5):
        manifest_path = OUT / f"c1p_{stage}" / f"c1p_{stage}_freeze_manifest.json"
        if not manifest_path.is_file():
            failures.append(f"C1P-{stage}: manifest missing")
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        statuses[stage] = manifest.get("status")
        if statuses[stage] not in ALLOWED[stage]:
            failures.append(f"C1P-{stage}: unexpected status {statuses[stage]}")
        if not all(manifest.get("checks", {}).values()):
            failures.append(f"C1P-{stage}: machine check failed")
        png_count = 0
        for item in manifest.get("frozen_artifacts", []):
            path = Path(item["path"])
            if not path.is_file() or _sha256(path) != item["sha256"]:
                failures.append(f"C1P-{stage}: artifact drift {path}")
                continue
            if path.suffix.lower() == ".png":
                png_count += 1
                if not _image_check(path)["passed"]:
                    failures.append(f"C1P-{stage}: invalid PNG {path}")
        if png_count < 1:
            failures.append(f"C1P-{stage}: no frozen PNG")

    final_c1 = json.loads(
        (ROOT / "output2" / "a4_pharlap" / "c1" / "final"
         / "a4_c1_validation.json").read_text(encoding="utf-8"))
    if final_c1.get("status") != FINAL_C1_STATUS:
        failures.append("formal C1 status was overwritten")
    final = json.loads(
        (OUT / "c1p_4" / "c1p_4_attribution.json").read_text(encoding="utf-8"))
    if final.get("c2_gate") != "manual_decision_required_formal_c1_still_failed":
        failures.append("C2 manual decision gate missing")
    print(json.dumps({"statuses": statuses, "failures": failures},
                     indent=2, ensure_ascii=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
