"""Strict no-recompute verification for the A4-C1-4 failure freeze."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import _image_check, _sha256


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "final"
MANIFEST = OUT / "a4_c1_freeze_manifest.json"
VALIDATION = OUT / "a4_c1_validation.json"
EXPECTED = "failed_c1_chapman_high_low_external_comparison"


def main() -> int:
    failures = []
    if not MANIFEST.is_file():
        print("A4-C1-4 verify: FAIL manifest missing")
        return 1
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    validation = json.loads(VALIDATION.read_text(encoding="utf-8"))
    if manifest.get("status") != EXPECTED or validation.get("status") != EXPECTED:
        failures.append("unexpected final status")
    if validation.get("c1_passed") is not False:
        failures.append("C1 failure was not preserved")
    if validation.get("c2_admitted") is not False:
        failures.append("C2 was incorrectly admitted")
    failed = validation.get("formal_failed_gates", [])
    if len(failed) != 2 or any(
        item.get("branch_id") != "high_ray"
        or item.get("gate") != "trajectory_rms_3d_km"
        for item in failed
    ):
        failures.append("formal failure signature drift")
    if not all(validation.get("checks", {}).values()):
        failures.append("completion check failed")

    for group in ("source_files", "frozen_artifacts"):
        for item in manifest.get(group, []):
            path = Path(item["path"])
            if not path.is_file() or _sha256(path) != item["sha256"]:
                failures.append(f"{group} drift: {path}")

    images = [
        Path(item["path"]) for item in manifest.get("frozen_artifacts", [])
        if Path(item["path"]).suffix.lower() == ".png"
    ]
    if len(images) != 3:
        failures.append("expected exactly three C1-4 visualizations")
    for path in images:
        if path.is_file() and not _image_check(path)["passed"]:
            failures.append(f"invalid visualization: {path}")

    report = ROOT / "word" / "a4_c1_chapman_external_validation_report.txt"
    text = report.read_text(encoding="utf-8") if report.is_file() else ""
    for required in (
        EXPECTED,
        "A4-C1最终状态为failed",
        "不允许进入A4-C2正式O/X验收",
        "不能判断AH3D或PHaRLAP哪一个更接近真实电离层",
    ):
        if required not in text:
            failures.append(f"report statement missing: {required}")

    print(json.dumps({
        "status": validation.get("status"),
        "source_files": len(manifest.get("source_files", [])),
        "frozen_artifacts": len(manifest.get("frozen_artifacts", [])),
        "visualizations": len(images),
        "failures": failures,
    }, indent=2, ensure_ascii=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
