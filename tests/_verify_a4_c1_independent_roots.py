"""Verify frozen A4-C1-2 diagnostic artifacts without ray tracing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "independent_roots"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((OUT / "c1_2_freeze_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "completed_c1_2_diagnostic_after_same_angle_failure"
    for item in manifest["frozen_artifacts"]:
        path = Path(item["path"])
        assert path.is_file() and sha256(path) == item["sha256"], path
    report = json.loads((OUT / "c1_2_validation.json").read_text(encoding="utf-8"))
    assert all(report["checks"].values()) and len(report["metrics"]) == 2
    assert all(item["numeric_gates_passed"] for item in report["metrics"])
    values = {item["branch_id"]: item for item in report["metrics"]}
    assert values["low_ray"]["launch_elevation_absolute_error_deg"] <= 0.1
    assert values["high_ray"]["launch_elevation_absolute_error_deg"] <= 0.1
    assert report["attribution_indicators"]["high_ray"]["same_angle_max_trajectory_rms_km"] > 2.0
    assert values["high_ray"]["trajectory_rms_3d_km"] < 2.0
    assert report["checks"]["does_not_override_c1_1_failure"]
    print(f"A4-C1-2 frozen verification: OK status={report['status']}")


if __name__ == "__main__":
    main()
