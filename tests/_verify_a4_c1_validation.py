"""Verify frozen A4-C1-0/1 artifacts without rerunning ray tracing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "output2" / "a4_pharlap" / "c1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(path: Path, statuses: set[str]) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    assert manifest["status"] in statuses, (path, manifest["status"])
    for item in manifest["frozen_artifacts"]:
        artifact = Path(item["path"])
        assert artifact.is_file() and sha256(artifact) == item["sha256"], artifact
    return manifest


def main() -> None:
    verify_manifest(BASE / "admission" / "c1_0_freeze_manifest.json",
                    {"passed_c1_0_frozen_admission"})
    verify_manifest(BASE / "same_angle" / "c1_1_freeze_manifest.json", {
        "passed_c1_1_same_angle_forward_comparison",
        "failed_c1_1_numeric_thresholds",
    })
    report = json.loads((BASE / "same_angle" / "c1_1_validation.json").read_text(
        encoding="utf-8"))
    assert len(report["metrics"]) == 4
    assert report["checks"]["all_repeatability"]
    assert report["checks"]["all_semantic_ground_single_turning"]
    assert report["checks"]["all_visualizations"]
    if report["status"] == "failed_c1_1_numeric_thresholds":
        failures = {(item["branch_id"], item["anchor_solver"]):
                    [key for key, passed in item["gates"].items() if not passed]
                    for item in report["metrics"] if not item["passed"]}
        assert failures == {
            ("high_ray", "PHaRLAP"): ["trajectory_rms_3d_km"],
            ("high_ray", "AH3D"): ["trajectory_rms_3d_km"],
        }
    print(f"A4-C1-0/1 frozen verification: OK status={report['status']}")


if __name__ == "__main__":
    main()
