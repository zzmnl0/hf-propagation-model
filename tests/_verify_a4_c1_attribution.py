"""Verify frozen A4-C1-3 attribution artifacts without ray tracing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "attribution"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((OUT / "c1_3_freeze_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed_pending_c1r"
    for item in manifest["frozen_artifacts"]:
        path = Path(item["path"])
        assert path.is_file() and sha256(path) == item["sha256"], path
    report = json.loads((OUT / "c1_3_validation.json").read_text(encoding="utf-8"))
    assert all(report["checks"].values())
    assert len(report["failed_formal_gates"]) == 2
    assert all(item["branch_id"] == "high_ray"
               and item["gate"] == "trajectory_rms_3d_km"
               for item in report["failed_formal_gates"])
    assert report["primary_attribution"] == "numerical_difference"
    assert report["subcause"] == "unresolved_interpolation_derivative_integration_or_implementation"
    categories = {item["category"]: item["status"] for item in report["categories"]}
    assert categories["input_difference"] == "excluded"
    assert categories["coordinate_earth_model_difference"] == "excluded"
    assert categories["event_mapping_difference"] == "excluded_as_primary"
    assert categories["numerical_difference"] == "supported_primary"
    assert report["evidence"]["high_same_angle_path_residual"][
        "high_pharlap_anchor"]["residual_at_u_0p50_km"] > 1.0
    print(f"A4-C1-3 frozen verification: OK status={report['status']}")


if __name__ == "__main__":
    main()
