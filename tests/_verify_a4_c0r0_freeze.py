"""Minimal acceptance check for the immutable A4-C0R-0 snapshot."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "output2" / "a4_pharlap" / "c0r0_freeze" / "c0r0_freeze_manifest.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> int:
    if not MANIFEST.is_file():
        print("A4-C0R-0: FAIL missing manifest")
        return 1
    report = json.loads(MANIFEST.read_text(encoding="utf-8"))
    c0 = report.get("frozen_c0_artifacts", [])
    legacy = report.get("legacy_coordinate_artifacts", [])
    inventory = report.get("migration_inventory", {})
    all_artifacts = c0 + legacy
    checks = {
        "status": report.get("status") == "passed_legacy_coordinate_artifacts_frozen",
        "c0_evidence_complete": len(c0) == 14,
        "legacy_marked": bool(legacy) and all(
            item.get("coordinate_semantics") == "legacy_coordinate_reference_only"
            for item in legacy
        ),
        "frozen_hashes_match": all(
            Path(item["frozen_path"]).is_file()
            and _sha256(Path(item["frozen_path"])) == item["sha256"]
            for item in all_artifacts
        ),
        "c0_sources_not_overwritten": all(
            Path(item["source_path"]).is_file()
            and _sha256(Path(item["source_path"])) == item["sha256"]
            for item in c0
        ),
        "inventory_tokens": set(inventory.get("tokens", ())) == {
            "trajectory_3d", "earth_flat", "TX_POS", "RX_POS", "z_cart",
            "legacy_parabolic",
        },
        "production_callers_recorded": inventory.get("caller_count", 0) > 0,
        "partial_migration_forbidden": inventory.get("partial_migration_allowed") is False,
        "entrypoint_switch_blocked": (
            report.get("forbidden_entrypoint_switch_until_inventory_closed") is True
        ),
    }
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-0: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
