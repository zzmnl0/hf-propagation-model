"""Read-only verification for A4-C1G frozen artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1g"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(stage: int, allowed: set[str]) -> tuple[dict, list[str]]:
    failures = []
    path = OUT / f"c1g_{stage}" / f"c1g_{stage}_freeze_manifest.json"
    if not path.is_file():
        return {}, [f"missing manifest: {path}"]
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") not in allowed:
        failures.append(f"unexpected stage {stage} status: {value.get('status')}")
    for name, passed in value.get("checks", {}).items():
        if name == "exactly_one_solver_closes_both_high":
            continue
        if passed is not True:
            failures.append(f"stage {stage} check failed: {name}")
    for item in value.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file():
            failures.append(f"missing artifact: {artifact}")
        elif sha256(artifact) != item["sha256"]:
            failures.append(f"artifact drift: {artifact}")
    return value, failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--through-stage", type=int, choices=(0, 1, 2, 6), default=6)
    args = parser.parse_args()
    failures = []
    _, stage0 = verify_manifest(0, {"passed_c1g0_frozen_admission"})
    failures.extend(stage0)
    if args.through_stage >= 1:
        manifest, stage1 = verify_manifest(1, {
            "identified_c1g1_spherical_forward_equation_candidate",
            "completed_c1g1_no_unique_spherical_candidate",
        })
        failures.extend(stage1)
        metrics = OUT / "c1g_1" / "c1g_1_spherical_benchmark_metrics.csv"
        if metrics.is_file():
            with metrics.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != 48:
                failures.append(f"expected 48 metric rows, got {len(rows)}")
        else:
            failures.append("missing C1G-1 metrics CSV")
        if manifest and "B=0 spherical Chapman" not in manifest.get("semantics", ""):
            failures.append("C1G-1 semantics missing B=0 spherical scope")
    if args.through_stage >= 2:
        _, stage2 = verify_manifest(2, {
            "identified_c1g2_pharlap_effective_no_field_medium_candidate",
            "completed_c1g2_effective_medium_not_closed",
        })
        failures.extend(stage2)
        metrics = OUT / "c1g_2" / "c1g_2_launch_effective_medium_metrics.csv"
        if metrics.is_file():
            with metrics.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != 4:
                failures.append(f"expected 4 C1G-2 rows, got {len(rows)}")
        else:
            failures.append("missing C1G-2 metrics CSV")
    if args.through_stage >= 6:
        manifest, stage6 = verify_manifest(6, {
            "identified_spherical_equation_or_invariant_defect",
            "identified_launch_or_geometry_semantics_defect",
            "pharlap_no_field_branch_difference",
            "multiple_or_unresolved_contributors",
        })
        failures.extend(stage6)
        result = OUT / "c1g_6" / "c1g_6_final_attribution.json"
        if result.is_file():
            value = json.loads(result.read_text(encoding="utf-8"))
            if value.get("formal_c1_status") != (
                    "failed_c1_chapman_high_low_external_comparison"):
                failures.append("formal C1 failure was not retained")
            if value.get("ah3d_production_repair") != (
                    "not_permitted_not_needed_by_c1g_evidence"):
                failures.append("unexpected AH3D production repair semantics")
            if set(value.get("skipped_by_early_stop", {})) != {
                    "C1G-3", "C1G-4", "C1G-5"}:
                failures.append("early-stop stage record is incomplete")
        else:
            failures.append("missing C1G-6 result")
        if manifest and "no production repair" not in manifest.get("semantics", ""):
            failures.append("C1G-6 no-repair semantics missing")
    print(f"A4-C1G read-only verification: {len(failures)} failure(s)")
    for failure in failures:
        print(f"FAIL: {failure}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
