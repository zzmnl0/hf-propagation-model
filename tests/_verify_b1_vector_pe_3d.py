"""Read-only verification of frozen B1 vector PE artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "b1_vector_pe_3d"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    failures = []
    metrics_path = OUT / "b1_metrics.json"
    manifest_path = OUT / "b1_manifest.json"
    if not metrics_path.is_file() or not manifest_path.is_file():
        print("FAIL missing B1 metrics or manifest")
        return 1
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    allowed = {
        "passed_b1_vector_wide_angle_pe_reference",
        "failed_b1_vector_pe_acceptance",
        "blocked_by_cpu_memory_or_solver_stability",
    }
    if metrics.get("final_status") not in allowed:
        failures.append("invalid final_status")
    if manifest.get("status") != metrics.get("final_status"):
        failures.append("manifest status mismatch")

    semantics = metrics.get("semantics", {})
    expected_semantics = {
        "solver_dimensionality": "full_3d_two_transverse_vector_pe",
        "propagation_direction": "one_way_forward",
        "validity": "research_reference_only",
        "calibration_status": "reference_uncalibrated",
        "external_c1_status": "failed_c1_chapman_high_low_external_comparison",
        "external_pharlap_status": (
            "external_pharlap_comparison_paused_due_effective_medium_inconsistency"
        ),
    }
    for key, expected in expected_semantics.items():
        if semantics.get(key) != expected:
            failures.append(f"semantic mismatch: {key}")

    for group in ("sources", "artifacts"):
        for relative, expected_hash in manifest.get(group, {}).items():
            path = ROOT / relative
            if not path.is_file():
                failures.append(f"missing {group[:-1]}: {relative}")
            elif _sha256(path) != expected_hash:
                failures.append(f"hash mismatch: {relative}")

    required_images = [
        "b1_tensor_eigen_acceptance.png",
        "b1_wide_angle_dispersion.png",
        "b1_bishop_frame_jacobian.png",
        "b1_gaussian_field_and_flux.png",
        "b1_wkb_mode_acceptance.png",
        "b1_pml_and_energy_budget.png",
        "b1_acceptance_matrix.png",
    ]
    for name in required_images:
        path = OUT / name
        if not path.is_file() or path.stat().st_size < 10_000:
            failures.append(f"invalid image: {name}")
            continue
        with Image.open(path) as image:
            if image.width < 600 or image.height < 300:
                failures.append(f"image dimensions too small: {name}")
            image.verify()

    gates = metrics.get("gates", {})
    required_nonresource = [
        "tensor",
        "wide_angle_30deg",
        "bishop_geometry",
        "uniform_mode_leakage",
        "light_gaussian",
        "step_halving",
        "light_wkb",
        "planned_w0_6lambda_25lambda_full_grid",
        "planned_grid_lambda4_lambda8",
        "planned_domain_expansion_50pct",
        "planned_pml_thickness_change",
        "planned_three_frequency_group_delay",
        "measured_pml_reflection_minus40dB",
        "collision_energy_budget",
    ]
    nonresource_pass = all(gates.get(name) is True for name in required_nonresource)
    refined_pass = gates.get(
        "planned_w0_6lambda_25lambda_lambda8_refined"
    ) is True
    resource_blocked = (
        metrics.get("resource_assessment", {}).get("resource_blocked") is True
    )
    status = metrics.get("final_status")
    if status == "passed_b1_vector_wide_angle_pe_reference":
        if not nonresource_pass or not refined_pass:
            failures.append("passed status without every required gate")
    elif status == "blocked_by_cpu_memory_or_solver_stability":
        if not nonresource_pass or refined_pass or not resource_blocked:
            failures.append("blocked status is inconsistent with gate matrix")

    for entry in metrics.get("tensor", []):
        if not float(entry["value"]) <= float(entry["limit"]):
            failures.append(f"tensor gate failed: {entry['metric']}")
    for entry in metrics.get("wide_angle", []):
        if not float(entry["relative_k_s_error"]) <= float(entry["limit"]):
            failures.append(f"wide-angle gate failed: {entry['angle_deg']}")

    production_files = [
        ROOT / "main.py",
        ROOT / "main2.py",
        ROOT / "models" / "unified_hybrid_pipeline.py",
    ]
    for path in production_files:
        if "vector_pe_3d" in path.read_text(encoding="utf-8"):
            failures.append(f"B1 imported by production entry: {path.name}")

    print(f"B1 artifact status: {status}")
    for failure in failures:
        print(f"FAIL {failure}")
    print(f"{len(failures)} fail(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
