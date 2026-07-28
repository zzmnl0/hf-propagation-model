"""Verify frozen A4-C1Q artifacts without recomputing any ray."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1q"
FINAL_C1 = (ROOT / "output2" / "a4_pharlap" / "c1" / "final"
            / "a4_c1_validation.json")
ALLOWED_FINAL = {
    "identified_and_repaired_formula_or_implementation",
    "identified_but_not_safely_repairable",
    "pharlap_effective_dynamics_not_recoverable",
    "multiple_unresolved_contributors",
}
ANGLES = np.array([
    21.94290326504532, 21.94562872092976,
    39.54254824570964, 39.58701138438118,
])
SOURCE_PATHS = {
    "runner": ROOT / "tests" / "_run_a4_c1q_attribution.py",
    "stage_runner": ROOT / "tests" / "_run_a4_c1q_attribution.py",
    "plan": ROOT / "word"
            / "a4_c1q_pharlap_truth_along_track_attribution_and_repair_plan.txt",
    "matlab_source": ROOT / "tests" / "matlab"
                     / "a4_c1q_pharlap_state_export.m",
    "refined_input": ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep1"
                     / "inputs" / "a4_c1_chapman_refined_1km.mat",
    "pharlap_mex": Path(r"D:\pharlap\pharlap_4.7.4\mex\raytrace_3d_sp.mexw64"),
    "ah_solver": ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
    "ah_dispersion": ROOT / "models" / "dispersion_ah.py",
    "ah_provider": ROOT / "models" / "ne_provider_geocentric.py",
    "config": ROOT / "config.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    failures = []
    manifests = []
    for stage in range(7):
        path = OUT / f"c1q_{stage}" / f"c1q_{stage}_freeze_manifest.json"
        if not path.is_file():
            failures.append(f"missing manifest: {path}")
            continue
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifests.append(manifest)
        for key, expected_hash in manifest.get("source_sha256", {}).items():
            source_path = SOURCE_PATHS.get(key)
            if source_path is not None and (
                    not source_path.is_file() or sha256(source_path) != expected_hash):
                failures.append(f"source drift: C1Q-{stage} {key}")
            if key.startswith("c1q") and key.endswith("_manifest"):
                upstream = int(key[3:key.index("_")])
                upstream_path = (OUT / f"c1q_{upstream}"
                                 / f"c1q_{upstream}_freeze_manifest.json")
                if not upstream_path.is_file() or sha256(upstream_path) != expected_hash:
                    failures.append(f"manifest chain drift: C1Q-{stage} {key}")
        for item in manifest.get("frozen_artifacts", []):
            artifact = Path(item["path"])
            if not artifact.is_file() or sha256(artifact) != item["sha256"]:
                failures.append(f"artifact drift: {artifact}")
            elif artifact.suffix.lower() == ".png":
                pixels = np.asarray(mpimg.imread(artifact), dtype=float)
                if (pixels.size == 0 or not np.isfinite(pixels).all()
                        or min(pixels.shape[:2]) < 200 or np.var(pixels) <= 1e-8):
                    failures.append(f"invalid visualization: {artifact}")
        if not all(manifest.get("checks", {}).values()):
            failures.append(f"failed checks: C1Q-{stage}")
    if len(manifests) == 7:
        expected = (
            "passed_c1q0_frozen_admission",
            "completed_c1q1_along_track_transport_closure",
            "completed_c1q2_state_medium_constant_audit",
            "completed_c1q3_b0_canonical_isolation",
            "completed_c1q4_single_factor_shadow_replay",
        )
        for stage, status in enumerate(expected):
            if manifests[stage].get("status") != status:
                failures.append(f"unexpected C1Q-{stage} status")
        if manifests[5].get("status") not in {
                "not_entered_c1q5_no_safe_primary_candidate",
                "blocked_c1q5_safe_candidate_requires_explicit_code_repair"}:
            failures.append("unexpected C1Q-5 status")
        final_report = json.loads(
            (OUT / "c1q_6" / "c1q_6_final_attribution.json").read_text(
                encoding="utf-8"))
        if final_report.get("final_attribution") not in ALLOWED_FINAL:
            failures.append("invalid C1Q final attribution")
        if "manual_decision_required" not in final_report.get("next_gate", ""):
            failures.append("C1Q does not stop at the manual gate")
        export = json.loads(
            (OUT / "c1q_0" / "c1q_0_pharlap_export_summary.json").read_text(
                encoding="utf-8"))
        if (export.get("tol") != [1e-12, 0.001, 1]
                or not np.allclose(export.get("angles_deg", []), ANGLES,
                                   rtol=0.0, atol=1e-12)
                or len(export.get("records", [])) != 4
                or any(item.get("ray_label") != 1
                       for item in export.get("records", []))):
            failures.append("strict PHaRLAP configuration or fixed cases drifted")
        shadow = json.loads(
            (OUT / "c1q_4" / "c1q_4_shadow.json").read_text(encoding="utf-8"))
        rows = shadow.get("metrics", [])
        if len(rows) != 24:
            failures.append("C1Q shadow matrix is incomplete")
        if shadow["documented_interpolation"]["production_eligible"]:
            failures.append("undocumented PHaRLAP interpolation entered production")
        required = {
            "trajectory_rms_3d_km", "trajectory_max_3d_km", "x_rms_km",
            "y_rms_km", "z_rms_km", "landing_dx_km",
            "x_error_closure_fraction",
        }
        if any(not np.isfinite([row[key] for key in required]).all() for row in rows):
            failures.append("non-finite C1Q shadow metric")
        safe = {"pharlap_launch", "exact_b0_rhs",
                "group_path_parameterized", "no_projection"}
        primary = []
        for variant in safe:
            high = [row for row in rows
                    if row["branch_id"] == "high_ray" and row["variant"] == variant]
            low = [row for row in rows
                   if row["branch_id"] == "low_ray" and row["variant"] == variant]
            if (len(high) == len(low) == 2
                    and all(row["x_error_closure_fraction"] >= 0.8
                            and row["trajectory_rms_3d_km"] <= 0.5 for row in high)
                    and all(row["trajectory_rms_3d_km"] <= 0.15 for row in low)):
                primary.append(variant)
        if sorted(primary) != sorted(shadow.get("safe_primary_candidates", [])):
            failures.append("shadow primary-candidate rule mismatch")
        effective = [row for row in rows
                     if row["branch_id"] == "high_ray"
                     and row["variant"] == "pharlap_effective_ne"]
        effective_primary = len(effective) == 2 and all(
            row["x_error_closure_fraction"] >= 0.8 for row in effective)
        if effective_primary != shadow.get("pharlap_effective_field_primary"):
            failures.append("effective-field attribution rule mismatch")
        medium = json.loads(
            (OUT / "c1q_2" / "c1q_2_state_medium.json").read_text(encoding="utf-8"))
        canonical = json.loads(
            (OUT / "c1q_3" / "c1q_3_canonical.json").read_text(encoding="utf-8"))
        if primary or effective_primary:
            expected_final = "identified_but_not_safely_repairable"
        elif (medium["attribution"]
              == "pharlap_internal_field_or_output_semantics_candidate"
              and canonical[
                  "effective_pharlap_n2_closes_momentum_better_than_common_provider"]):
            expected_final = "multiple_unresolved_contributors"
        else:
            expected_final = "pharlap_effective_dynamics_not_recoverable"
        if final_report.get("final_attribution") != expected_final:
            failures.append("final attribution does not follow frozen rules")
    formal = json.loads(FINAL_C1.read_text(encoding="utf-8"))
    if formal.get("status") != "failed_c1_chapman_high_low_external_comparison":
        failures.append("formal C1 failure status changed")
    print(json.dumps({"failures": failures}, indent=2, ensure_ascii=True))
    print(f"A4-C1Q verification: {len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
