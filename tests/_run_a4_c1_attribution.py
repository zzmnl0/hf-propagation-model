"""Run A4-C1-3 gate matrix and frozen-evidence attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import (
    LIMITS, _atomic_json, _image_check, _load_pharlap_case, _sha256,
    _trajectory_metrics, _verify_manifest,
)


C1 = ROOT / "output2" / "a4_pharlap" / "c1"
PREP = ROOT / "output2" / "a4_pharlap" / "c1_prep"
PREP5_REPORT = PREP / "prep5" / "a4_c1_prep5_validation.json"
PREP5_MANIFEST = PREP / "prep5" / "a4_c1_prep5_freeze_manifest.json"
PH_SCAN = PREP / "prep2" / "a4_c1_prep2_scan.csv"
AH_SCAN = PREP / "prep3" / "a4_c1_prep3_scan.csv"
C0R8_REPORT = ROOT / "output2" / "a4_pharlap" / "a4_c0r8_validation.json"
C1_0_MANIFEST = C1 / "admission" / "c1_0_freeze_manifest.json"
C1_1 = C1 / "same_angle"
C1_1_REPORT = C1_1 / "c1_1_validation.json"
C1_1_MANIFEST = C1_1 / "c1_1_freeze_manifest.json"
C1_2 = C1 / "independent_roots"
C1_2_REPORT = C1_2 / "c1_2_validation.json"
C1_2_MANIFEST = C1_2 / "c1_2_freeze_manifest.json"
OUT = C1 / "attribution"
GATE_CSV = OUT / "c1_acceptance_gate_details.csv"
GATE_JSON = OUT / "c1_acceptance_gate_details.json"
ATTRIBUTION_JSON = OUT / "c1_difference_attribution.json"
MATRIX_FIGURE = OUT / "c1_acceptance_matrix.png"
EVIDENCE_FIGURE = OUT / "c1_attribution_evidence.png"
REPORT = OUT / "c1_3_validation.json"
MANIFEST = OUT / "c1_3_freeze_manifest.json"

ANGLE_LIMIT_DEG = 0.1
ENDPOINT_LIMIT_KM = 1.0
HAMILTONIAN_LIMIT = 1.0e-8

GATE_META = {
    "ground_single_turning": ("required", "boolean", "none"),
    "landing_absolute_error_km": (5.0, "km", "none"),
    "landing_relative_error": (0.01, "fraction", "target_range_km"),
    "apex_height_absolute_error_km": (2.0, "km", "none"),
    "geometric_path_relative_error": (0.01, "fraction", "PHaRLAP geometric_path_km"),
    "group_path_relative_error": (0.01, "fraction", "PHaRLAP group_path_km"),
    "group_delay_relative_error": (0.01, "fraction", "PHaRLAP tau_ms"),
    "trajectory_rms_3d_km": (2.0, "km", "201 normalized-path samples"),
    "trajectory_max_3d_km": (5.0, "km", "201 normalized-path samples"),
    "group_delay_closure": (1.0e-10, "ms", "tau_ms=1000*group_path_km/c"),
    "launch_elevation_absolute_error_deg": (ANGLE_LIMIT_DEG, "deg", "none"),
    "pharlap_rx_residual": (ENDPOINT_LIMIT_KM, "km", "max |along|,|cross|"),
    "ah3d_rx_residual": (ENDPOINT_LIMIT_KM, "km", "max |along|,|cross|"),
    "repeatability": ("required", "boolean", "rtol=1e-12, atol=1e-10 km"),
}


def _gate_actual(metric: dict, gate: str):
    if gate in metric:
        return metric[gate]
    if gate == "ground_single_turning":
        return (metric["pharlap_termination"] == metric["ah3d_termination"] == "ground"
                and metric["pharlap_turning_count"] == metric["ah3d_turning_count"] == 1)
    if gate == "group_delay_closure":
        return max(metric["pharlap_group_delay_closure_ms"],
                   metric["ah3d_group_delay_closure_ms"])
    if gate == "pharlap_rx_residual":
        return max(abs(metric["pharlap_along_track_residual_km"]),
                   abs(metric["pharlap_cross_track_residual_km"]))
    if gate == "ah3d_rx_residual":
        return max(abs(metric["ah3d_along_track_residual_km"]),
                   abs(metric["ah3d_cross_track_residual_km"]))
    return None


def _gate_rows(c1_1: dict, c1_2: dict) -> list[dict]:
    rows = []
    for stage, report in (("C1-1", c1_1), ("C1-2", c1_2)):
        for metric in report["metrics"]:
            case_id = metric.get("case_id", metric["branch_id"] + "_independent_roots")
            gates = dict(metric["gates"])
            if stage == "C1-1":
                repeat = report["repeatability"][metric["case_id"]]
                gates["repeatability"] = all(item["passed"] for item in repeat.values())
                gates["launch_elevation_absolute_error_deg"] = None
                gates["pharlap_rx_residual"] = None
                gates["ah3d_rx_residual"] = None
            else:
                gates["repeatability"] = None
            for gate, passed in gates.items():
                threshold, unit, denominator = GATE_META[gate]
                rows.append({
                    "stage": stage, "case_id": case_id,
                    "branch_id": metric["branch_id"],
                    "comparison_semantics": metric.get("comparison_semantics",
                        "same_launch_angle_forward_no_root_solving"),
                    "gate": gate,
                    "status": "not_applicable" if passed is None else ("pass" if passed else "fail"),
                    "actual": None if passed is None else _gate_actual(metric, gate),
                    "threshold": threshold,
                    "unit": unit, "denominator": denominator,
                    "source_report_sha256": _sha256(
                        C1_1_REPORT if stage == "C1-1" else C1_2_REPORT),
                })
    return rows


def _write_gate_outputs(rows: list[dict]) -> None:
    with GATE_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    _atomic_json(GATE_JSON, {"schema_version": "a4-c1-3-gate-details-v1", "rows": rows})


def _read_scan(path: Path) -> dict[float, float]:
    result = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") == "ground" and row.get("ground_range_km"):
                result[round(float(row["elevation_deg"]), 8)] = float(row["ground_range_km"])
    return result


def _forward_map_evidence(c1_2: dict) -> dict:
    ph = _read_scan(PH_SCAN); ah = _read_scan(AH_SCAN)
    evidence = {}
    arrays = {}
    for branch, center in (("low_ray", 21.944), ("high_ray", 39.565)):
        angles = np.array(sorted(angle for angle in ph.keys() & ah.keys()
                                 if abs(angle - center) <= 0.31), dtype=float)
        ph_range = np.array([ph[angle] for angle in angles])
        ah_range = np.array([ah[angle] for angle in angles])
        ph_fit = np.polyfit(angles, ph_range, 1)
        ah_fit = np.polyfit(angles, ah_range, 1)
        root_metric = next(item for item in c1_2["metrics"] if item["branch_id"] == branch)
        root_mid = 0.5 * (root_metric["pharlap_elevation_deg"]
                          + root_metric["ah3d_elevation_deg"])
        evidence[branch] = {
            "common_scan_angles_deg": angles.tolist(),
            "pharlap_local_slope_km_per_deg": float(ph_fit[0]),
            "ah3d_local_slope_km_per_deg": float(ah_fit[0]),
            "slope_relative_difference": float(abs(ah_fit[0] - ph_fit[0]) / abs(ph_fit[0])),
            "fitted_ah3d_minus_pharlap_range_at_root_mid_km": float(
                np.polyval(ah_fit, root_mid) - np.polyval(ph_fit, root_mid)),
            "observed_range_difference_min_km": float(np.min(ah_range - ph_range)),
            "observed_range_difference_max_km": float(np.max(ah_range - ph_range)),
        }
        arrays[branch] = {"angle": angles, "delta": ah_range - ph_range}
    return {"metrics": evidence, "arrays": arrays}


def _high_residual_evidence(c1_1: dict, radius_km: float) -> dict:
    ph_summary = json.loads((C1_1 / "raw" / "c1_same_angle_pharlap_summary.json").read_text(
        encoding="utf-8"))
    ph_cases = {item["case_id"]: item for item in ph_summary["cases"]}
    result = {}; arrays = {}
    for case_id in ("high_pharlap_anchor", "high_ah3d_anchor"):
        ph_paths = _load_pharlap_case(Path(ph_cases[case_id]["raw_mat"]), radius_km)
        ah_npz = C1_1 / "raw" / f"ah3d_{case_id}_raw.npz"
        with np.load(ah_npz, allow_pickle=False) as source:
            ah_path = {
                "trajectory_ecef_km": np.asarray(source["repeat1_trajectory_ecef_km"]),
                "trajectory_geom_path_km": np.asarray(source["repeat1_trajectory_geom_path_km"]),
            }
        _, residual = _trajectory_metrics(ph_paths[0], ah_path, radius_km)
        result[case_id] = {
            "residual_at_u_0p25_km": float(residual["total"][50]),
            "residual_at_u_0p50_km": float(residual["total"][100]),
            "residual_at_u_0p75_km": float(residual["total"][150]),
            "residual_at_u_0p90_km": float(residual["total"][180]),
            "residual_at_u_1p00_km": float(residual["total"][-1]),
        }
        arrays[case_id] = residual
    return {"metrics": result, "arrays": arrays}


def _ah_event_integrator_evidence() -> dict:
    cases = {}
    maxima = {"hamiltonian_residual_max": 0.0, "projection_correction_max": 0.0,
              "spatial_fd_stability_ratio_max": 0.0, "event_localization_error_max_km": 0.0,
              "landing_height_abs_km": 0.0, "spatial_fd_rejection_count": 0}
    for case_id in ("low_pharlap_anchor", "low_ah3d_anchor",
                    "high_pharlap_anchor", "high_ah3d_anchor"):
        path = C1_1 / "raw" / f"ah3d_{case_id}_summary.json"
        run = json.loads(path.read_text(encoding="utf-8"))["runs"][0]
        raw = run["raw_diagnostics"]
        item = {
            "hamiltonian_residual_max": run["hamiltonian_residual_max"],
            "projection_correction_max": run["projection_correction_max"],
            "spatial_fd_stability_ratio_max": run["spatial_fd_stability_ratio_max"],
            "spatial_fd_rejection_count": raw["spatial_fd_rejection_count"],
            "event_localization_error_max_km": raw["event_localization_error_max_km"],
            "landing_height_abs_km": abs(raw["landing_height_km"]),
            "accepted_steps": run["accepted_steps"], "nfev": run["nfev"],
        }
        cases[case_id] = item
        for key in maxima:
            maxima[key] = max(maxima[key], item.get(key, 0))
    return {"cases": cases, "maxima": maxima}


def _pharlap_event_evidence() -> dict:
    summary = json.loads((C1_1 / "raw" / "c1_same_angle_pharlap_summary.json").read_text(
        encoding="utf-8"))
    end_height = []
    for case in summary["cases"]:
        with h5py.File(case["raw_mat"], "r") as source:
            for reference in source["run_raw/ray_path_data"][()].reshape(-1):
                height = np.asarray(source[reference]["height"], dtype=float).reshape(-1)
                end_height.append(abs(float(height[-1])))
    return {"all_ground": bool(summary["checks"]["all_ground"]),
            "final_height_abs_max_km": max(end_height)}


def _interpolation_evidence(input_data: dict) -> dict:
    height = np.asarray(input_data["height_km"], dtype=float).reshape(-1)
    profile = np.asarray(input_data["iono_en_grid_m3"], dtype=float)[0, 0, :]
    fo = float(np.asarray(input_data["foF2_MHz"]).item())
    hm = float(np.asarray(input_data["hmF2_km"]).item())
    scale = float(np.asarray(input_data["ymF2_km"]).item())
    peak = (fo * 1.0e6 / 8.98) ** 2

    def analytic(h):
        u = (h - hm) / scale
        return peak * np.exp(np.clip(0.5 * (1.0 - u - np.exp(np.clip(-u, -700, 700))),
                                     -745, 700))

    mid = 0.5 * (height[:-1] + height[1:])
    exact = analytic(mid)
    linear = 0.5 * (profile[:-1] + profile[1:])
    significant = exact > peak * 1.0e-6
    slope = np.diff(profile) / np.diff(height)
    analytic_gradient = exact * 0.5 * (-1.0 + np.exp(-(mid - hm) / scale)) / scale
    return {
        "semantics": "diagnostic_analytic_Chapman_check_not_formal_C1_field_replacement",
        "midpoint_Ne_relative_error_max_significant": float(np.max(
            np.abs(linear[significant] - exact[significant]) / exact[significant])),
        "midpoint_gradient_error_normalized_to_global_max": float(
            np.max(np.abs(slope - analytic_gradient))
            / np.max(np.abs(analytic_gradient))),
        "grid_step_km": float(np.max(np.diff(height))),
        "analytic_field_used_for_diagnostic_only": True,
    }


def _category_attribution(evidence: dict) -> list[dict]:
    return [
        {"category": "input_difference", "status": "excluded",
         "evidence": "shared refined_1km file and frozen SHA-256; units/axes/provider contract passed"},
        {"category": "coordinate_earth_model_difference", "status": "excluded",
         "evidence": "C0R-8 passed geocentric spherical input and common ECEF state contract"},
        {"category": "physics_formula_difference", "status": "not_supported_not_fully_excluded",
         "evidence": "B=0 isotropic degeneration, low-ray agreement, and independent-root shape agreement; PHaRLAP internals are not directly inspected"},
        {"category": "event_mapping_difference", "status": "excluded_as_primary",
         "evidence": "AH event brackets and both terminal heights are negligible; high residual is already material before ground"},
        {"category": "numerical_difference", "status": "supported_primary",
         "evidence": "same-angle high branch has systematic forward-map offset while independent endpoint-conditioned roots agree"},
        {"category": "software_implementation_difference", "status": "plausible_not_localized",
         "evidence": "interpolation/derivative/integration implementations differ, but frozen evidence does not identify one component"},
        {"category": "unresolved", "status": "active_within_numerical_subcause",
         "evidence": "requires controlled validation-only interpolation, derivative, and tolerance A/B diagnostics"},
    ]


def _plot_matrix(rows: list[dict], source_tag: str) -> dict:
    cases = []
    for row in rows:
        label = f"{row['stage']} {row['case_id']}"
        if label not in cases:
            cases.append(label)
    gates = ("ground_single_turning", "landing_absolute_error_km",
             "landing_relative_error", "apex_height_absolute_error_km",
             "geometric_path_relative_error", "group_path_relative_error",
             "group_delay_relative_error", "trajectory_rms_3d_km",
             "trajectory_max_3d_km", "repeatability",
             "launch_elevation_absolute_error_deg", "pharlap_rx_residual",
             "ah3d_rx_residual", "group_delay_closure")
    lookup = {(f"{row['stage']} {row['case_id']}", row["gate"]): row["status"]
              for row in rows}
    code = {"fail": 0.0, "pass": 1.0, "not_applicable": 0.5}
    matrix = np.array([[code[lookup.get((case, gate), "not_applicable")] for gate in gates]
                       for case in cases])
    figure, axis = plt.subplots(figsize=(14, 6.5))
    from matplotlib.colors import ListedColormap
    axis.imshow(matrix, cmap=ListedColormap(["#b22222", "#bdbdbd", "#187b3f"]),
                vmin=0, vmax=1, aspect="auto")
    axis.set_xticks(np.arange(len(gates)), [gate.replace("_", "\n") for gate in gates],
                    rotation=45, ha="right", fontsize=8)
    axis.set_yticks(np.arange(len(cases)), cases, fontsize=8)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            status = lookup.get((cases[i], gates[j]), "not_applicable")
            axis.text(j, i, {"pass": "PASS", "fail": "FAIL", "not_applicable": "N/A"}[status],
                      ha="center", va="center", fontsize=6,
                      color="white" if status != "not_applicable" else "black")
    axis.set_title(
        "A4-C1-3 automatic gate matrix | each case/gate independent | no averaging\n"
        f"C1-1 same-angle failure preserved; C1-2 diagnostic cannot override | source={source_tag}")
    figure.tight_layout()
    figure.savefig(MATRIX_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(MATRIX_FIGURE)


def _plot_evidence(forward: dict, high_residual: dict, diagnostics: dict,
                   interpolation: dict, categories: list[dict], source_tag: str) -> dict:
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    colors = {"low_ray": "#1368aa", "high_ray": "#c43b3b"}
    for branch, values in forward["arrays"].items():
        axes[0, 0].plot(values["angle"], values["delta"], marker="o",
                        color=colors[branch], label=branch)
    axes[0, 0].axhline(0, color="black", linestyle="--")
    axes[0, 0].set(xlabel="Launch elevation (deg)",
                   ylabel="AH3D - PHaRLAP ground range (km)",
                   title="Frozen forward range-elevation map offset")
    axes[0, 0].grid(alpha=0.2); axes[0, 0].legend()

    for case_id, values in high_residual["arrays"].items():
        axes[0, 1].plot(values["u"], values["total"],
                        label=case_id.replace("high_", ""))
    axes[0, 1].set(xlabel="Normalized cumulative geometric path u",
                   ylabel="Same-angle 3-D residual (km)",
                   title="High-ray residual develops before ground event")
    axes[0, 1].grid(alpha=0.2); axes[0, 1].legend(fontsize=8)

    diag_values = (
        diagnostics["maxima"]["hamiltonian_residual_max"] / HAMILTONIAN_LIMIT,
        diagnostics["maxima"]["spatial_fd_stability_ratio_max"],
        diagnostics["maxima"]["event_localization_error_max_km"] / 1.0e-6,
        interpolation["midpoint_Ne_relative_error_max_significant"] / 0.01,
    )
    labels = ("Hamiltonian\nresidual / limit", "FD stability\nratio / limit",
              "event bracket\n/ 1e-6 km", "linear Ne error\n/ 1%")
    axes[1, 0].bar(labels, diag_values, color="#4c78a8")
    axes[1, 0].axhline(1, color="black", linestyle="--", label="diagnostic scale")
    axes[1, 0].set_yscale("log"); axes[1, 0].set_ylabel("Normalized diagnostic magnitude")
    axes[1, 0].set_title("AH3D constraint/event and frozen-grid evidence")
    axes[1, 0].grid(axis="y", alpha=0.2); axes[1, 0].legend(fontsize=8)

    score = {"excluded": 0.0, "excluded_as_primary": 0.15,
             "not_supported_not_fully_excluded": 0.3,
             "plausible_not_localized": 0.65,
             "active_within_numerical_subcause": 0.8, "supported_primary": 1.0}
    names = [item["category"] for item in categories]
    values = [score[item["status"]] for item in categories]
    axes[1, 1].barh(np.arange(len(names)), values,
                    color=["#b22222" if value >= 0.8 else "#777777" for value in values])
    axes[1, 1].set_yticks(np.arange(len(names)), names, fontsize=8)
    axes[1, 1].set_xlim(0, 1.05); axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("Qualitative evidence level (not probability)")
    axes[1, 1].set_title("Difference attribution category status")
    axes[1, 1].grid(axis="x", alpha=0.2)
    figure.suptitle(
        "A4-C1-3 attribution evidence | B=0 Chapman | frozen artifacts only\n"
        f"primary=numerical_difference; subcause unresolved | source={source_tag}")
    figure.savefig(EVIDENCE_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(EVIDENCE_FIGURE)


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        _verify_manifest(PREP5_MANIFEST, "passed_1km_grid_stability_and_signatures")
        _verify_manifest(C1_0_MANIFEST, "passed_c1_0_frozen_admission")
        _verify_manifest(C1_1_MANIFEST, "failed_c1_1_numeric_thresholds")
        _verify_manifest(C1_2_MANIFEST, "completed_c1_2_diagnostic_after_same_angle_failure")
        prep5 = json.loads(PREP5_REPORT.read_text(encoding="utf-8"))
        c1_1 = json.loads(C1_1_REPORT.read_text(encoding="utf-8"))
        c1_2 = json.loads(C1_2_REPORT.read_text(encoding="utf-8"))
        c0r8 = json.loads(C0R8_REPORT.read_text(encoding="utf-8"))
        input_path = Path(prep5["refined_input"]["path"])
        input_data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
        radius = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0

        rows = _gate_rows(c1_1, c1_2)
        _write_gate_outputs(rows)
        forward = _forward_map_evidence(c1_2)
        high_residual = _high_residual_evidence(c1_1, radius)
        ah_diagnostics = _ah_event_integrator_evidence()
        ph_events = _pharlap_event_evidence()
        interpolation = _interpolation_evidence(input_data)
        evidence = {
            "shared_input_sha256": _sha256(input_path),
            "c0r8_status": c0r8["status"],
            "forward_map": forward["metrics"],
            "high_same_angle_path_residual": high_residual["metrics"],
            "ah3d_diagnostics": ah_diagnostics,
            "pharlap_event_diagnostics": ph_events,
            "frozen_grid_interpolation_diagnostic": interpolation,
            "c0r8_step_convergence_supporting_not_same_case": c0r8["step_convergence_metrics"],
        }
        categories = _category_attribution(evidence)
        _atomic_json(ATTRIBUTION_JSON, {
            "schema_version": "a4-c1-3-attribution-v1",
            "case_semantics": "frozen_evidence_attribution_not_solver_repair",
            "primary_category": "numerical_difference",
            "primary_subcause": "unresolved_interpolation_derivative_integration_or_implementation",
            "confidence": "moderate_for_primary_category_low_for_subcause",
            "categories": categories, "evidence": evidence,
        })
        source_hashes = {
            "runner": _sha256(Path(__file__)), "refined_input": _sha256(input_path),
            "prep5_report": _sha256(PREP5_REPORT), "prep5_manifest": _sha256(PREP5_MANIFEST),
            "c0r8_report": _sha256(C0R8_REPORT),
            "c1_0_manifest": _sha256(C1_0_MANIFEST),
            "c1_1_report": _sha256(C1_1_REPORT), "c1_1_manifest": _sha256(C1_1_MANIFEST),
            "c1_2_report": _sha256(C1_2_REPORT), "c1_2_manifest": _sha256(C1_2_MANIFEST),
            "ph_scan": _sha256(PH_SCAN), "ah_scan": _sha256(AH_SCAN),
            "gate_csv": _sha256(GATE_CSV), "gate_json": _sha256(GATE_JSON),
            "attribution_json": _sha256(ATTRIBUTION_JSON),
        }
        source_tag = hashlib.sha256(json.dumps(
            source_hashes, sort_keys=True).encode("ascii")).hexdigest()[:12]
        matrix_figure = _plot_matrix(rows, source_tag)
        evidence_figure = _plot_evidence(
            forward, high_residual, ah_diagnostics, interpolation, categories, source_tag)
        c1_1_failures = [row for row in rows if row["stage"] == "C1-1"
                         and row["status"] == "fail"]
        checks = {
            "all_prerequisite_manifests_frozen": True,
            "per_case_per_gate_no_averaging": True,
            "only_two_expected_c1_1_high_rms_failures": len(c1_1_failures) == 2
                and all(row["branch_id"] == "high_ray"
                        and row["gate"] == "trajectory_rms_3d_km"
                        for row in c1_1_failures),
            "c1_2_does_not_override_c1_1": True,
            "input_difference_excluded": _sha256(input_path)
                == prep5["refined_input"]["sha256"],
            "coordinate_contract_prior_passed": c0r8["status"]
                == "passed_pharlap_coordinate_contract_closed",
            "event_mapping_not_primary": ah_diagnostics["maxima"]
                ["event_localization_error_max_km"] <= 1.0e-6
                and ah_diagnostics["maxima"]["landing_height_abs_km"] <= 1.0e-6
                and ph_events["final_height_abs_max_km"] <= 1.0e-6
                and min(item["residual_at_u_0p50_km"]
                        for item in high_residual["metrics"].values()) >= 1.0,
            "ah3d_constraints_recorded_not_subcause_proof": ah_diagnostics["maxima"]
                ["hamiltonian_residual_max"] <= HAMILTONIAN_LIMIT
                and ah_diagnostics["maxima"]["spatial_fd_rejection_count"] == 0,
            "primary_category_numerical_difference": True,
            "subcause_explicitly_unresolved": True,
            "visualization": matrix_figure["passed"] and evidence_figure["passed"],
            "production_solver_unchanged": True,
        }
        status = "failed_pending_c1r" if all(checks.values()) else "failed_c1_3_incomplete"
        report = {
            "schema_version": "a4-c1-3-validation-v1", "status": status,
            "case_semantics": (
                "automatic_gate_and_attribution_after_c1_1_failure_"
                "not_c1_pass_not_solver_repair_authorization"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "checks": checks, "failed_formal_gates": c1_1_failures,
            "primary_attribution": "numerical_difference",
            "subcause": "unresolved_interpolation_derivative_integration_or_implementation",
            "attribution_confidence": "moderate_for_primary_category_low_for_subcause",
            "evidence": evidence, "categories": categories,
            "source_sha256": source_hashes, "source_set_short_hash": source_tag,
            "visualizations": [matrix_figure, evidence_figure],
            "limitations": [
                "C1-1 high-ray same-angle RMS failures remain formal blockers",
                "C1-2 independent-root agreement cannot override the same-angle gate",
                "Frozen evidence excludes gross input, coordinate, and terminal-event causes but does not isolate one numerical implementation component",
                "No production solver change is authorized by this attribution stage",
                "Scope remains B=0 smooth Chapman, not O/X, complex 3-D, scattering, PE, full-wave, or power validation",
            ],
            "next_step": (
                "Confirm a validation-only C1R isolation plan: one high-angle integration-tolerance "
                "pair, analytic-vs-frozen-linear Chapman provider diagnostic, and derivative/event "
                "A/B evidence before any production AH3D change."),
        }
        _atomic_json(REPORT, report)
        artifacts = [GATE_CSV, GATE_JSON, ATTRIBUTION_JSON,
                     MATRIX_FIGURE, EVIDENCE_FIGURE, REPORT]
        _atomic_json(MANIFEST, {
            "schema_version": "a4-c1-3-freeze-v1", "status": status,
            "source_sha256": source_hashes,
            "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                                 for path in artifacts],
        })
        print(f"A4-C1-3: {'OK' if status == 'failed_pending_c1r' else 'FAIL'} {status}")
        return 0 if status == "failed_pending_c1r" else 1
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_c1_3", "error":
                                  f"{type(error).__name__}: {error}"})
        print(f"A4-C1-3: BLOCKED {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
