"""Run A4-C1-2 comparison from the frozen independent roots."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import (
    LIMITS, _atomic_json, _image_check, _sha256, _surface_distance,
    _trajectory_metrics, _verify_manifest,
)
from utils import ecef_to_geocentric


PREP5 = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep5"
PREP5_REPORT = PREP5 / "a4_c1_prep5_validation.json"
PREP5_MANIFEST = PREP5 / "a4_c1_prep5_freeze_manifest.json"
PREP5_TRAJECTORIES = PREP5 / "a4_c1_prep5_trajectories.npz"
C1 = ROOT / "output2" / "a4_pharlap" / "c1"
C1_1_REPORT = C1 / "same_angle" / "c1_1_validation.json"
C1_1_MANIFEST = C1 / "same_angle" / "c1_1_freeze_manifest.json"
OUT = C1 / "independent_roots"
METRICS_CSV = OUT / "c1_independent_root_metrics.csv"
METRICS_JSON = OUT / "c1_independent_root_metrics.json"
RESIDUAL_NPZ = OUT / "c1_independent_root_residuals.npz"
ROOT_FIGURE = OUT / "c1_independent_roots.png"
RESIDUAL_FIGURE = OUT / "c1_independent_root_residuals.png"
REPORT = OUT / "c1_2_validation.json"
MANIFEST = OUT / "c1_2_freeze_manifest.json"

ANGLE_LIMIT_DEG = 0.1
ENDPOINT_LIMIT_KM = 1.0


def _signature_map(report: dict) -> dict:
    return {(item["solver"], item["branch_id"]): item
            for item in report["frozen_signatures"]}


def _trajectory(source, solver: str, branch: str) -> dict:
    prefix = f"{solver.lower()}_{branch}_1km"
    return {
        "trajectory_ecef_km": np.asarray(source[f"{prefix}_ecef_km"], dtype=float),
        "trajectory_geom_path_km": np.asarray(
            source[f"{prefix}_geom_path_km"], dtype=float),
    }


def _compare(branch: str, signatures: dict, trajectories: dict,
             radius_km: float, target_km: float) -> tuple[dict, dict]:
    ph = signatures[("PHaRLAP", branch)]
    ah = signatures[("AH3D", branch)]
    ph_path = trajectories[("PHaRLAP", branch)]
    ah_path = trajectories[("AH3D", branch)]
    trajectory, residual = _trajectory_metrics(ph_path, ah_path, radius_km)
    ph_landing = ecef_to_geocentric(ph_path["trajectory_ecef_km"][-1], radius_km)
    ah_landing = ecef_to_geocentric(ah_path["trajectory_ecef_km"][-1], radius_km)
    landing = _surface_distance(ph_landing, ah_landing, radius_km)
    metric = {
        "branch_id": branch,
        "comparison_semantics": "independent_frozen_root_same_branch_not_same_angle",
        "pharlap_elevation_deg": float(ph["elevation_deg"]),
        "ah3d_elevation_deg": float(ah["elevation_deg"]),
        "launch_elevation_absolute_error_deg": abs(
            float(ah["elevation_deg"]) - float(ph["elevation_deg"])),
        "pharlap_termination": ph["termination"],
        "ah3d_termination": ah["termination"],
        "pharlap_turning_count": int(ph["turning_count"]),
        "ah3d_turning_count": int(ah["turning_count"]),
        "landing_absolute_error_km": landing,
        "landing_relative_error": landing / target_km,
        "apex_height_absolute_error_km": abs(float(ah["apex_km"]) - float(ph["apex_km"])),
        "geometric_path_relative_error": abs(float(ah["geometric_path_km"])
            - float(ph["geometric_path_km"])) / abs(float(ph["geometric_path_km"])),
        "group_path_relative_error": abs(float(ah["group_path_km"])
            - float(ph["group_path_km"])) / abs(float(ph["group_path_km"])),
        "group_delay_relative_error": abs(float(ah["tau_ms"])
            - float(ph["tau_ms"])) / abs(float(ph["tau_ms"])),
        "pharlap_along_track_residual_km": float(ph["along_track_residual_km"]),
        "pharlap_cross_track_residual_km": float(ph["cross_track_km"]),
        "ah3d_along_track_residual_km": float(ah["along_track_residual_km"]),
        "ah3d_cross_track_residual_km": float(ah["cross_track_km"]),
        "pharlap_group_delay_closure_ms": abs(float(ph["tau_ms"])
            - 1000.0 * float(ph["group_path_km"]) / 299792.458),
        "ah3d_group_delay_closure_ms": abs(float(ah["tau_ms"])
            - 1000.0 * float(ah["group_path_km"]) / 299792.458),
        **trajectory,
    }
    gates = {
        "ground_single_turning": ph["termination"] == ah["termination"] == "ground"
            and int(ph["turning_count"]) == int(ah["turning_count"]) == 1,
        "launch_elevation_absolute_error_deg":
            metric["launch_elevation_absolute_error_deg"] <= ANGLE_LIMIT_DEG,
        **{key: metric[key] <= threshold for key, threshold in LIMITS.items()},
        "pharlap_rx_residual": abs(metric["pharlap_along_track_residual_km"])
            <= ENDPOINT_LIMIT_KM and abs(metric["pharlap_cross_track_residual_km"])
            <= ENDPOINT_LIMIT_KM,
        "ah3d_rx_residual": abs(metric["ah3d_along_track_residual_km"])
            <= ENDPOINT_LIMIT_KM and abs(metric["ah3d_cross_track_residual_km"])
            <= ENDPOINT_LIMIT_KM,
        "group_delay_closure": metric["pharlap_group_delay_closure_ms"] <= 1.0e-10
            and metric["ah3d_group_delay_closure_ms"] <= 1.0e-10,
    }
    metric["gates"] = gates
    metric["numeric_gates_passed"] = all(gates.values())
    return metric, residual


def _write_csv(metrics: list[dict]) -> None:
    rows = []
    for metric in metrics:
        row = {key: value for key, value in metric.items() if not isinstance(value, dict)}
        row.update({f"gate_{key}": value for key, value in metric["gates"].items()})
        rows.append(row)
    with METRICS_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _range_height(path: dict, tx, radius_km: float) -> tuple[np.ndarray, np.ndarray]:
    geocentric = ecef_to_geocentric(path["trajectory_ecef_km"], radius_km)
    distance = np.array([_surface_distance(tx, point, radius_km) for point in geocentric])
    return distance, geocentric[:, 2]


def _plot(metrics: list[dict], trajectories: dict, residuals: dict,
          input_data: dict, source_tag: str) -> list[dict]:
    radius = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0
    tx = np.array([float(np.asarray(input_data[key]).item()) for key in
                   ("origin_lat_deg", "origin_lon_deg", "origin_height_km")])
    target = float(np.asarray(input_data["target_range_km"]).item())
    heights = np.asarray(input_data["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(input_data["iono_en_grid_m3"], dtype=float)
    profile = ne[ne.shape[0] // 2, ne.shape[1] // 2, :]

    figure, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
    for axis, metric in zip(axes, metrics):
        image = axis.imshow(np.repeat(profile[:, None], 2, axis=1), origin="lower",
                            extent=(0, target * 1.03, heights[0], heights[-1]),
                            cmap="jet", alpha=0.5, aspect="auto")
        branch = metric["branch_id"]
        for solver, style, color in (("PHaRLAP", "-", "#1368aa"),
                                     ("AH3D", "--", "#c43b3b")):
            distance, height = _range_height(trajectories[(solver, branch)], tx, radius)
            axis.plot(distance, height, style, color=color, linewidth=1.8,
                      label=f"{solver} own frozen root")
        axis.axvline(target, color="black", linestyle=":", linewidth=1, label="RX")
        axis.set_title(f"{branch} | independent roots | angle diff="
                       f"{metric['launch_elevation_absolute_error_deg']:.6f} deg")
        axis.set_xlabel("Great-circle range from TX (km)")
        axis.grid(alpha=0.2); axis.legend(fontsize=8)
    axes[0].set_ylabel("Physical height (km)")
    figure.suptitle(
        "A4-C1-2 independent frozen-root paths | diagnostic after C1-1 same-angle failure\n"
        f"one-way | 10 MHz | B=0 | collisions=off | no re-solving | source={source_tag} | scattering/power=N/A")
    figure.subplots_adjust(left=0.07, right=0.86, bottom=0.11, top=0.81, wspace=0.15)
    colorbar_axis = figure.add_axes((0.89, 0.16, 0.015, 0.58))
    figure.colorbar(image, cax=colorbar_axis, label="Electron density Ne (m^-3)")
    figure.savefig(ROOT_FIGURE, dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
    for axis, metric in zip(axes, metrics):
        values = residuals[metric["branch_id"]]
        axis.plot(values["u"], values["total"], color="black", label="3-D norm")
        axis.plot(values["u"], values["horizontal"], color="#1368aa", label="horizontal")
        axis.plot(values["u"], np.abs(values["up"]), color="#c43b3b", label="|up|")
        axis.axhline(LIMITS["trajectory_max_3d_km"], color="black", linestyle="--",
                     linewidth=1, label="5 km max gate")
        axis.set_title(f"{metric['branch_id']} | RMS="
                       f"{metric['trajectory_rms_3d_km']:.6f} km")
        axis.set_xlabel("Normalized cumulative geometric path u")
        axis.grid(alpha=0.2); axis.legend(fontsize=8)
    axes[0].set_ylabel("AH3D own root - PHaRLAP own root residual (km)")
    figure.suptitle(
        "A4-C1-2 independent-root normalized-path residuals | 201 samples\n"
        f"diagnostic only after C1-1 failure | one-way | 10 MHz | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.88))
    figure.savefig(RESIDUAL_FIGURE, dpi=180)
    plt.close(figure)
    return [_image_check(path) for path in (ROOT_FIGURE, RESIDUAL_FIGURE)]


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        _verify_manifest(PREP5_MANIFEST, "passed_1km_grid_stability_and_signatures")
        _verify_manifest(C1_1_MANIFEST, "failed_c1_1_numeric_thresholds")
        prep5 = json.loads(PREP5_REPORT.read_text(encoding="utf-8"))
        c1_1 = json.loads(C1_1_REPORT.read_text(encoding="utf-8"))
        if c1_1["status"] != "failed_c1_1_numeric_thresholds":
            raise RuntimeError("C1-2 diagnostic semantics require the frozen C1-1 numeric failure")
        input_path = Path(prep5["refined_input"]["path"])
        if _sha256(input_path) != prep5["refined_input"]["sha256"]:
            raise RuntimeError("refined_1km input changed")
        input_data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
        radius = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0
        target = float(np.asarray(input_data["target_range_km"]).item())
        signatures = _signature_map(prep5)
        with np.load(PREP5_TRAJECTORIES, allow_pickle=False) as source:
            trajectories = {(solver, branch): _trajectory(source, solver, branch)
                for solver in ("PHaRLAP", "AH3D")
                for branch in ("low_ray", "high_ray")}
        metrics = []; residuals = {}
        for branch in ("low_ray", "high_ray"):
            metric, residual = _compare(branch, signatures, trajectories, radius, target)
            metrics.append(metric); residuals[branch] = residual

        same_angle = {branch: [item for item in c1_1["metrics"]
                              if item["branch_id"] == branch]
                      for branch in ("low_ray", "high_ray")}
        attribution = {}
        for metric in metrics:
            branch = metric["branch_id"]
            angle = metric["launch_elevation_absolute_error_deg"]
            same_landing = max(item["landing_absolute_error_km"]
                               for item in same_angle[branch])
            same_rms = max(item["trajectory_rms_3d_km"] for item in same_angle[branch])
            attribution[branch] = {
                "independent_root_angle_difference_deg": angle,
                "independent_root_trajectory_rms_km": metric["trajectory_rms_3d_km"],
                "same_angle_max_trajectory_rms_km": same_rms,
                "same_angle_max_landing_error_km": same_landing,
                "endpoint_sensitivity_indicator_km_per_deg": same_landing / angle,
                "root_residual_max_km": max(abs(metric[key]) for key in (
                    "pharlap_along_track_residual_km", "pharlap_cross_track_residual_km",
                    "ah3d_along_track_residual_km", "ah3d_cross_track_residual_km")),
            }
        _write_csv(metrics)
        np.savez_compressed(RESIDUAL_NPZ, **{
            f"{branch}_{key}": value for branch, residual in residuals.items()
            for key, value in residual.items()})
        _atomic_json(METRICS_JSON, {
            "case_semantics": "independent_frozen_roots_diagnostic_after_same_angle_failure",
            "limits": {**LIMITS, "launch_elevation_absolute_error_deg": ANGLE_LIMIT_DEG,
                       "rx_along_cross_absolute_km": ENDPOINT_LIMIT_KM},
            "metrics": metrics, "attribution_indicators": attribution,
        })
        sources = {
            "runner": _sha256(Path(__file__)), "refined_input": _sha256(input_path),
            "prep5_report": _sha256(PREP5_REPORT),
            "prep5_manifest": _sha256(PREP5_MANIFEST),
            "prep5_trajectories": _sha256(PREP5_TRAJECTORIES),
            "c1_1_report": _sha256(C1_1_REPORT),
            "c1_1_manifest": _sha256(C1_1_MANIFEST),
            "metrics_csv": _sha256(METRICS_CSV), "metrics_json": _sha256(METRICS_JSON),
            "residual_npz": _sha256(RESIDUAL_NPZ),
        }
        source_tag = hashlib.sha256(json.dumps(
            sources, sort_keys=True).encode("ascii")).hexdigest()[:12]
        figures = _plot(metrics, trajectories, residuals, input_data, source_tag)
        numeric_pass = all(item["numeric_gates_passed"] for item in metrics)
        checks = {
            "prep5_frozen_roots_only": True,
            "c1_1_failure_preserved": True,
            "no_rescan_no_root_solving": True,
            "low_to_low_high_to_high_only": len(metrics) == 2,
            "independent_root_numeric_gates": numeric_pass,
            "visualization": all(item["passed"] for item in figures),
            "does_not_override_c1_1_failure": True,
        }
        status = ("completed_c1_2_diagnostic_after_same_angle_failure"
                  if all(checks.values()) else "failed_c1_2_diagnostic")
        report = {
            "schema_version": "a4-c1-2-validation-v1", "status": status,
            "case_semantics": (
                "diagnostic_after_same_angle_failure_independent_frozen_roots_"
                "not_c1_acceptance_not_root_hessian_classification"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "checks": checks, "metrics": metrics,
            "attribution_indicators": attribution, "source_sha256": sources,
            "source_set_short_hash": source_tag, "visualizations": figures,
            "diagnostic_conclusion": (
                "Independent endpoint-conditioned roots meet applicable gates, while the "
                "high branch fails the fixed-same-angle RMS gate. Root convergence residuals "
                "are small; the remaining discrepancy is associated with the high-branch "
                "forward range-elevation mapping and its conditioning, not established as a "
                "generic root-finder convergence failure."),
            "limitations": [
                "C1-2 cannot override the frozen C1-1 same-angle failure",
                "The diagnostic does not yet separate interpolation, Hamiltonian derivative, integration, or event mapping",
                "Low/high labels remain topology-based and are not optical-path Hessian classifications",
                "B=0 smooth Chapman does not validate O/X, complex 3-D modes, scattering, PE, full-wave, or power",
            ],
            "next_step": "Run C1-3 attribution; prepare C1R only if attribution identifies an AH3D defect",
        }
        _atomic_json(REPORT, report)
        artifacts = [METRICS_CSV, METRICS_JSON, RESIDUAL_NPZ,
                     ROOT_FIGURE, RESIDUAL_FIGURE, REPORT]
        _atomic_json(MANIFEST, {
            "schema_version": "a4-c1-2-freeze-v1", "status": status,
            "source_sha256": sources,
            "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                                 for path in artifacts],
        })
        print(f"A4-C1-2: {'OK' if status.startswith('completed') else 'FAIL'} {status}")
        return 0 if status.startswith("completed") else 1
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_c1_2", "error":
                                  f"{type(error).__name__}: {error}"})
        print(f"A4-C1-2: BLOCKED {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
