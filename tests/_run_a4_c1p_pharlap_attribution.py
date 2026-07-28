"""Run A4-C1P PHaRLAP-side numerical-attribution diagnostics."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
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

import config as cfg
from models.dispersion_ah import C_KM_S
from models.ray_tracer_ah_3d_ecef import ah3d_ecef_rhs
from tests._run_a4_c1_prep3 import _load_context, _model
from tests._run_a4_c1_validation import (
    CASES,
    _image_check,
    _refined_input,
    _repeatability,
    _sha256,
    _surface_distance,
    _trajectory_metrics,
    _verify_manifest,
)
from utils import ecef_to_geocentric, geocentric_to_ecef


C1 = ROOT / "output2" / "a4_pharlap" / "c1"
OUT = C1 / "c1p"
FINAL_C1_STATUS = "failed_c1_chapman_high_low_external_comparison"
EXPECTED_INPUT_HASH = "0a978478b52f5363c93cbc77b0ee1868a4f50e488c82c829eb7794f0cf1235a7"

C1_MANIFESTS = (
    (C1 / "admission" / "c1_0_freeze_manifest.json", "passed_c1_0_frozen_admission"),
    (C1 / "same_angle" / "c1_1_freeze_manifest.json", "failed_c1_1_numeric_thresholds"),
    (C1 / "independent_roots" / "c1_2_freeze_manifest.json",
     "completed_c1_2_diagnostic_after_same_angle_failure"),
    (C1 / "attribution" / "c1_3_freeze_manifest.json", "failed_pending_c1r"),
    (C1 / "final" / "a4_c1_freeze_manifest.json", FINAL_C1_STATUS),
)
C1R_ALLOWED = {
    0: {"passed_c1r0_failure_signature_frozen"},
    1: {"passed_c1r1_integration_not_primary"},
    2: {"passed_c1r2_frozen_field_not_primary",
        "completed_c1r2_field_effect_unresolved"},
    3: {"passed_c1r3_spatial_derivative_not_primary",
        "completed_c1r3_spatial_derivative_candidate",
        "completed_c1r3_spatial_derivative_unresolved"},
    4: {"completed_c1r4_remaining_implementation_attribution"},
    5: {"completed_c1r5_ah3d_defect_candidate_requires_repair_plan",
        "completed_c1r5_software_difference_or_unresolved"},
}

MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SOURCE = ROOT / "tests" / "matlab" / "a4_c1p_pharlap_scan.m"
MATLAB_SCRIPT = MATLAB_DIR / MATLAB_SOURCE.name
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
AH_SOLVER = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
AH_PROVIDER = ROOT / "models" / "ne_provider_geocentric.py"
CONFIG = ROOT / "config.py"

PROFILE_ROWS = (
    ("baseline", 1e-8, 0.01, 10.0),
    ("loose_tol", 1e-6, 0.01, 10.0),
    ("tight_tol", 1e-10, 0.01, 10.0),
    ("ultratight_tol", 1e-12, 0.01, 10.0),
    ("max_step_5", 1e-10, 0.01, 5.0),
    ("max_step_2", 1e-10, 0.01, 2.0),
    ("max_step_1", 1e-10, 0.01, 1.0),
    ("min_step_0p1", 1e-12, 0.1, 1.0),
    ("min_step_0p01", 1e-12, 0.01, 1.0),
    ("reference", 1e-12, 0.001, 1.0),
)
PROFILE_IDS = tuple(item[0] for item in PROFILE_ROWS)
CASE_IDS = tuple(item[0] for item in CASES)
CASE_META = {item[0]: item for item in CASES}
REFERENCE_PROFILE = "reference"
REFERENCE_NEIGHBOR = "min_step_0p01"


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _source_tag(sources: dict) -> str:
    return hashlib.sha256(json.dumps(
        sources, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()[:12]


def _base_sources(input_path: Path) -> dict:
    return {
        "runner": _sha256(Path(__file__)),
        "matlab_source": _sha256(MATLAB_SOURCE),
        "refined_input": _sha256(input_path),
        "pharlap_mex": _sha256(MEX),
        "ah_solver": _sha256(AH_SOLVER),
        "ah_provider": _sha256(AH_PROVIDER),
        "config": _sha256(CONFIG),
    }


def _verify_any_manifest(path: Path, allowed: set[str]) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") not in allowed:
        raise RuntimeError(f"unexpected manifest status: {path}")
    for item in manifest.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"frozen artifact drift: {artifact}")
    return manifest


def _manifest(stage: int, status: str, sources: dict, artifacts: list[Path],
              checks: dict) -> dict:
    stage_dir = OUT / f"c1p_{stage}"
    payload = {
        "schema_version": f"a4-c1p-{stage}-freeze-v1",
        "status": status,
        "created_utc": _utc(),
        "semantics": (
            "PHaRLAP-side validation-only numerical attribution; one-way; "
            "10 MHz; B=0; collision=0; fixed angles; no root solving; "
            "formal C1 remains failed"
        ),
        "checks": checks,
        "source_sha256": sources,
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in artifacts
        ],
    }
    _write_json(stage_dir / f"c1p_{stage}_freeze_manifest.json", payload)
    return payload


def _require_stage(stage: int, allowed: set[str]) -> dict:
    path = OUT / f"c1p_{stage}" / f"c1p_{stage}_freeze_manifest.json"
    return _verify_any_manifest(path, allowed)


def _plot_status(path: Path, title: str, checks: dict, source_tag: str) -> dict:
    labels = list(checks)
    values = np.array([bool(checks[key]) for key in labels], dtype=int)
    figure, axis = plt.subplots(figsize=(11.5, max(5.0, 0.42 * len(labels))))
    axis.barh(np.arange(len(labels)), values,
              color=np.where(values, "#187b3f", "#b22222"))
    axis.set_yticks(np.arange(len(labels)),
                    [label.replace("_", " ") for label in labels])
    axis.set_xlim(0.0, 1.05)
    axis.set_xticks((0, 1), ("FAIL", "PASS"))
    axis.grid(axis="x", alpha=0.2)
    axis.set_title(f"{title}\nsource={source_tag}")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    result = _image_check(path)
    if not result["passed"]:
        raise RuntimeError(f"invalid visualization: {path}")
    return result


def run_c1p0() -> str:
    stage_dir = OUT / "c1p_0"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, refined = _refined_input()
    sources = _base_sources(input_path)
    manifests = [
        _verify_manifest(path, status) for path, status in C1_MANIFESTS
    ]
    c1r = {}
    for stage, allowed in C1R_ALLOWED.items():
        path = C1 / "c1r" / f"c1r_{stage}" / f"c1r_{stage}_freeze_manifest.json"
        c1r[stage] = _verify_any_manifest(path, allowed)

    c1_0_sources = manifests[0]["source_sha256"]
    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    collision = np.asarray(data["collision_freq_Hz"], dtype=float)
    frozen_angles = np.array([item[3] for item in CASES], dtype=float)
    expected_angles = np.array([
        21.94290326504532, 21.94562872092976,
        39.54254824570964, 39.58701138438118,
    ])
    checks = {
        "c1_manifests_frozen": len(manifests) == len(C1_MANIFESTS),
        "c1r_manifests_frozen": len(c1r) == 6,
        "formal_c1_failure_preserved": manifests[-1]["status"] == FINAL_C1_STATUS,
        "refined_input_hash_exact": sources["refined_input"]
            == refined["sha256"] == EXPECTED_INPUT_HASH,
        "pharlap_mex_hash_frozen": sources["pharlap_mex"]
            == c1_0_sources["pharlap_mex"],
        "ah3d_sources_frozen": all(
            sources[key] == c1_0_sources[key]
            for key in ("ah_solver", "ah_provider", "config")),
        "physics_10mhz_b0_collision0_onehop": (
            float(np.asarray(data["frequency_MHz"]).item()) == 10.0
            and int(np.asarray(data["OX_mode"]).item()) == 0
            and int(np.asarray(data["nhops"]).item()) == 1
            and np.all(collision == 0.0)),
        "spherical_earth_6371km": (
            float(np.asarray(data["earth_radius_m"]).item()) == 6371000.0),
        "four_frozen_angles_exact": np.array_equal(frozen_angles, expected_angles),
        "matlab_and_pharlap_available": MATLAB.is_file() and MEX.is_file(),
    }
    source_tag = _source_tag(sources)
    figure_path = stage_dir / "c1p_0_frozen_admission.png"
    figure = _plot_status(
        figure_path,
        "A4-C1P-0 frozen admission | one-way | 10 MHz | B=0 | collision=0",
        checks, source_tag)
    checks["visualization"] = figure["passed"]
    status = ("passed_c1p0_frozen_admission"
              if all(checks.values()) else "blocked_c1p0_source_or_input_drift")
    report_path = stage_dir / "c1p_0_admission.json"
    _write_json(report_path, {
        "schema_version": "a4-c1p-0-admission-v1",
        "status": status,
        "created_utc": _utc(),
        "checks": checks,
        "frozen_angles_deg": frozen_angles,
        "profile_matrix": PROFILE_ROWS,
        "source_sha256": sources,
        "source_set_short_hash": source_tag,
        "formal_c1_status": manifests[-1]["status"],
        "c1r_statuses": {str(stage): value["status"] for stage, value in c1r.items()},
        "visualization": figure,
    })
    _manifest(0, status, sources, [report_path, figure_path], checks)
    if not status.startswith("passed"):
        raise RuntimeError(status)
    return status


def _deploy_matlab() -> None:
    MATLAB_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MATLAB_SOURCE, MATLAB_SCRIPT)
    if _sha256(MATLAB_SOURCE) != _sha256(MATLAB_SCRIPT):
        raise RuntimeError("deployed MATLAB script differs from source")


def _run_matlab(input_path: Path, raw_dir: Path, summary_path: Path,
                log_path: Path) -> subprocess.CompletedProcess:
    quote = lambda value: str(value).replace("'", "''")
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); a4_c1p_pharlap_scan("
        f"'{quote(input_path)}','{quote(raw_dir)}','{quote(summary_path)}')")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run(
        [str(MATLAB), "-batch", expression], env=env,
        text=True, capture_output=True, timeout=3600, errors="replace")
    log_path.write_text((result.stdout or "") + "\n" + (result.stderr or ""),
                        encoding="utf-8")
    return result


def _h5_vector(group, name: str) -> np.ndarray:
    return np.asarray(group[name], dtype=float).reshape(-1)


def _load_pharlap_raw(path: Path, radius_km: float) -> tuple[list[dict], list[dict]]:
    paths = []
    states = []
    with h5py.File(path, "r") as source:
        path_node = source["run_raw/ray_path_data"]
        state_node = source["run_raw/ray_state_vec"]
        path_groups = ([path_node] if isinstance(path_node, h5py.Group) else
                       [source[reference] for reference in path_node[()].reshape(-1)])
        state_groups = ([state_node] if isinstance(state_node, h5py.Group) else
                        [source[reference] for reference in state_node[()].reshape(-1)])
        if len(path_groups) != len(state_groups):
            raise RuntimeError(f"path/state repeat mismatch: {path}")
        for path_group, state_group in zip(path_groups, state_groups):
            lat = _h5_vector(path_group, "lat")
            lon = _h5_vector(path_group, "lon")
            height = _h5_vector(path_group, "height")
            geom = _h5_vector(path_group, "geometric_distance")
            geom -= geom[0]
            paths.append({
                "trajectory_ecef_km": geocentric_to_ecef(
                    lat, lon, height, radius_km),
                "trajectory_geocentric_lat_lon_height": np.column_stack(
                    (lat, lon, height)),
                "trajectory_physical_height_km": height,
                "trajectory_geom_path_km": geom,
            })
            states.append({
                "position_ecef_km": np.column_stack([
                    _h5_vector(state_group, key)
                    for key in ("pos_x", "pos_y", "pos_z")]) / 1000.0,
                "phase_index_vector": np.column_stack([
                    _h5_vector(state_group, key)
                    for key in ("dir_x", "dir_y", "dir_z")]),
                "group_path_km": _h5_vector(state_group, "group_path") / 1000.0,
                "geometric_path_km": (
                    _h5_vector(state_group, "geometrical_path") / 1000.0),
                "phase_path_km": _h5_vector(state_group, "phase_path") / 1000.0,
                "indep_var_raw": _h5_vector(state_group, "indep_var"),
                "ode_step_size_raw": _h5_vector(state_group, "ODE_step_size"),
            })
    return paths, states


def _record_map(summary: dict) -> dict:
    records = summary["records"]
    if isinstance(records, dict):
        records = [records]
    return {(item["profile_id"], item["case_id"]): item for item in records}


def _ah_path(case_id: str) -> dict:
    path = C1 / "same_angle" / "raw" / f"ah3d_{case_id}_raw.npz"
    with np.load(path, allow_pickle=False) as source:
        return {key: np.asarray(source[f"repeat1_{key}"], dtype=float)
                for key in (
                    "trajectory_ecef_km",
                    "trajectory_geocentric_lat_lon_height",
                    "trajectory_physical_height_km",
                    "trajectory_geom_path_km",
                )}


def _c1_metrics() -> dict:
    report = json.loads(
        (C1 / "same_angle" / "c1_1_validation.json").read_text(encoding="utf-8"))
    if report["status"] != "failed_c1_1_numeric_thresholds":
        raise RuntimeError("formal C1 status changed")
    return {item["case_id"]: item for item in report["metrics"]}


def _load_scan(stage_dir: Path, radius: float) -> tuple[dict, dict, dict]:
    summary = json.loads(
        (stage_dir / "c1p_1_pharlap_scan_summary.json").read_text(encoding="utf-8"))
    records = _record_map(summary)
    paths = {}
    states = {}
    for profile_id in PROFILE_IDS:
        for case_id in CASE_IDS:
            raw = Path(records[(profile_id, case_id)]["raw_mat"])
            paths[(profile_id, case_id)], states[(profile_id, case_id)] = (
                _load_pharlap_raw(raw, radius))
    return records, paths, states


def _step_stats(state: dict) -> dict:
    step = np.abs(np.asarray(state["ode_step_size_raw"], dtype=float)) / 1000.0
    finite = step[np.isfinite(step)]
    return {
        "step_min_km": float(np.min(finite)),
        "step_median_km": float(np.median(finite)),
        "step_p95_km": float(np.percentile(finite, 95)),
        "step_max_km": float(np.max(finite)),
        "indep_var_end_raw": float(state["indep_var_raw"][-1]),
    }


def _scan_repeatability(runs: list[dict], paths: list[dict]) -> dict:
    scalar_keys = (
        "ground_range_km", "along_track_residual_km", "apex_km",
        "geometric_path_km", "group_path_km", "tau_ms",
        "landing_lat_deg", "landing_lon_deg",
    )
    checks = {
        "termination": runs[0]["status"] == runs[1]["status"],
        "turning_count": int(runs[0]["turning_count"]) == int(runs[1]["turning_count"]),
        "point_count": int(runs[0]["path_point_count"]) == int(runs[1]["path_point_count"]),
        "critical_scalars": all(np.isclose(
            float(runs[0][key]), float(runs[1][key]),
            rtol=1.0e-12, atol=1.0e-10) for key in scalar_keys),
        "trajectory": (
            np.asarray(paths[0]["trajectory_ecef_km"]).shape
            == np.asarray(paths[1]["trajectory_ecef_km"]).shape
            and np.allclose(
                paths[0]["trajectory_ecef_km"], paths[1]["trajectory_ecef_km"],
                rtol=1.0e-12, atol=1.0e-10)),
    }
    return {"passed": all(checks.values()), "checks": checks}


def _plot_convergence(stage_dir: Path, rows: list[dict], source_tag: str) -> list[Path]:
    by_key = {(row["profile_id"], row["case_id"]): row for row in rows}
    x = np.arange(len(PROFILE_IDS))
    short = [name.replace("_step_", "\n").replace("_tol", "\ntol")
             for name in PROFILE_IDS]

    first = stage_dir / "c1p_1_trajectory_convergence.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, case_id in zip(axes.flat, CASE_IDS):
        values = [by_key[(profile, case_id)]["trajectory_rms_to_reference_km"]
                  for profile in PROFILE_IDS]
        axis.plot(x, values, marker="o", color="#1368aa")
        axis.axhline(0.1, color="#b22222", linestyle="--",
                     label="0.1 km not-primary limit")
        axis.set_yscale("symlog", linthresh=1e-6)
        axis.set_title(f"{case_id} | {CASE_META[case_id][1]}")
        axis.set_ylabel("PHaRLAP RMS to strict reference (km)")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    for axis in axes[1]:
        axis.set_xticks(x, short, rotation=35, ha="right")
    figure.suptitle(
        "A4-C1P-1 PHaRLAP integration convergence | fixed angles | no root solving\n"
        f"one-way | 10 MHz | B=0 | collision=0 | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(first, dpi=180)
    plt.close(figure)

    second = stage_dir / "c1p_1_path_quantity_convergence.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, case_id in zip(axes.flat, CASE_IDS):
        landing = [by_key[(profile, case_id)]["landing_delta_to_reference_km"]
                   for profile in PROFILE_IDS]
        apex = [by_key[(profile, case_id)]["apex_delta_to_reference_km"]
                for profile in PROFILE_IDS]
        group = [by_key[(profile, case_id)]["group_path_delta_to_reference_km"]
                 for profile in PROFILE_IDS]
        axis.plot(x, landing, marker="o", label="landing")
        axis.plot(x, apex, marker="s", label="apex")
        axis.plot(x, group, marker="^", label="group path")
        axis.set_yscale("symlog", linthresh=1e-6)
        axis.set_title(f"{case_id} | {CASE_META[case_id][1]}")
        axis.set_ylabel("Absolute difference to strict reference (km)")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    for axis in axes[1]:
        axis.set_xticks(x, short, rotation=35, ha="right")
    figure.suptitle(
        "A4-C1P-1 landing/apex/path convergence | PHaRLAP only\n"
        f"one-way | 10 MHz | B=0 | collision=0 | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(second, dpi=180)
    plt.close(figure)

    third = stage_dir / "c1p_1_step_diagnostics.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, case_id in zip(axes.flat, CASE_IDS):
        median = [by_key[(profile, case_id)]["step_median_km"]
                  for profile in PROFILE_IDS]
        p95 = [by_key[(profile, case_id)]["step_p95_km"]
               for profile in PROFILE_IDS]
        axis.plot(x, median, marker="o", label="median")
        axis.plot(x, p95, marker="s", label="p95")
        axis.set_title(f"{case_id} | {CASE_META[case_id][1]}")
        axis.set_ylabel("PHaRLAP ODE step (km)")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    for axis in axes[1]:
        axis.set_xticks(x, short, rotation=35, ha="right")
    figure.suptitle(
        "A4-C1P-1 PHaRLAP ODE step diagnostics | ray_state_vec\n"
        f"one-way | 10 MHz | B=0 | collision=0 | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(third, dpi=180)
    plt.close(figure)
    return [first, second, third]


def run_c1p1() -> str:
    _require_stage(0, {"passed_c1p0_frozen_admission"})
    stage_dir = OUT / "c1p_1"
    raw_dir = stage_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    _deploy_matlab()
    sources = {**_base_sources(input_path),
               "matlab_deployed": _sha256(MATLAB_SCRIPT),
               "c1p0_manifest": _sha256(
                   OUT / "c1p_0" / "c1p_0_freeze_manifest.json")}
    if sources["matlab_source"] != sources["matlab_deployed"]:
        raise RuntimeError("MATLAB deployment hash mismatch")
    summary_path = stage_dir / "c1p_1_pharlap_scan_summary.json"
    log_path = stage_dir / "c1p_1_pharlap_scan.log"
    reusable = summary_path.is_file() and len(list(raw_dir.glob("*.mat"))) == 40
    summary = (json.loads(summary_path.read_text(encoding="utf-8"))
               if reusable else {})
    if summary.get("status") != "completed_c1p1_pharlap_fixed_angle_scan":
        process = _run_matlab(input_path, raw_dir, summary_path, log_path)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if (process.returncode != 0
                or summary.get("status") != "completed_c1p1_pharlap_fixed_angle_scan"):
            raise RuntimeError(
                f"PHaRLAP scan failed: returncode={process.returncode}")
    expected_matrix = np.array([item[1:] for item in PROFILE_ROWS], dtype=float)
    observed_matrix = np.asarray(summary["tol_matrix"], dtype=float)
    if (tuple(summary["profile_ids"]) != PROFILE_IDS
            or not np.array_equal(observed_matrix, expected_matrix)
            or tuple(summary["case_ids"]) != CASE_IDS
            or not np.array_equal(np.asarray(summary["angles_deg"]), np.array(
                [item[3] for item in CASES]))):
        raise RuntimeError("PHaRLAP scan configuration drift")

    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    records, paths, states = _load_scan(stage_dir, radius)
    c1_metrics = _c1_metrics()
    rows = []
    repeatability = {}
    for profile_id in PROFILE_IDS:
        for case_id in CASE_IDS:
            record = records[(profile_id, case_id)]
            run = record["runs"][0] if isinstance(record["runs"], list) else record["runs"]
            current = paths[(profile_id, case_id)][0]
            reference = paths[(REFERENCE_PROFILE, case_id)][0]
            metric, _ = _trajectory_metrics(reference, current, radius)
            reference_run_record = records[(REFERENCE_PROFILE, case_id)]
            reference_run = (reference_run_record["runs"][0]
                             if isinstance(reference_run_record["runs"], list)
                             else reference_run_record["runs"])
            landing = _surface_distance(
                np.asarray(current["trajectory_geocentric_lat_lon_height"])[-1],
                np.asarray(reference["trajectory_geocentric_lat_lon_height"])[-1],
                radius)
            step = _step_stats(states[(profile_id, case_id)][0])
            rows.append({
                "profile_id": profile_id,
                "case_id": case_id,
                "branch_id": CASE_META[case_id][1],
                "elevation_deg": CASE_META[case_id][3],
                "tolerance": PROFILE_ROWS[PROFILE_IDS.index(profile_id)][1],
                "min_step_km": PROFILE_ROWS[PROFILE_IDS.index(profile_id)][2],
                "max_step_km": PROFILE_ROWS[PROFILE_IDS.index(profile_id)][3],
                "termination": run["status"],
                "turning_count": int(run["turning_count"]),
                "path_point_count": int(run["path_point_count"]),
                "trajectory_rms_to_reference_km": metric["trajectory_rms_3d_km"],
                "trajectory_max_to_reference_km": metric["trajectory_max_3d_km"],
                "landing_delta_to_reference_km": landing,
                "apex_delta_to_reference_km": abs(
                    float(run["apex_km"]) - float(reference_run["apex_km"])),
                "geometric_path_delta_to_reference_km": abs(
                    float(run["geometric_path_km"])
                    - float(reference_run["geometric_path_km"])),
                "group_path_delta_to_reference_km": abs(
                    float(run["group_path_km"])
                    - float(reference_run["group_path_km"])),
                "group_delay_delta_to_reference_ms": abs(
                    float(run["tau_ms"]) - float(reference_run["tau_ms"])),
                **step,
            })
            if profile_id in {"baseline", REFERENCE_PROFILE}:
                runs = record["runs"] if isinstance(record["runs"], list) else [record["runs"]]
                repeatability[(profile_id, case_id)] = _scan_repeatability(
                    runs, paths[(profile_id, case_id)])

    by_key = {(row["profile_id"], row["case_id"]): row for row in rows}
    cross = {}
    for case_id in CASE_IDS:
        ah = _ah_path(case_id)
        baseline = paths[("baseline", case_id)][0]
        reference = paths[(REFERENCE_PROFILE, case_id)][0]
        baseline_cross, _ = _trajectory_metrics(baseline, ah, radius)
        reference_cross, _ = _trajectory_metrics(reference, ah, radius)
        neighbor = paths[(REFERENCE_NEIGHBOR, case_id)][0]
        reference_neighbor, _ = _trajectory_metrics(reference, neighbor, radius)
        neighbor_landing = _surface_distance(
            np.asarray(reference["trajectory_geocentric_lat_lon_height"])[-1],
            np.asarray(neighbor["trajectory_geocentric_lat_lon_height"])[-1], radius)
        e_base = float(c1_metrics[case_id]["trajectory_rms_3d_km"])
        if abs(baseline_cross["trajectory_rms_3d_km"] - e_base) > 1.0e-9:
            raise RuntimeError(f"baseline C1 metric drift: {case_id}")
        e_ref = float(reference_cross["trajectory_rms_3d_km"])
        cross[case_id] = {
            "E_base_km": e_base,
            "E_ref_km": e_ref,
            "closure_fraction": (e_base - e_ref) / e_base,
            "pharlap_shift_km": by_key[("baseline", case_id)][
                "trajectory_rms_to_reference_km"],
            "reference_neighbor_rms_km": reference_neighbor["trajectory_rms_3d_km"],
            "reference_neighbor_landing_km": neighbor_landing,
            "reference_converged": (
                reference_neighbor["trajectory_rms_3d_km"] <= 0.02
                and neighbor_landing <= 0.05),
        }

    figures = _plot_convergence(stage_dir, rows, _source_tag(sources))
    image_checks = [_image_check(path) for path in figures]
    all_ground = all(row["termination"] == "ground"
                     and row["turning_count"] == 1 for row in rows)
    repeats_ok = all(item["passed"] for item in repeatability.values())
    checks = {
        "c1p0_frozen": True,
        "ten_profiles_exact": len(PROFILE_IDS) == 10,
        "four_fixed_angles_only": len(CASE_IDS) == 4,
        "forty_profile_case_records": len(rows) == 40,
        "all_ground_single_turning": all_ground,
        "baseline_and_reference_repeatable": repeats_ok,
        "baseline_reproduces_formal_c1": True,
        "all_visualizations": all(item["passed"] for item in image_checks),
        "formal_c1_failure_unchanged": (
            json.loads((C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS),
    }
    status = ("completed_c1p1_pharlap_convergence_scan"
              if all(checks.values()) else "blocked_c1p1_pharlap_scan")
    csv_path = stage_dir / "c1p_1_convergence_metrics.csv"
    report_path = stage_dir / "c1p_1_convergence.json"
    _write_csv(csv_path, rows)
    _write_json(report_path, {
        "schema_version": "a4-c1p-1-convergence-v1",
        "status": status,
        "created_utc": _utc(),
        "checks": checks,
        "profile_matrix": PROFILE_ROWS,
        "cross_solver_attribution_inputs": cross,
        "repeatability": {
            f"{key[0]}::{key[1]}": value for key, value in repeatability.items()},
        "visualizations": image_checks,
        "source_sha256": sources,
        "source_set_short_hash": _source_tag(sources),
    })
    artifacts = [summary_path, log_path, csv_path, report_path, *figures]
    artifacts.extend(sorted(raw_dir.glob("*.mat")))
    _manifest(1, status, sources, artifacts, checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def _state_audit(state: dict, context: dict) -> tuple[dict, dict]:
    position = np.asarray(state["position_ecef_km"], dtype=float)
    p = np.asarray(state["phase_index_vector"], dtype=float)
    ne = np.array([
        context["provider"].get_Ne_ecef(point, context["radius"])
        for point in position
    ])
    f_hz = context["frequency_MHz"] * 1.0e6
    n2 = 1.0 - cfg.K_FP ** 2 * ne / (f_hz * f_hz)
    p2 = np.einsum("ij,ij->i", p, p)
    residual = np.abs(p2 - n2) / np.maximum(np.abs(n2), 1.0e-12)
    group = np.asarray(state["group_path_km"], dtype=float)
    geom = np.asarray(state["geometric_path_km"], dtype=float)
    phase = np.asarray(state["phase_path_km"], dtype=float)
    indep_km = np.asarray(state["indep_var_raw"], dtype=float) / 1000.0
    group_offset = group - group[0]
    geom_offset = geom - geom[0]
    indep_offset = indep_km - indep_km[0]
    group_closure = indep_offset - group_offset
    geom_closure = indep_offset - geom_offset
    metrics = {
        "point_count": int(position.shape[0]),
        "dispersion_residual_max": float(np.max(residual)),
        "dispersion_residual_rms": float(np.sqrt(np.mean(residual ** 2))),
        "group_path_monotone": bool(np.all(np.diff(group) >= -1.0e-12)),
        "geometric_path_monotone": bool(np.all(np.diff(geom) >= -1.0e-12)),
        "phase_path_monotone": bool(np.all(np.diff(phase) >= -1.0e-12)),
        "indep_var_monotone": bool(np.all(np.diff(indep_km) >= -1.0e-12)),
        "indep_group_closure_max_km": float(np.max(np.abs(group_closure))),
        "indep_geom_closure_rms_km": float(np.sqrt(np.mean(geom_closure ** 2))),
    }
    arrays = {
        "height_km": ecef_to_geocentric(position, context["radius"])[:, 2],
        "u": geom_offset / max(float(geom_offset[-1]), 1.0e-30),
        "Ne_m3": ne,
        "n2": n2,
        "p2": p2,
        "dispersion_residual": residual,
        "group_closure_km": group_closure,
        "geom_closure_km": geom_closure,
    }
    return metrics, arrays


def run_c1p2() -> str:
    _require_stage(1, {"completed_c1p1_pharlap_convergence_scan"})
    stage_dir = OUT / "c1p_2"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = _load_context(input_path)
    records, _, states = _load_scan(OUT / "c1p_1", context["radius"])
    sources = {**_base_sources(input_path),
               "c1p1_manifest": _sha256(
                   OUT / "c1p_1" / "c1p_1_freeze_manifest.json")}
    metrics = []
    arrays = {}
    for profile_id in ("baseline", REFERENCE_PROFILE):
        for case_id in CASE_IDS:
            metric, detail = _state_audit(
                states[(profile_id, case_id)][0], context)
            row = {
                "profile_id": profile_id,
                "case_id": case_id,
                "branch_id": CASE_META[case_id][1],
                **metric,
            }
            metrics.append(row)
            arrays.update({
                f"{profile_id}_{case_id}_{key}": value
                for key, value in detail.items()
            })
    csv_path = stage_dir / "c1p_2_state_metrics.csv"
    npz_path = stage_dir / "c1p_2_state_samples.npz"
    _write_csv(csv_path, metrics)
    np.savez_compressed(npz_path, **arrays)

    figure_path = stage_dir / "c1p_2_state_dispersion_audit.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, case_id in zip(axes.flat, CASE_IDS):
        for profile_id, style in (("baseline", "-"), (REFERENCE_PROFILE, "--")):
            axis.semilogy(
                arrays[f"{profile_id}_{case_id}_u"],
                np.maximum(arrays[f"{profile_id}_{case_id}_dispersion_residual"], 1e-18),
                style, label=profile_id)
        axis.set_title(f"{case_id} | {CASE_META[case_id][1]}")
        axis.set_ylabel("||p|^2 - n^2| / |n^2|")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[1, 0].set_xlabel("Normalized geometric path u")
    axes[1, 1].set_xlabel("Normalized geometric path u")
    figure.suptitle(
        "A4-C1P-2 PHaRLAP ray_state_vec dispersion audit\n"
        f"one-way | 10 MHz | B=0 | collision=0 | p=dir | source={_source_tag(sources)}")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    image = _image_check(figure_path)
    checks = {
        "c1p1_frozen": True,
        "baseline_and_reference_all_four_cases": len(metrics) == 8,
        "all_path_quantities_monotone": all(
            row[key] for row in metrics for key in (
                "group_path_monotone", "geometric_path_monotone",
                "phase_path_monotone", "indep_var_monotone")),
        "indep_var_equals_group_path_b0": all(
            row["indep_group_closure_max_km"] <= 1.0e-10 for row in metrics),
        "all_values_finite": all(
            np.isfinite(value).all() for value in arrays.values()),
        "visualization": image["passed"],
        "formal_c1_failure_unchanged": (
            json.loads((C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS),
    }
    status = ("completed_c1p2_state_semantics_audit"
              if all(checks.values()) else "blocked_c1p2_state_audit")
    report_path = stage_dir / "c1p_2_state_audit.json"
    _write_json(report_path, {
        "schema_version": "a4-c1p-2-state-audit-v1",
        "status": status,
        "created_utc": _utc(),
        "checks": checks,
        "metrics": metrics,
        "interpretation": (
            "B=0 PHaRLAP indep_var is tested against group path; "
            "dispersion is reconstructed from frozen Ne and p=dir"),
        "source_sha256": sources,
        "visualization": image,
    })
    _manifest(2, status, sources,
              [csv_path, npz_path, report_path, figure_path], checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def _central_derivative(values: np.ndarray, path: np.ndarray) -> np.ndarray:
    result = np.full_like(values, np.nan, dtype=float)
    denominator = path[2:] - path[:-2]
    valid = denominator > 1.0e-12
    result[1:-1][valid] = (
        values[2:][valid] - values[:-2][valid]) / denominator[valid, None]
    return result


def _angle_deg(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    first_norm = np.linalg.norm(first, axis=1)
    second_norm = np.linalg.norm(second, axis=1)
    valid = (first_norm > 1.0e-30) & (second_norm > 1.0e-30)
    result = np.full(first.shape[0], np.nan)
    cosine = np.einsum("ij,ij->i", first[valid], second[valid])
    cosine /= first_norm[valid] * second_norm[valid]
    result[valid] = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return result


def _rhs_audit(state: dict, context: dict) -> tuple[dict, dict]:
    position = np.asarray(state["position_ecef_km"], dtype=float)
    p = np.asarray(state["phase_index_vector"], dtype=float)
    path = np.asarray(state["geometric_path_km"], dtype=float)
    path -= path[0]
    dr_numeric = _central_derivative(position, path)
    dp_numeric = _central_derivative(p, path)
    model = _model(context)
    sample_count = min(81, max(position.shape[0] - 2, 0))
    indices = np.unique(np.linspace(
        1, position.shape[0] - 2, sample_count, dtype=int))
    rhs_dr = np.full((indices.size, 3), np.nan)
    rhs_dp = np.full((indices.size, 3), np.nan)
    valid = np.ones(indices.size, dtype=bool)
    for output_index, state_index in enumerate(indices):
        physical_state = np.r_[
            position[state_index], model.k0 * p[state_index], 0.0]
        try:
            rhs = ah3d_ecef_rhs(path[state_index], physical_state, model)
            rhs_dr[output_index] = rhs[:3]
            rhs_dp[output_index] = rhs[3:6] / model.k0
        except (ValueError, FloatingPointError):
            valid[output_index] = False
    numeric_dr = dr_numeric[indices]
    numeric_dp = dp_numeric[indices]
    tangent_angle = _angle_deg(numeric_dr, rhs_dr)
    rhs_dp_norm = np.linalg.norm(rhs_dp, axis=1)
    numeric_dp_norm = np.linalg.norm(numeric_dp, axis=1)
    maximum = float(np.nanmax(rhs_dp_norm)) if np.any(np.isfinite(rhs_dp_norm)) else 0.0
    active = valid & np.isfinite(tangent_angle) & np.isfinite(rhs_dp_norm)
    active &= rhs_dp_norm > max(maximum * 1.0e-4, 1.0e-12)
    dp_relative = np.full(indices.size, np.nan)
    dp_direction = np.full(indices.size, np.nan)
    dp_relative[active] = (
        np.linalg.norm(numeric_dp[active] - rhs_dp[active], axis=1)
        / np.maximum(rhs_dp_norm[active], 1.0e-30))
    dp_direction[active] = _angle_deg(numeric_dp[active], rhs_dp[active])
    finite_tangent = tangent_angle[valid & np.isfinite(tangent_angle)]
    finite_dp = dp_relative[np.isfinite(dp_relative)]
    finite_direction = dp_direction[np.isfinite(dp_direction)]
    if finite_tangent.size == 0 or finite_dp.size == 0:
        raise RuntimeError("insufficient finite local RHS samples")
    geocentric = ecef_to_geocentric(position[indices], context["radius"])
    metrics = {
        "sample_count": int(indices.size),
        "valid_rhs_count": int(np.count_nonzero(valid)),
        "active_bending_count": int(finite_dp.size),
        "tangent_median_deg": float(np.median(finite_tangent)),
        "tangent_p95_deg": float(np.percentile(finite_tangent, 95)),
        "dp_ds_relative_median": float(np.median(finite_dp)),
        "dp_ds_relative_p95": float(np.percentile(finite_dp, 95)),
        "dp_ds_direction_median_deg": float(np.median(finite_direction)),
        "dp_ds_direction_p95_deg": float(np.percentile(finite_direction, 95)),
    }
    arrays = {
        "u": path[indices] / max(float(path[-1]), 1.0e-30),
        "height_km": geocentric[:, 2],
        "tangent_angle_deg": tangent_angle,
        "dp_ds_relative": dp_relative,
        "dp_ds_direction_deg": dp_direction,
        "rhs_dp_norm_per_km": rhs_dp_norm,
        "numeric_dp_norm_per_km": numeric_dp_norm,
        "active": active.astype(np.int8),
    }
    return metrics, arrays


def run_c1p3() -> str:
    _require_stage(2, {"completed_c1p2_state_semantics_audit"})
    stage_dir = OUT / "c1p_3"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = _load_context(input_path)
    _, _, states = _load_scan(OUT / "c1p_1", context["radius"])
    sources = {**_base_sources(input_path),
               "c1p2_manifest": _sha256(
                   OUT / "c1p_2" / "c1p_2_freeze_manifest.json")}
    metrics = []
    arrays = {}
    for profile_id in ("baseline", REFERENCE_PROFILE):
        for case_id in CASE_IDS:
            metric, detail = _rhs_audit(states[(profile_id, case_id)][0], context)
            metrics.append({
                "profile_id": profile_id,
                "case_id": case_id,
                "branch_id": CASE_META[case_id][1],
                **metric,
            })
            arrays.update({
                f"{profile_id}_{case_id}_{key}": value
                for key, value in detail.items()
            })
    strict = {row["case_id"]: row for row in metrics
              if row["profile_id"] == REFERENCE_PROFILE}
    low_rows = [strict[case_id] for case_id in CASE_IDS
                if CASE_META[case_id][1] == "low_ray"]
    high_rows = [strict[case_id] for case_id in CASE_IDS
                 if CASE_META[case_id][1] == "high_ray"]
    low_tangent = max(row["tangent_p95_deg"] for row in low_rows)
    low_dp = max(row["dp_ds_relative_median"] for row in low_rows)
    per_high_candidate = {
        row["case_id"]: (
            row["active_bending_count"] >= 20
            and (row["tangent_p95_deg"] > max(0.01, 10.0 * low_tangent)
                 or row["dp_ds_relative_median"] > max(0.1, 10.0 * low_dp)))
        for row in high_rows
    }
    candidate = all(per_high_candidate.values())

    csv_path = stage_dir / "c1p_3_local_rhs_metrics.csv"
    npz_path = stage_dir / "c1p_3_local_rhs_samples.npz"
    _write_csv(csv_path, metrics)
    np.savez_compressed(npz_path, **arrays)
    figure_path = stage_dir / "c1p_3_local_rhs_comparison.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9))
    colors = {
        "low_pharlap_anchor": "#1368aa",
        "low_ah3d_anchor": "#46a3c7",
        "high_pharlap_anchor": "#b22222",
        "high_ah3d_anchor": "#e0704f",
    }
    for case_id in CASE_IDS:
        key = f"{REFERENCE_PROFILE}_{case_id}"
        u = arrays[f"{key}_u"]
        height = arrays[f"{key}_height_km"]
        tangent = arrays[f"{key}_tangent_angle_deg"]
        dp_rel = arrays[f"{key}_dp_ds_relative"]
        axes[0, 0].semilogy(u, np.maximum(tangent, 1e-12),
                           label=case_id, color=colors[case_id])
        axes[0, 1].semilogy(u, np.maximum(dp_rel, 1e-12),
                           label=case_id, color=colors[case_id])
        axes[1, 0].semilogy(height, np.maximum(tangent, 1e-12),
                           color=colors[case_id])
        axes[1, 1].semilogy(height, np.maximum(dp_rel, 1e-12),
                           color=colors[case_id])
    axes[0, 0].set(xlabel="Normalized path u", ylabel="Tangent angle error (deg)",
                   title="dr/ds direction versus AH3D RHS")
    axes[0, 1].set(xlabel="Normalized path u", ylabel="dp/ds relative error",
                   title="dp/ds magnitude/vector residual")
    axes[1, 0].set(xlabel="Physical height (km)", ylabel="Tangent angle error (deg)",
                   title="Tangent residual by height")
    axes[1, 1].set(xlabel="Physical height (km)", ylabel="dp/ds relative error",
                   title="Bending-region residual by height")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7)
    figure.suptitle(
        "A4-C1P-3 PHaRLAP secants versus AH3D B=0 Hamiltonian RHS\n"
        f"strict PHaRLAP profile | fixed angles | source={_source_tag(sources)}")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    image = _image_check(figure_path)
    checks = {
        "c1p2_frozen": True,
        "baseline_and_reference_four_cases": len(metrics) == 8,
        "all_rhs_samples_finite_enough": all(
            row["valid_rhs_count"] >= 20 and row["active_bending_count"] >= 20
            for row in metrics),
        "candidate_rule_applied_to_both_high_controls": (
            len(per_high_candidate) == 2),
        "visualization": image["passed"],
        "formal_c1_failure_unchanged": (
            json.loads((C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS),
    }
    status = ("completed_c1p3_local_rhs_audit"
              if all(checks.values()) else "blocked_c1p3_local_rhs_audit")
    report_path = stage_dir / "c1p_3_local_rhs_audit.json"
    _write_json(report_path, {
        "schema_version": "a4-c1p-3-local-rhs-v1",
        "status": status,
        "created_utc": _utc(),
        "checks": checks,
        "metrics": metrics,
        "low_control_thresholds": {
            "tangent_p95_deg_max": low_tangent,
            "dp_ds_relative_median_max": low_dp,
        },
        "high_case_candidate_flags": per_high_candidate,
        "state_or_rhs_semantics_candidate": candidate,
        "source_sha256": sources,
        "visualization": image,
    })
    _manifest(3, status, sources,
              [csv_path, npz_path, report_path, figure_path], checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def _write_final_report(path: Path, status: str, integration: str,
                        final_attribution: str, high: dict,
                        rhs_candidate: bool) -> None:
    lines = [
        "A4-C1P PHaRLAP侧数值归因实施报告",
        "",
        f"日期：{_utc()}",
        f"机器状态：{status}",
        f"最终归因：{final_attribution}",
        "",
        "一、实施范围",
        "",
        "1. 使用冻结Chapman 1 km网格、B=0、collision=0、10 MHz和四个C1公共角。",
        "2. 完成10组PHaRLAP tolerance/min/max-step扫描、状态语义审计和局部RHS比较。",
        "3. 未重新求根，未修改AH3D、PHaRLAP MEX、原C1/C1R产物或2 km门槛。",
        "",
        "二、积分归因",
        "",
        f"PHaRLAP积分分类：{integration}。",
    ]
    for case_id, item in high.items():
        lines.append(
            f"{case_id}: E_base={item['E_base_km']:.9g} km, "
            f"E_ref={item['E_ref_km']:.9g} km, "
            f"pharlap_shift={item['pharlap_shift_km']:.9g} km, "
            f"closure={100.0 * item['closure_fraction']:.6g}%, "
            f"reference_converged={item['reference_converged']}。")
    lines.extend([
        "",
        "三、状态与局部RHS",
        "",
        "1. PHaRLAP ray_state_vec按公开语义转换：ECEF位置为m，dir向量模为相折射率。",
        "2. B=0下indep_var与group_path闭合；色散、路径单调性和局部RHS统计均已冻结。",
        f"3. state_or_rhs_semantics_candidate={rhs_candidate}。",
        "",
        "四、结论与门控",
        "",
        f"1. 当前差异源分类为{final_attribution}。",
        "2. 原A4-C1状态仍为failed_c1_chapman_high_low_external_comparison。",
        "3. C1P完成不自动开放C2；后续由用户审阅本报告后决定修复、第三方验证或有限诊断。",
        "4. 本结论仅适用于B=0平滑Chapman单程高/低射线，不外推到O/X、复杂3D扰动、散射、PE或全波。",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def run_c1p4() -> str:
    _require_stage(3, {"completed_c1p3_local_rhs_audit"})
    stage_dir = OUT / "c1p_4"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = {**_base_sources(input_path),
               **{f"c1p{stage}_manifest": _sha256(
                   OUT / f"c1p_{stage}" / f"c1p_{stage}_freeze_manifest.json")
                  for stage in range(4)}}
    convergence = json.loads(
        (OUT / "c1p_1" / "c1p_1_convergence.json").read_text(encoding="utf-8"))
    rhs = json.loads(
        (OUT / "c1p_3" / "c1p_3_local_rhs_audit.json").read_text(encoding="utf-8"))
    cross = convergence["cross_solver_attribution_inputs"]
    high = {case_id: cross[case_id] for case_id in CASE_IDS
            if CASE_META[case_id][1] == "high_ray"}
    not_primary = all(
        item["pharlap_shift_km"] <= 0.1
        and abs(item["closure_fraction"]) <= 0.05
        for item in high.values())
    primary = (
        all(item["reference_converged"] for item in high.values())
        and any(item["closure_fraction"] >= 0.5 for item in high.values()))
    if primary:
        integration = "pharlap_integration_primary_candidate"
    elif not_primary:
        integration = "pharlap_integration_not_primary"
    else:
        integration = "pharlap_integration_contributor"
    rhs_candidate = bool(rhs["state_or_rhs_semantics_candidate"])
    if integration == "pharlap_integration_primary_candidate":
        final_attribution = integration
    elif rhs_candidate:
        final_attribution = "state_or_rhs_semantics_candidate"
    elif integration == "pharlap_integration_contributor":
        final_attribution = integration
    else:
        final_attribution = "software_difference_or_unresolved"
    status = f"completed_c1p4_{final_attribution}"

    figure_path = stage_dir / "c1p_4_attribution_matrix.png"
    labels = (
        "PHaRLAP integration\nprimary",
        "PHaRLAP integration\ncontributor",
        "State/RHS\ncandidate",
        "Software difference\nor unresolved",
        "Formal C1\nremains failed",
        "C2 automatic\nadmission",
    )
    values = np.array([
        integration == "pharlap_integration_primary_candidate",
        integration == "pharlap_integration_contributor",
        rhs_candidate,
        final_attribution == "software_difference_or_unresolved",
        1,
        0,
    ], dtype=int)
    figure, axis = plt.subplots(figsize=(12, 5.5))
    colors = ["#187b3f" if value else "#c7c7c7" for value in values]
    colors[-2] = "#b22222"
    colors[-1] = "#b22222"
    axis.bar(np.arange(len(labels)), values, color=colors)
    axis.set_xticks(np.arange(len(labels)), labels)
    axis.set_ylim(0, 1.15)
    axis.set_yticks((0, 1), ("NO", "YES"))
    axis.grid(axis="y", alpha=0.2)
    axis.set_title(
        "A4-C1P-4 PHaRLAP-side attribution and manual decision gate\n"
        f"result={final_attribution} | source={_source_tag(sources)}")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    image = _image_check(figure_path)
    checks = {
        "c1p0_to_c1p3_frozen": True,
        "two_high_cases_classified_independently": len(high) == 2,
        "attribution_rule_executed": final_attribution in {
            "pharlap_integration_primary_candidate",
            "pharlap_integration_contributor",
            "state_or_rhs_semantics_candidate",
            "software_difference_or_unresolved",
        },
        "formal_c1_failure_preserved": (
            json.loads((C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS),
        "c2_not_automatically_admitted": True,
        "visualization": image["passed"],
    }
    machine_path = stage_dir / "c1p_4_attribution.json"
    report_path = ROOT / "word" / "a4_c1p_pharlap_numerical_attribution_report.txt"
    _write_json(machine_path, {
        "schema_version": "a4-c1p-4-attribution-v1",
        "status": status,
        "created_utc": _utc(),
        "checks": checks,
        "integration_attribution": integration,
        "state_or_rhs_semantics_candidate": rhs_candidate,
        "final_attribution": final_attribution,
        "high_case_metrics": high,
        "formal_c1_status": FINAL_C1_STATUS,
        "c2_gate": "manual_decision_required_formal_c1_still_failed",
        "source_sha256": sources,
        "visualization": image,
    })
    _write_final_report(
        report_path, status, integration, final_attribution, high, rhs_candidate)
    _manifest(4, status, sources,
              [machine_path, report_path, figure_path], checks)
    return status


RUNNERS = {
    0: run_c1p0,
    1: run_c1p1,
    2: run_c1p2,
    3: run_c1p3,
    4: run_c1p4,
}


def main() -> int:
    stages = range(5)
    if len(sys.argv) > 1:
        if len(sys.argv) != 2 or not sys.argv[1].startswith("--stage="):
            raise SystemExit("usage: _run_a4_c1p_pharlap_attribution.py [--stage=0|1|2|3|4|all]")
        value = sys.argv[1].split("=", 1)[1]
        stages = range(5) if value == "all" else (int(value),)
    for stage in stages:
        print(f"A4-C1P-{stage}: {RUNNERS[stage]()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
