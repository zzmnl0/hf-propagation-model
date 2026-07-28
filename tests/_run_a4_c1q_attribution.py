"""Run A4-C1Q PHaRLAP-truth along-track attribution diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
import models.ray_tracer_ah_3d_ecef as ahcore
from models.dispersion_ah import C_KM_S
from models.ne_provider_geocentric import GeocentricNeProvider3D
from tests import _run_a4_c1_prep3 as prep3
from tests._run_a4_c1_validation import (
    CASES, _image_check, _refined_input, _resample, _sha256, _trajectory_metrics,
)
from tests._run_a4_c1p_pharlap_attribution import (
    C1, FINAL_C1_STATUS, _verify_any_manifest,
)
from utils import ecef_to_geocentric, geocentric_to_ecef


OUT = C1 / "c1q"
MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SOURCE = ROOT / "tests" / "matlab" / "a4_c1q_pharlap_state_export.m"
MATLAB_SCRIPT = MATLAB_DIR / MATLAB_SOURCE.name
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
AH_SOLVER = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
AH_DISPERSION = ROOT / "models" / "dispersion_ah.py"
AH_PROVIDER = ROOT / "models" / "ne_provider_geocentric.py"
CONFIG = ROOT / "config.py"
PLAN = ROOT / "word" / "a4_c1q_pharlap_truth_along_track_attribution_and_repair_plan.txt"
CASE_IDS = tuple(item[0] for item in CASES)
CASE_META = {item[0]: item for item in CASES}
STRICT_TOL = (1e-12, 0.001, 1.0)
FINAL_STATES = {
    "identified_and_repaired_formula_or_implementation",
    "identified_but_not_safely_repairable",
    "pharlap_effective_dynamics_not_recoverable",
    "multiple_unresolved_contributors",
}


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value) -> None:
    def json_safe(item):
        if isinstance(item, dict):
            return {key: json_safe(entry) for key, entry in item.items()}
        if isinstance(item, (list, tuple)):
            return [json_safe(entry) for entry in item]
        if isinstance(item, np.ndarray):
            return json_safe(item.tolist())
        if isinstance(item, np.generic):
            return json_safe(item.item())
        if isinstance(item, float) and not np.isfinite(item):
            return None
        return item

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        json_safe(value), indent=2, ensure_ascii=True, allow_nan=False,
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


def _sources(input_path: Path) -> dict:
    return {
        "runner": _sha256(Path(__file__)),
        "plan": _sha256(PLAN),
        "matlab_source": _sha256(MATLAB_SOURCE),
        "refined_input": _sha256(input_path),
        "pharlap_mex": _sha256(MEX),
        "ah_solver": _sha256(AH_SOLVER),
        "ah_dispersion": _sha256(AH_DISPERSION),
        "ah_provider": _sha256(AH_PROVIDER),
        "config": _sha256(CONFIG),
    }


def _source_tag(sources: dict) -> str:
    return hashlib.sha256(json.dumps(
        sources, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()[:12]


def _manifest(stage: int, status: str, sources: dict, artifacts: list[Path],
              checks: dict) -> None:
    stage_dir = OUT / f"c1q_{stage}"
    sources = {**sources, "stage_runner": _sha256(Path(__file__))}
    _write_json(stage_dir / f"c1q_{stage}_freeze_manifest.json", {
        "schema_version": f"a4-c1q-{stage}-freeze-v1",
        "status": status,
        "created_utc": _utc(),
        "semantics": (
            "PHaRLAP_truth_for_project_validation; same_angle_forward_trace; "
            "one_way; 10 MHz; B=0; collision=0; no root solving; no empirical fit"
        ),
        "checks": checks,
        "source_sha256": sources,
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in artifacts
        ],
    })


def _require(stage: int, allowed: set[str]) -> dict:
    return _verify_any_manifest(
        OUT / f"c1q_{stage}" / f"c1q_{stage}_freeze_manifest.json", allowed)


def _image(path: Path) -> dict:
    result = _image_check(path)
    if not result["passed"]:
        raise RuntimeError(f"invalid C1Q visualization: {path}")
    return result


def _h5_vector(group, name: str) -> np.ndarray:
    return np.asarray(group[name], dtype=float).reshape(-1)


def _h5_struct(source, name: str):
    node = source[name]
    if isinstance(node, h5py.Group):
        return node
    references = node[()].reshape(-1)
    if references.size != 1:
        raise RuntimeError(f"{name} is not a scalar MATLAB structure")
    return source[references[0]]


def _load_pharlap(path: Path, radius_km: float) -> dict:
    with h5py.File(path, "r") as source:
        path_group = _h5_struct(source, "ray_path_data")
        state_group = _h5_struct(source, "ray_state_vec")
        lat = _h5_vector(path_group, "lat")
        lon = _h5_vector(path_group, "lon")
        height = _h5_vector(path_group, "height")
        position = np.column_stack([
            _h5_vector(state_group, key) for key in ("pos_x", "pos_y", "pos_z")
        ]) / 1000.0
        phase_index = np.column_stack([
            _h5_vector(state_group, key) for key in ("dir_x", "dir_y", "dir_z")
        ])
        result = {
            "position_ecef_km": position,
            "path_position_ecef_km": geocentric_to_ecef(
                lat, lon, height, radius_km),
            "phase_index_vector": phase_index,
            "group_path_km": _h5_vector(state_group, "group_path") / 1000.0,
            "geom_path_km": _h5_vector(state_group, "geometrical_path") / 1000.0,
            "phase_path_km": _h5_vector(state_group, "phase_path") / 1000.0,
            "indep_var_km": _h5_vector(state_group, "indep_var") / 1000.0,
            "ode_step_km": _h5_vector(state_group, "ODE_step_size") / 1000.0,
            "lat_deg": lat,
            "lon_deg": lon,
            "height_km": height,
            "electron_density_m3": _h5_vector(path_group, "electron_density") * 1e6,
            "refractive_index": _h5_vector(path_group, "refractive_index"),
            "group_refractive_index": _h5_vector(
                path_group, "group_refractive_index"),
            "path_group_km": _h5_vector(path_group, "group_range"),
            "path_geom_km": _h5_vector(path_group, "geometric_distance"),
        }
    sizes = {np.asarray(value).shape[0] for value in result.values()
             if np.asarray(value).ndim > 0}
    if len(sizes) != 1:
        raise RuntimeError(f"PHaRLAP path/state length mismatch: {path}")
    for key in ("group_path_km", "geom_path_km", "phase_path_km",
                "indep_var_km", "path_group_km", "path_geom_km"):
        result[key] = np.asarray(result[key], dtype=float)
        result[key] -= result[key][0]
    return result


def _c1q0_records() -> dict:
    summary = json.loads(
        (OUT / "c1q_0" / "c1q_0_pharlap_export_summary.json").read_text(
            encoding="utf-8"))
    return {item["case_id"]: item for item in summary["records"]}


def _ph_cases(radius: float) -> dict:
    return {
        case_id: _load_pharlap(Path(record["raw_mat"]), radius)
        for case_id, record in _c1q0_records().items()
    }


def _ah_case(case_id: str, context: dict) -> dict:
    raw = C1 / "same_angle" / "raw" / f"ah3d_{case_id}_raw.npz"
    with np.load(raw, allow_pickle=False) as source:
        position = np.asarray(source["repeat1_trajectory_ecef_km"], dtype=float)
        wavevector = np.asarray(
            source["repeat1_wavevector_ecef_km^-1"], dtype=float)
        geom = np.asarray(
            source["repeat1_trajectory_geom_path_km"], dtype=float)
    model = prep3._model(context)
    ne = np.array([
        context["provider"].get_Ne_ecef(point, context["radius"])
        for point in position
    ])
    f_hz = context["frequency_MHz"] * 1e6
    n2 = 1.0 - cfg.K_FP**2 * ne / f_hz**2
    n = np.sqrt(np.maximum(n2, 1e-30))
    group = cumulative_trapezoid(1.0 / n, geom, initial=0.0)
    geocentric = ecef_to_geocentric(position, context["radius"])
    return {
        "position_ecef_km": position,
        "phase_index_vector": wavevector / model.k0,
        "wavevector_ecef_km^-1": wavevector,
        "geom_path_km": geom,
        "group_path_km": group,
        "height_km": geocentric[:, 2],
        "lat_deg": geocentric[:, 0],
        "lon_deg": geocentric[:, 1],
        "electron_density_m3": ne,
        "n2": n2,
    }


def _link_basis(context: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tx = geocentric_to_ecef(*context["tx"], context["radius"])
    rx = geocentric_to_ecef(*context["rx"], context["radius"])
    tx_hat = tx / np.linalg.norm(tx)
    rx_hat = rx / np.linalg.norm(rx)
    normal = np.cross(tx_hat, rx_hat)
    normal /= np.linalg.norm(normal)
    tangent = np.cross(normal, tx_hat)
    return tx_hat, tangent, normal


def _xyz(position: np.ndarray, context: dict) -> np.ndarray:
    tx_hat, tangent, normal = _link_basis(context)
    unit = position / np.linalg.norm(position, axis=1)[:, None]
    x = context["radius"] * np.arctan2(unit @ tangent, unit @ tx_hat)
    y = context["radius"] * np.arcsin(np.clip(unit @ normal, -1.0, 1.0))
    z = np.linalg.norm(position, axis=1) - context["radius"]
    return np.column_stack((x, y, z))


def _path_dict(case: dict) -> dict:
    position = case["position_ecef_km"]
    geocentric = np.column_stack((
        case["lat_deg"], case["lon_deg"], case["height_km"]))
    return {
        "trajectory_ecef_km": position,
        "trajectory_geocentric_lat_lon_height": geocentric,
        "trajectory_physical_height_km": case["height_km"],
        "trajectory_geom_path_km": case["geom_path_km"],
    }


def _local_poly_derivative(x: np.ndarray, values: np.ndarray,
                           window: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    vector = values.ndim == 2
    source = values if vector else values[:, None]
    output = np.full_like(source, np.nan)
    half = window // 2
    for index in range(half, len(x) - half):
        sl = slice(index - half, index + half + 1)
        local_x = x[sl] - x[index]
        if np.ptp(local_x) <= 1e-12:
            continue
        for axis in range(source.shape[1]):
            coefficient = np.polynomial.polynomial.polyfit(
                local_x, source[sl, axis], min(3, window - 1))
            output[index, axis] = coefficient[1]
    return output if vector else output[:, 0]


def _rms(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.sqrt(np.mean(finite**2))) if finite.size else float("nan")


def _max(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.max(np.abs(finite))) if finite.size else float("nan")


def _deploy_and_run_matlab(input_path: Path, stage_dir: Path) -> subprocess.CompletedProcess:
    MATLAB_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MATLAB_SOURCE, MATLAB_SCRIPT)
    if _sha256(MATLAB_SOURCE) != _sha256(MATLAB_SCRIPT):
        raise RuntimeError("deployed C1Q MATLAB script hash mismatch")
    raw = stage_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    summary = stage_dir / "c1q_0_pharlap_export_summary.json"
    quote = lambda value: str(value).replace("'", "''")
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); a4_c1q_pharlap_state_export("
        f"'{quote(input_path)}','{quote(raw)}','{quote(summary)}')")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run(
        [str(MATLAB), "-batch", expression], env=env,
        capture_output=True, text=True, timeout=1800, errors="replace")
    (stage_dir / "c1q_0_matlab.log").write_text(
        (result.stdout or "") + "\n" + (result.stderr or ""), encoding="utf-8")
    return result


def run_c1q0() -> str:
    stage_dir = OUT / "c1q_0"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, refined = _refined_input()
    sources = _sources(input_path)
    c1p_statuses = {}
    allowed = {
        0: {"passed_c1p0_frozen_admission"},
        1: {"completed_c1p1_pharlap_convergence_scan"},
        2: {"completed_c1p2_state_semantics_audit"},
        3: {"completed_c1p3_local_rhs_audit"},
        4: {"completed_c1p4_software_difference_or_unresolved"},
    }
    for stage, statuses in allowed.items():
        manifest = _verify_any_manifest(
            C1 / "c1p" / f"c1p_{stage}" / f"c1p_{stage}_freeze_manifest.json",
            statuses)
        c1p_statuses[str(stage)] = manifest["status"]
    final = json.loads(
        (C1 / "final" / "a4_c1_validation.json").read_text(encoding="utf-8"))
    result = _deploy_and_run_matlab(input_path, stage_dir)
    summary_path = stage_dir / "c1q_0_pharlap_export_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    context = prep3._load_context(input_path)
    records = {item["case_id"]: item for item in summary.get("records", [])}
    comparisons = {}
    c1p_summary = json.loads(
        (C1 / "c1p" / "c1p_1" / "c1p_1_pharlap_scan_summary.json").read_text(
            encoding="utf-8"))
    c1p_records = {
        (item["profile_id"], item["case_id"]): item
        for item in c1p_summary["records"]
    }
    for case_id in CASE_IDS:
        current = _load_pharlap(Path(records[case_id]["raw_mat"]), context["radius"])
        reference_path = Path(c1p_records[("reference", case_id)]["raw_mat"])
        from tests._run_a4_c1p_pharlap_attribution import _load_pharlap_raw
        reference, _ = _load_pharlap_raw(reference_path, context["radius"])
        metric, _ = _trajectory_metrics(reference[0], _path_dict(current),
                                         context["radius"])
        comparisons[case_id] = metric["trajectory_rms_3d_km"]
    checks = {
        "c1p_all_stages_frozen": len(c1p_statuses) == 5,
        "formal_c1_failure_preserved": final["status"] == FINAL_C1_STATUS,
        "refined_input_frozen": _sha256(input_path) == refined["sha256"],
        "strict_profile_exact": tuple(summary.get("tol", [])) == STRICT_TOL,
        "four_fixed_angles": set(records) == set(CASE_IDS),
        "all_ground": all(int(item["ray_label"]) == 1 for item in records.values()),
        "matlab_completed": result.returncode == 0 and summary.get("status")
            == "completed_c1q0_strict_pharlap_state_export",
        "matches_c1p_strict_reference": all(value <= 1e-5
                                             for value in comparisons.values()),
    }
    figure_path = stage_dir / "c1q_0_frozen_admission.png"
    figure, axis = plt.subplots(figsize=(11, 5.5))
    labels = list(checks)
    values = [int(checks[key]) for key in labels]
    axis.barh(np.arange(len(labels)), values,
              color=["#187b3f" if value else "#b22222" for value in values])
    axis.set_yticks(np.arange(len(labels)), [item.replace("_", " ") for item in labels])
    axis.set_xticks((0, 1), ("FAIL", "PASS"))
    axis.set_xlim(0, 1.05)
    axis.set_title(f"A4-C1Q-0 frozen admission | source={_source_tag(sources)}")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    image = _image(figure_path)
    checks["visualization"] = image["passed"]
    status = ("passed_c1q0_frozen_admission"
              if all(checks.values()) else "blocked_c1q0_frozen_admission")
    report = stage_dir / "c1q_0_admission.json"
    _write_json(report, {
        "schema_version": "a4-c1q-0-admission-v1",
        "status": status,
        "checks": checks,
        "c1p_statuses": c1p_statuses,
        "strict_reference_rms_km": comparisons,
        "source_sha256": sources,
    })
    artifacts = [summary_path, stage_dir / "c1q_0_matlab.log", report, figure_path]
    artifacts.extend(Path(item["raw_mat"]) for item in records.values())
    _manifest(0, status, sources, artifacts, checks)
    if not status.startswith("passed"):
        raise RuntimeError(status)
    return status


def _height_branch_x_rms(ph: dict, ah: dict, context: dict) -> tuple[float, float]:
    result = []
    for source in (ph, ah):
        source["xyz"] = _xyz(source["position_ecef_km"], context)
    for descending in (False, True):
        arrays = []
        for source in (ph, ah):
            apex = int(np.argmax(source["height_km"]))
            if descending:
                height = source["height_km"][apex:][::-1]
                x = source["xyz"][apex:, 0][::-1]
            else:
                height = source["height_km"][:apex + 1]
                x = source["xyz"][:apex + 1, 0]
            unique, index = np.unique(height, return_index=True)
            arrays.append((unique, x[index]))
        lower = max(item[0][0] for item in arrays)
        upper = min(item[0][-1] for item in arrays)
        sample = np.linspace(lower, upper, 401)
        difference = np.interp(sample, arrays[1][0], arrays[1][1]) - np.interp(
            sample, arrays[0][0], arrays[0][1])
        result.append(_rms(difference))
    return result[0], result[1]


def run_c1q1() -> str:
    previous = _require(0, {"passed_c1q0_frozen_admission"})
    stage_dir = OUT / "c1q_1"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    ph_cases = _ph_cases(context["radius"])
    rows, arrays = [], {}
    phases = ("launch", "ascending", "apex", "descending", "ground_tail")
    for case_id, branch, anchor, angle in CASES:
        ph = ph_cases[case_id]
        ah = _ah_case(case_id, context)
        u = np.linspace(0.0, 1.0, 1001)
        ph_pos = _resample(ph["geom_path_km"], ph["position_ecef_km"],
                           u * ph["geom_path_km"][-1])
        ah_pos = _resample(ah["geom_path_km"], ah["position_ecef_km"],
                           u * ah["geom_path_km"][-1])
        delta = _xyz(ah_pos, context) - _xyz(ph_pos, context)
        common_end = min(ph["geom_path_km"][-1], ah["geom_path_km"][-1])
        common_s = np.linspace(0.0, common_end, 2001)
        ph_common = _xyz(_resample(
            ph["geom_path_km"], ph["position_ecef_km"], common_s), context)
        ah_common = _xyz(_resample(
            ah["geom_path_km"], ah["position_ecef_km"], common_s), context)
        common_dx = ah_common[:, 0] - ph_common[:, 0]
        ph_tail = _xyz(ph["position_ecef_km"][[-1]], context)[0, 0] - ph_common[-1, 0]
        ah_tail = _xyz(ah["position_ecef_km"][[-1]], context)[0, 0] - ah_common[-1, 0]
        endpoint_dx = delta[-1, 0]
        transport_dx = common_dx[-1] - common_dx[0]
        tail_dx = ah_tail - ph_tail
        closure = endpoint_dx - (transport_dx + tail_dx)
        apex_u = 0.5 * (
            ph["geom_path_km"][np.argmax(ph["height_km"])] / ph["geom_path_km"][-1]
            + ah["geom_path_km"][np.argmax(ah["height_km"])] / ah["geom_path_km"][-1])
        boundaries_u = np.array([
            0.0, 0.05, max(0.06, apex_u - 0.05),
            min(0.94, apex_u + 0.05), 0.95, 1.0])
        boundaries_u = np.maximum.accumulate(boundaries_u)
        boundary_index = np.rint(boundaries_u * (len(common_s) - 1)).astype(int)
        contributions = {}
        for index, phase in enumerate(phases):
            contribution = (
                common_dx[boundary_index[index + 1]]
                - common_dx[boundary_index[index]])
            if phase == "ground_tail":
                contribution += tail_dx
            contributions[phase] = float(contribution)
        group_end = min(ph["group_path_km"][-1], ah["group_path_km"][-1])
        group = np.linspace(0.0, group_end, 1001)
        group_delta = _xyz(_resample(
            ah["group_path_km"], ah["position_ecef_km"], group), context) - _xyz(
            _resample(ph["group_path_km"], ph["position_ecef_km"], group), context)
        asc_rms, desc_rms = _height_branch_x_rms(ph, ah, context)
        row = {
            "case_id": case_id, "branch_id": branch, "anchor_solver": anchor,
            "elevation_deg": angle,
            "normalized_x_rms_km": _rms(delta[:, 0]),
            "normalized_y_rms_km": _rms(delta[:, 1]),
            "normalized_z_rms_km": _rms(delta[:, 2]),
            "landing_dx_km": float(endpoint_dx),
            "common_transport_dx_km": float(transport_dx),
            "unequal_tail_dx_km": float(tail_dx),
            "budget_closure_km": float(closure),
            "group_x_rms_km": _rms(group_delta[:, 0]),
            "group_y_rms_km": _rms(group_delta[:, 1]),
            "group_z_rms_km": _rms(group_delta[:, 2]),
            "same_height_ascending_x_rms_km": asc_rms,
            "same_height_descending_x_rms_km": desc_rms,
            **{f"{phase}_dx_contribution_km": contributions[phase]
               for phase in phases},
        }
        rows.append(row)
        arrays.update({
            f"{case_id}_u": u,
            f"{case_id}_normalized_delta_xyz_km": delta,
            f"{case_id}_common_s_km": common_s,
            f"{case_id}_common_dx_km": common_dx,
            f"{case_id}_group_km": group,
            f"{case_id}_group_delta_xyz_km": group_delta,
        })
    csv_path = stage_dir / "c1q_1_transport_metrics.csv"
    npz_path = stage_dir / "c1q_1_transport_samples.npz"
    _write_csv(csv_path, rows)
    np.savez_compressed(npz_path, **arrays)
    residual_figure = stage_dir / "c1q_1_xyz_residuals.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, row in zip(axes.flat, rows):
        case_id = row["case_id"]
        delta = arrays[f"{case_id}_normalized_delta_xyz_km"]
        for index, label, color in (
                (0, "x along-link", "#1368aa"), (1, "y cross-track", "#16836b"),
                (2, "z height", "#c43b3b")):
            axis.plot(arrays[f"{case_id}_u"], delta[:, index],
                      label=label, color=color)
        axis.set_title(f"{row['branch_id']} | {row['anchor_solver']} anchor")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[1, 0].set_xlabel("Normalized geometric path u")
    axes[1, 1].set_xlabel("Normalized geometric path u")
    axes[0, 0].set_ylabel("AH3D - PHaRLAP (km)")
    axes[1, 0].set_ylabel("AH3D - PHaRLAP (km)")
    figure.suptitle("A4-C1Q-1 link-local trajectory residuals")
    figure.tight_layout()
    figure.savefig(residual_figure, dpi=180)
    plt.close(figure)
    budget_figure = stage_dir / "c1q_1_x_error_budget.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, row in zip(axes.flat, rows):
        case_id = row["case_id"]
        s = arrays[f"{case_id}_common_s_km"]
        dx = arrays[f"{case_id}_common_dx_km"]
        axis.plot(s / s[-1], dx - dx[0], color="#1368aa", label="common-path transport")
        axis.axhline(row["landing_dx_km"], color="#c43b3b", linestyle="--",
                     label="landing dx")
        axis.set_title(f"{row['branch_id']} | {row['anchor_solver']} anchor")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[1, 0].set_xlabel("Normalized common geometric path")
    axes[1, 1].set_xlabel("Normalized common geometric path")
    axes[0, 0].set_ylabel("Cumulative dx (km)")
    axes[1, 0].set_ylabel("Cumulative dx (km)")
    figure.suptitle("A4-C1Q-1 cumulative along-link error transport")
    figure.tight_layout()
    figure.savefig(budget_figure, dpi=180)
    plt.close(figure)
    checks = {
        "four_cases": len(rows) == 4,
        "finite_metrics": all(np.isfinite([
            value for key, value in row.items()
            if key not in {"case_id", "branch_id", "anchor_solver"}
        ]).all() for row in rows),
        "budget_closure": all(
            abs(row["budget_closure_km"]) <= max(
                0.01, 0.01 * abs(row["landing_dx_km"])) for row in rows),
        "high_error_dominantly_x": all(
            row["normalized_x_rms_km"] > row["normalized_z_rms_km"]
            and row["normalized_x_rms_km"] > 1000 * row["normalized_y_rms_km"]
            for row in rows if row["branch_id"] == "high_ray"),
        "visualizations": all(_image(path)["passed"]
                              for path in (residual_figure, budget_figure)),
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
    }
    status = ("completed_c1q1_along_track_transport_closure"
              if all(checks.values()) else "blocked_c1q1_transport_closure")
    report = stage_dir / "c1q_1_transport.json"
    _write_json(report, {"schema_version": "a4-c1q-1-transport-v1",
                         "status": status, "checks": checks, "metrics": rows})
    sources = {**previous["source_sha256"],
               "c1q0_manifest": _sha256(
                   OUT / "c1q_0" / "c1q_0_freeze_manifest.json")}
    _manifest(1, status, sources,
              [csv_path, npz_path, report, residual_figure, budget_figure], checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def run_c1q2() -> str:
    previous = _require(1, {"completed_c1q1_along_track_transport_closure"})
    stage_dir = OUT / "c1q_2"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    ph_cases = _ph_cases(context["radius"])
    f_hz = context["frequency_MHz"] * 1e6
    rows, arrays = [], {}
    for case_id, branch, anchor, angle in CASES:
        ph = ph_cases[case_id]
        ah = _ah_case(case_id, context)
        p2 = np.einsum("ij,ij->i", ph["phase_index_vector"],
                       ph["phase_index_vector"])
        n2_index = ph["refractive_index"]**2
        n2_path_ne = 1.0 - cfg.K_FP**2 * ph["electron_density_m3"] / f_hz**2
        provider_ne = np.array([
            context["provider"].get_Ne_ecef(point, context["radius"])
            for point in ph["position_ecef_km"]
        ])
        n2_provider = 1.0 - cfg.K_FP**2 * provider_ne / f_hz**2
        mask = ph["electron_density_m3"] > 0.01 * np.max(
            ph["electron_density_m3"])
        valid_k = mask & (ph["electron_density_m3"] > 0.0) & (p2 < 1.0)
        k_eff = f_hz * np.sqrt(np.maximum(
            1.0 - p2[valid_k], 0.0) / ph["electron_density_m3"][valid_k])
        position_closure = np.linalg.norm(
            ph["position_ecef_km"] - ph["path_position_ecef_km"], axis=1)
        launch_delta_p = np.linalg.norm(
            ah["phase_index_vector"][0] - ph["phase_index_vector"][0])
        row = {
            "case_id": case_id, "branch_id": branch, "anchor_solver": anchor,
            "elevation_deg": angle,
            "state_vs_path_index_rms": _rms(p2 - n2_index),
            "state_vs_path_index_max": _max(p2 - n2_index),
            "path_index_vs_path_ne_rms": _rms(n2_index - n2_path_ne),
            "path_ne_vs_provider_relative_rms": _rms(
                (ph["electron_density_m3"][mask] - provider_ne[mask])
                / np.maximum(provider_ne[mask], 1.0)),
            "state_vs_provider_n2_rms": _rms(p2 - n2_provider),
            "effective_K_FP_median": float(np.median(k_eff)),
            "effective_K_FP_relative_to_project": float(
                np.median(k_eff) / cfg.K_FP - 1.0),
            "effective_K_FP_relative_spread": float(
                np.std(k_eff) / max(np.mean(k_eff), 1e-30)),
            "state_path_position_rms_km": _rms(position_closure),
            "state_path_position_max_km": _max(position_closure),
            "launch_p_vector_difference": float(launch_delta_p),
            "launch_n_ah3d": float(np.linalg.norm(ah["phase_index_vector"][0])),
            "launch_n_pharlap": float(np.linalg.norm(ph["phase_index_vector"][0])),
        }
        rows.append(row)
        u = ph["geom_path_km"] / ph["geom_path_km"][-1]
        arrays.update({
            f"{case_id}_u": u,
            f"{case_id}_p2_minus_index2": p2 - n2_index,
            f"{case_id}_index2_minus_path_ne_n2": n2_index - n2_path_ne,
            f"{case_id}_p2_minus_provider_n2": p2 - n2_provider,
            f"{case_id}_path_ne_relative_provider": (
                ph["electron_density_m3"] - provider_ne)
                / np.maximum(provider_ne, 1.0),
        })
    state_consistent = all(
        row["state_vs_path_index_rms"] <= 1e-10 for row in rows)
    field_difference = all(
        row["path_ne_vs_provider_relative_rms"] > 1e-6 for row in rows)
    constant_stable = all(
        row["effective_K_FP_relative_spread"] <= 1e-4 for row in rows)
    if state_consistent and field_difference:
        attribution = "pharlap_internal_field_or_output_semantics_candidate"
    elif state_consistent and constant_stable and any(
            abs(row["effective_K_FP_relative_to_project"]) > 1e-6
            for row in rows):
        attribution = "physical_constant_candidate"
    else:
        attribution = "state_medium_constant_unresolved"
    csv_path = stage_dir / "c1q_2_state_medium_metrics.csv"
    npz_path = stage_dir / "c1q_2_state_medium_samples.npz"
    _write_csv(csv_path, rows)
    np.savez_compressed(npz_path, **arrays)
    figure_path = stage_dir / "c1q_2_state_medium_closure.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, row in zip(axes.flat, rows):
        case_id = row["case_id"]
        u = arrays[f"{case_id}_u"]
        axis.semilogy(u, np.maximum(np.abs(
            arrays[f"{case_id}_p2_minus_index2"]), 1e-18),
            label="|p2-index2|")
        axis.semilogy(u, np.maximum(np.abs(
            arrays[f"{case_id}_index2_minus_path_ne_n2"]), 1e-18),
            label="|index2-n2(path Ne)|")
        axis.semilogy(u, np.maximum(np.abs(
            arrays[f"{case_id}_p2_minus_provider_n2"]), 1e-18),
            label="|p2-n2(provider)|")
        axis.set_title(f"{row['branch_id']} | {row['anchor_solver']} anchor")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[1, 0].set_xlabel("Normalized geometric path u")
    axes[1, 1].set_xlabel("Normalized geometric path u")
    axes[0, 0].set_ylabel("Absolute closure residual")
    axes[1, 0].set_ylabel("Absolute closure residual")
    figure.suptitle("A4-C1Q-2 PHaRLAP state, medium and constant closure")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    checks = {
        "four_cases": len(rows) == 4,
        "finite_metrics": all(np.isfinite([
            value for key, value in row.items()
            if key not in {"case_id", "branch_id", "anchor_solver"}
        ]).all() for row in rows),
        "path_state_same_length_and_position": all(
            row["state_path_position_max_km"] <= 1e-4 for row in rows),
        "visualization": _image(figure_path)["passed"],
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
    }
    status = ("completed_c1q2_state_medium_constant_audit"
              if all(checks.values()) else "blocked_c1q2_state_medium_audit")
    report = stage_dir / "c1q_2_state_medium.json"
    _write_json(report, {
        "schema_version": "a4-c1q-2-state-medium-v1",
        "status": status, "checks": checks, "attribution": attribution,
        "metrics": rows,
    })
    sources = {**previous["source_sha256"],
               "c1q1_manifest": _sha256(
                   OUT / "c1q_1" / "c1q_1_freeze_manifest.json")}
    _manifest(2, status, sources, [csv_path, npz_path, report, figure_path], checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def _branch_derivative(height: np.ndarray, values: np.ndarray,
                       window: int = 9) -> np.ndarray:
    result = np.full_like(values, np.nan, dtype=float)
    apex = int(np.argmax(height))
    for indices in (np.arange(0, apex + 1), np.arange(apex, len(height))):
        local = _local_poly_derivative(height[indices], values[indices], window)
        result[indices] = local
    return result


def _provider_gradient_n2(position: np.ndarray, context: dict) -> np.ndarray:
    f_hz = context["frequency_MHz"] * 1e6
    gradients = np.array([
        context["provider"].get_grad_Ne_ecef(
            point, context["radius"], step_km=0.25)
        for point in position
    ])
    return -cfg.K_FP**2 * gradients / f_hz**2


def _canonical_metrics(case: dict, context: dict, solver: str) -> tuple[dict, dict]:
    position = case["position_ecef_km"]
    p = case["phase_index_vector"]
    group = case["group_path_km"]
    group = group - group[0]
    gradient_provider = _provider_gradient_n2(position, context)
    predicted_r = position[0] + cumulative_trapezoid(p, group, axis=0, initial=0.0)
    r_closure = np.linalg.norm(position - predicted_r, axis=1)
    predicted_p = p[0] + 0.5 * cumulative_trapezoid(
        gradient_provider, group, axis=0, initial=0.0)
    p_provider_closure = np.linalg.norm(p - predicted_p, axis=1)
    window_metrics = {}
    dr_windows = {}
    dp_windows = {}
    for window in (5, 7, 9):
        dr = _local_poly_derivative(group, position, window)
        dp = _local_poly_derivative(group, p, window)
        dr_windows[window] = dr
        dp_windows[window] = dp
        window_metrics[f"dr_window{window}_rms"] = _rms(
            np.linalg.norm(dr - p, axis=1))
        window_metrics[f"dp_provider_window{window}_rms"] = _rms(
            np.linalg.norm(dp - 0.5 * gradient_provider, axis=1))
    detail = {
        "u": group / group[-1],
        "r_integral_closure_km": r_closure,
        "p_provider_integral_closure": p_provider_closure,
        "dr_window_spread": np.fmax(
            np.linalg.norm(dr_windows[5] - dr_windows[9], axis=1),
            np.linalg.norm(dr_windows[7] - dr_windows[9], axis=1)),
        "dp_window_spread": np.fmax(
            np.linalg.norm(dp_windows[5] - dp_windows[9], axis=1),
            np.linalg.norm(dp_windows[7] - dp_windows[9], axis=1)),
    }
    metric = {
        "solver": solver,
        "r_integral_closure_rms_km": _rms(r_closure),
        "r_integral_closure_max_km": _max(r_closure),
        "p_provider_integral_closure_rms": _rms(p_provider_closure),
        "p_provider_integral_closure_max": _max(p_provider_closure),
        "dr_window_spread_rms": _rms(detail["dr_window_spread"]),
        "dp_window_spread_rms": _rms(detail["dp_window_spread"]),
        **window_metrics,
    }
    if solver == "PHaRLAP":
        n2_effective = case["refractive_index"]**2
        dn2_dh = _branch_derivative(case["height_km"], n2_effective)
        radial = position / np.linalg.norm(position, axis=1)[:, None]
        gradient_effective = dn2_dh[:, None] * radial
        predicted_effective = p[0] + 0.5 * cumulative_trapezoid(
            np.nan_to_num(gradient_effective), group, axis=0, initial=0.0)
        p_effective_closure = np.linalg.norm(p - predicted_effective, axis=1)
        dp9 = dp_windows[9]
        metric.update({
            "p_effective_integral_closure_rms": _rms(p_effective_closure),
            "p_effective_integral_closure_max": _max(p_effective_closure),
            "dp_effective_window9_rms": _rms(np.linalg.norm(
                dp9 - 0.5 * gradient_effective, axis=1)),
        })
        detail.update({
            "p_effective_integral_closure": p_effective_closure,
            "gradient_effective": gradient_effective,
        })
    return metric, detail


def run_c1q3() -> str:
    previous = _require(2, {"completed_c1q2_state_medium_constant_audit"})
    stage_dir = OUT / "c1q_3"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    ph_cases = _ph_cases(context["radius"])
    rows, arrays = [], {}
    for case_id, branch, anchor, angle in CASES:
        for solver, case in (
                ("PHaRLAP", ph_cases[case_id]), ("AH3D", _ah_case(case_id, context))):
            metric, detail = _canonical_metrics(case, context, solver)
            row = {
                "case_id": case_id, "branch_id": branch,
                "anchor_solver": anchor, "elevation_deg": angle, **metric,
            }
            if solver != "PHaRLAP":
                row.update({
                    "p_effective_integral_closure_rms": np.nan,
                    "p_effective_integral_closure_max": np.nan,
                    "dp_effective_window9_rms": np.nan,
                })
            rows.append(row)
            arrays.update({
                f"{solver}_{case_id}_{key}": value
                for key, value in detail.items()
            })
    csv_path = stage_dir / "c1q_3_canonical_metrics.csv"
    npz_path = stage_dir / "c1q_3_canonical_samples.npz"
    _write_csv(csv_path, rows)
    np.savez_compressed(npz_path, **arrays)
    figure_path = stage_dir / "c1q_3_canonical_integral_closure.png"
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for axis, (case_id, branch, anchor, _) in zip(axes.flat, CASES):
        for solver, color in (("PHaRLAP", "#1368aa"), ("AH3D", "#c43b3b")):
            u = arrays[f"{solver}_{case_id}_u"]
            axis.semilogy(u, np.maximum(
                arrays[f"{solver}_{case_id}_r_integral_closure_km"], 1e-15),
                color=color, label=f"{solver} r integral")
            axis.semilogy(u, np.maximum(
                arrays[f"{solver}_{case_id}_p_provider_integral_closure"], 1e-15),
                color=color, linestyle="--", label=f"{solver} p/provider")
        axis.semilogy(
            arrays[f"PHaRLAP_{case_id}_u"],
            np.maximum(arrays[
                f"PHaRLAP_{case_id}_p_effective_integral_closure"], 1e-15),
            color="#16836b", linestyle=":", label="PHaRLAP p/effective n2")
        axis.set_title(f"{branch} | {anchor} anchor")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    axes[1, 0].set_xlabel("Normalized group path")
    axes[1, 1].set_xlabel("Normalized group path")
    axes[0, 0].set_ylabel("Integral closure")
    axes[1, 0].set_ylabel("Integral closure")
    figure.suptitle("A4-C1Q-3 B=0 canonical-equation integral closure")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    ph_rows = [row for row in rows if row["solver"] == "PHaRLAP"]
    effective_better = all(
        row["p_effective_integral_closure_rms"]
        < row["p_provider_integral_closure_rms"] for row in ph_rows)
    checks = {
        "eight_solver_case_rows": len(rows) == 8,
        "finite_required_metrics": all(np.isfinite([
            row["r_integral_closure_rms_km"],
            row["r_integral_closure_max_km"],
            row["p_provider_integral_closure_rms"],
            row["p_provider_integral_closure_max"],
            row["dr_window_spread_rms"],
            row["dp_window_spread_rms"],
        ]).all() for row in rows),
        "three_derivative_windows": all(
            all(f"dr_window{window}_rms" in row for window in (5, 7, 9))
            for row in rows),
        "visualization": _image(figure_path)["passed"],
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
    }
    status = ("completed_c1q3_b0_canonical_isolation"
              if all(checks.values()) else "blocked_c1q3_canonical_isolation")
    report = stage_dir / "c1q_3_canonical.json"
    _write_json(report, {
        "schema_version": "a4-c1q-3-canonical-v1",
        "status": status, "checks": checks,
        "effective_pharlap_n2_closes_momentum_better_than_common_provider":
            effective_better,
        "metrics": rows,
    })
    sources = {**previous["source_sha256"],
               "c1q2_manifest": _sha256(
                   OUT / "c1q_2" / "c1q_2_freeze_manifest.json")}
    _manifest(3, status, sources, [csv_path, npz_path, report, figure_path], checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


class _PharlapEffectiveProvider(GeocentricNeProvider3D):
    name = "c1q_validation_only_pharlap_effective_profile"

    def __init__(self, radius: float, base: GeocentricNeProvider3D,
                 ph_cases: dict):
        self.earth_radius_km = float(radius)
        self.base = base
        heights = np.concatenate([case["height_km"] for case in ph_cases.values()])
        densities = np.concatenate([
            case["electron_density_m3"] for case in ph_cases.values()])
        order = np.argsort(heights)
        heights = heights[order]
        densities = densities[order]
        keys = np.round(heights, 4)
        unique = np.unique(keys)
        profile_height = np.array([np.mean(heights[keys == key]) for key in unique])
        profile_density = np.array([np.mean(densities[keys == key]) for key in unique])
        self.lower = float(profile_height[0])
        self.upper = float(profile_height[-1])
        self.profile = PchipInterpolator(
            profile_height, np.maximum(profile_density, 0.0), extrapolate=False)

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        if self.lower <= height_km <= self.upper:
            return max(float(self.profile(height_km)), 0.0)
        return self.base.get_Ne_geocentric(lat_deg, lon_deg, height_km, t=t)

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                 height_km: float) -> float:
        return self.base.domain_margin_geocentric(lat_deg, lon_deg, height_km)

    def get_metadata(self) -> dict:
        frozen = {
            "name": self.name,
            "earth_radius_km": self.earth_radius_km,
            "coordinate_system": "geocentric_spherical",
            "profile_range_km": [self.lower, self.upper],
            "semantics": "validation_only_direct_PHARLAP_output_profile_no_production_use",
        }
        return {**super().get_metadata(), **frozen,
                "provider_hash": hashlib.sha256(json.dumps(
                    frozen, sort_keys=True).encode("ascii")).hexdigest()}


@contextmanager
def _patched(**replacements):
    originals = {name: getattr(ahcore, name) for name in replacements}
    try:
        for name, value in replacements.items():
            setattr(ahcore, name, value)
        yield
    finally:
        for name, value in originals.items():
            setattr(ahcore, name, value)


def _exact_b0_rhs(path_km, state, model, *, medium_time_s=0.0):
    del path_km
    ne = model.provider.get_Ne_ecef(
        state[:3], model.earth_radius_km, t=medium_time_s)
    f_hz = model.freq_MHz * 1e6
    n2 = 1.0 - cfg.K_FP**2 * ne / f_hz**2
    if n2 <= 0.0:
        raise ValueError("exact B0 diagnostic entered evanescence")
    n = np.sqrt(n2)
    grad_ne = model.provider.get_grad_Ne_ecef(
        state[:3], model.earth_radius_km, t=medium_time_s, step_km=0.25)
    grad_n2 = -cfg.K_FP**2 * grad_ne / f_hz**2
    direction = state[3:6] / np.linalg.norm(state[3:6])
    return np.r_[direction, model.k0 * grad_n2 / (2.0 * n),
                 1.0 / (C_KM_S * n)]


def _trace_current(context: dict, elevation: float, provider=None,
                   variant: str = "current", ph: dict | None = None) -> tuple[dict, dict]:
    local = {**context, "provider": context["provider"] if provider is None else provider}
    replacements = {}
    if variant == "exact_b0_rhs":
        replacements["ah3d_ecef_rhs"] = _exact_b0_rhs
    elif variant == "no_projection":
        replacements["project_dispersion_ah3d_ecef"] = (
            lambda state, model, **kwargs: (
                np.asarray(state, dtype=float).copy(), 0.0, True))
    elif variant == "pharlap_launch":
        if ph is None:
            raise ValueError("PHaRLAP launch variant requires a PHaRLAP state")
        original_launch = ahcore.launch_state_ah3d_ecef

        def ph_launch(model, **kwargs):
            state = original_launch(model, **kwargs)
            state[:3] = ph["position_ecef_km"][0]
            state[3:6] = model.k0 * ph["phase_index_vector"][0]
            state[6] = 0.0
            return state
        replacements["launch_state_ah3d_ecef"] = ph_launch
    with _patched(**replacements):
        return prep3._trace(prep3._model(local), local, elevation)


def _trace_group_parameterized(context: dict, elevation: float) -> tuple[dict, dict]:
    model = prep3._model(context)
    settings = cfg.UPDATE2["AH3D_INTEGRATION"][prep3.PROFILE]
    initial = ahcore.launch_state_ah3d_ecef(
        model, tx_lat_deg=context["tx"][0], tx_lon_deg=context["tx"][1],
        tx_height_km=context["tx"][2], elevation_deg=elevation,
        bearing_deg=context["bearing_deg"])
    initial_state = np.r_[initial[:3], initial[3:6] / model.k0, 0.0]
    f_hz = context["frequency_MHz"] * 1e6

    def rhs(group_path, state):
        ne = context["provider"].get_Ne_ecef(state[:3], context["radius"])
        n2 = 1.0 - cfg.K_FP**2 * ne / f_hz**2
        grad_ne = context["provider"].get_grad_Ne_ecef(
            state[:3], context["radius"], step_km=0.25)
        grad_n2 = -cfg.K_FP**2 * grad_ne / f_hz**2
        return np.r_[state[3:6], 0.5 * grad_n2, np.sqrt(max(n2, 1e-30))]

    def ground(group_path, state):
        if group_path < 1e-4:
            return 1.0
        return np.linalg.norm(state[:3]) - context["radius"]
    ground.direction = -1
    ground.terminal = True
    solution = solve_ivp(
        rhs, (0.0, 4000.0), initial_state, method="DOP853",
        rtol=settings["rtol"],
        atol=np.r_[
            np.full(3, settings["atol_r_km"]),
            np.full(3, settings["atol_k"] / model.k0),
            settings["atol_r_km"],
        ],
        max_step=settings["ds_max_km"], events=ground)
    if not solution.success or not solution.t_events[0].size:
        raise RuntimeError("group-parameterized C1Q shadow trace did not reach ground")
    position = solution.y[:3].T
    p = solution.y[3:6].T
    geom = solution.y[6]
    geocentric = ecef_to_geocentric(position, context["radius"])
    ray = {
        "trajectory_ecef_km": position,
        "trajectory_geocentric_lat_lon_height": geocentric,
        "trajectory_physical_height_km": geocentric[:, 2],
        "trajectory_geom_path_km": geom,
        "wavevector_ecef_km^-1": model.k0 * p,
        "group_path_km": float(solution.t[-1]),
        "geom_path_km": float(geom[-1]),
        "apex_height_km": float(np.max(geocentric[:, 2])),
        "termination_status": "ground",
        "hamiltonian_residual_max": float(np.max(np.abs(
            np.einsum("ij,ij->i", p, p) - np.array([
                1.0 - cfg.K_FP**2 * context["provider"].get_Ne_ecef(
                    point, context["radius"]) / f_hz**2 for point in position])))),
    }
    record = {
        "ground_range_km": float(_xyz(position, context)[-1, 0]),
        "apex_km": ray["apex_height_km"],
        "geometric_path_km": ray["geom_path_km"],
        "group_path_km": ray["group_path_km"],
        "status": "ground",
    }
    return record, ray


def _variant_path(ray: dict) -> dict:
    return {
        "trajectory_ecef_km": np.asarray(ray["trajectory_ecef_km"]),
        "trajectory_geocentric_lat_lon_height": np.asarray(
            ray["trajectory_geocentric_lat_lon_height"]),
        "trajectory_physical_height_km": np.asarray(
            ray["trajectory_physical_height_km"]),
        "trajectory_geom_path_km": np.asarray(ray["trajectory_geom_path_km"]),
    }


def _variant_metrics(ph: dict, ray: dict, context: dict) -> dict:
    metric, _ = _trajectory_metrics(
        _path_dict(ph), _variant_path(ray), context["radius"])
    u = np.linspace(0.0, 1.0, 1001)
    ph_position = _resample(
        ph["geom_path_km"], ph["position_ecef_km"], u * ph["geom_path_km"][-1])
    ah_position = _resample(
        ray["trajectory_geom_path_km"], ray["trajectory_ecef_km"],
        u * ray["trajectory_geom_path_km"][-1])
    delta = _xyz(ah_position, context) - _xyz(ph_position, context)
    ph_xyz = _xyz(ph["position_ecef_km"], context)
    ah_xyz = _xyz(ray["trajectory_ecef_km"], context)
    return {
        "trajectory_rms_3d_km": metric["trajectory_rms_3d_km"],
        "trajectory_max_3d_km": metric["trajectory_max_3d_km"],
        "x_rms_km": _rms(delta[:, 0]),
        "y_rms_km": _rms(delta[:, 1]),
        "z_rms_km": _rms(delta[:, 2]),
        "landing_dx_km": float(ah_xyz[-1, 0] - ph_xyz[-1, 0]),
        "apex_error_km": float(abs(
            np.max(ray["trajectory_physical_height_km"]) - np.max(ph["height_km"]))),
        "geom_path_relative_error": float(abs(
            ray["trajectory_geom_path_km"][-1] - ph["geom_path_km"][-1])
            / ph["geom_path_km"][-1]),
        "group_path_relative_error": float(abs(
            ray.get("group_path_km", np.nan) - ph["group_path_km"][-1])
            / ph["group_path_km"][-1]),
    }


def run_c1q4() -> str:
    previous = _require(3, {"completed_c1q3_b0_canonical_isolation"})
    stage_dir = OUT / "c1q_4"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    ph_cases = _ph_cases(context["radius"])
    effective_provider = _PharlapEffectiveProvider(
        context["radius"], context["provider"], ph_cases)
    variants = (
        "baseline", "pharlap_launch", "exact_b0_rhs", "group_path_parameterized",
        "no_projection", "pharlap_effective_ne",
    )
    rows = []
    raw_dir = stage_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for case_id, branch, anchor, angle in CASES:
        ph = ph_cases[case_id]
        baseline_case = _ah_case(case_id, context)
        baseline_ray = {
            "trajectory_ecef_km": baseline_case["position_ecef_km"],
            "trajectory_geocentric_lat_lon_height": np.column_stack((
                baseline_case["lat_deg"], baseline_case["lon_deg"],
                baseline_case["height_km"])),
            "trajectory_physical_height_km": baseline_case["height_km"],
            "trajectory_geom_path_km": baseline_case["geom_path_km"],
            "group_path_km": baseline_case["group_path_km"][-1],
        }
        baseline = _variant_metrics(ph, baseline_ray, context)
        for variant in variants:
            if variant == "baseline":
                ray = baseline_ray
                record = {"status": "ground"}
            elif variant == "group_path_parameterized":
                record, ray = _trace_group_parameterized(context, angle)
            elif variant == "pharlap_effective_ne":
                record, ray = _trace_current(
                    context, angle, provider=effective_provider, variant="current")
            else:
                record, ray = _trace_current(
                    context, angle, variant=variant, ph=ph)
            metric = _variant_metrics(ph, ray, context)
            closure = (
                (baseline["x_rms_km"] - metric["x_rms_km"])
                / max(baseline["x_rms_km"], 1e-30))
            row = {
                "case_id": case_id, "branch_id": branch,
                "anchor_solver": anchor, "elevation_deg": angle,
                "variant": variant, **metric,
                "x_error_closure_fraction": float(closure),
                "termination": ray.get("termination_status", record.get("status")),
                "hamiltonian_residual_max": float(
                    ray.get("hamiltonian_residual_max", np.nan)),
            }
            rows.append(row)
            arrays = {
                key: np.asarray(value) for key, value in ray.items()
                if isinstance(value, np.ndarray)
            }
            np.savez_compressed(raw_dir / f"{case_id}_{variant}.npz", **arrays)
    documented_interpolation = {
        "variant": "pharlap_documented_interpolation",
        "status": "not_run_no_public_pharlap_4p7p4_interpolation_kernel",
        "production_eligible": False,
    }
    high = [row for row in rows if row["branch_id"] == "high_ray"]
    low = [row for row in rows if row["branch_id"] == "low_ray"]
    safe_variants = ("pharlap_launch", "exact_b0_rhs",
                     "group_path_parameterized", "no_projection")
    primary = []
    for variant in safe_variants:
        high_rows = [row for row in high if row["variant"] == variant]
        low_rows = [row for row in low if row["variant"] == variant]
        if (len(high_rows) == 2 and all(
                row["x_error_closure_fraction"] >= 0.8
                and row["trajectory_rms_3d_km"] <= 0.5 for row in high_rows)
                and all(row["trajectory_rms_3d_km"] <= 0.15 for row in low_rows)):
            primary.append(variant)
    effective_rows = [
        row for row in high if row["variant"] == "pharlap_effective_ne"]
    effective_primary = len(effective_rows) == 2 and all(
        row["x_error_closure_fraction"] >= 0.8 for row in effective_rows)
    if len(primary) == 1:
        attribution = "safe_primary_candidate"
    elif len(primary) > 1:
        attribution = "multiple_safe_candidates"
    elif effective_primary:
        attribution = "pharlap_effective_field_candidate_not_production_safe"
    else:
        attribution = "no_single_factor_primary_candidate"
    csv_path = stage_dir / "c1q_4_shadow_metrics.csv"
    _write_csv(csv_path, rows)
    figure_path = stage_dir / "c1q_4_shadow_closure_matrix.png"
    figure, axes = plt.subplots(1, 2, figsize=(15, 6))
    labels = [variant.replace("_", "\n") for variant in variants]
    x = np.arange(len(variants))
    width = 0.18
    for index, (case_id, branch, anchor, _) in enumerate(CASES):
        case_rows = {row["variant"]: row for row in rows if row["case_id"] == case_id}
        axes[0].bar(x + (index - 1.5) * width,
                    [case_rows[variant]["x_error_closure_fraction"]
                     for variant in variants], width,
                    label=f"{branch}/{anchor}")
        axes[1].bar(x + (index - 1.5) * width,
                    [case_rows[variant]["trajectory_rms_3d_km"]
                     for variant in variants], width,
                    label=f"{branch}/{anchor}")
    axes[0].axhline(0.8, color="black", linestyle="--", label="80% primary gate")
    axes[0].set_ylabel("x error closure fraction")
    axes[1].axhline(0.5, color="black", linestyle="--", label="high RMS gate")
    axes[1].set_ylabel("3D trajectory RMS (km)")
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_title("Single-factor along-track closure")
    axes[1].set_title("Single-factor trajectory RMS")
    figure.suptitle("A4-C1Q-4 single-factor shadow replay")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    checks = {
        "six_variants_four_cases": len(rows) == 24,
        "all_production_eligible_variants_ground": all(
            row["termination"] == "ground" for row in rows
            if row["variant"] != "pharlap_effective_ne"),
        "effective_profile_outcome_explicit": all(
            row["termination"] in {"ground", "physics_domain_rejected"}
            for row in rows if row["variant"] == "pharlap_effective_ne"),
        "finite_required_metrics": all(np.isfinite([
            row["trajectory_rms_3d_km"], row["trajectory_max_3d_km"],
            row["x_rms_km"], row["y_rms_km"], row["z_rms_km"],
            row["landing_dx_km"], row["x_error_closure_fraction"],
        ]).all() for row in rows),
        "no_empirical_production_fix": not documented_interpolation["production_eligible"],
        "visualization": _image(figure_path)["passed"],
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
    }
    status = ("completed_c1q4_single_factor_shadow_replay"
              if all(checks.values()) else "blocked_c1q4_shadow_replay")
    report = stage_dir / "c1q_4_shadow.json"
    _write_json(report, {
        "schema_version": "a4-c1q-4-shadow-v1",
        "status": status, "checks": checks, "attribution": attribution,
        "safe_primary_candidates": primary,
        "pharlap_effective_field_primary": effective_primary,
        "documented_interpolation": documented_interpolation,
        "metrics": rows,
    })
    sources = {**previous["source_sha256"],
               "c1q3_manifest": _sha256(
                   OUT / "c1q_3" / "c1q_3_freeze_manifest.json")}
    artifacts = [csv_path, report, figure_path, *sorted(raw_dir.glob("*.npz"))]
    _manifest(4, status, sources, artifacts, checks)
    if not status.startswith("completed"):
        raise RuntimeError(status)
    return status


def run_c1q5() -> str:
    previous = _require(4, {"completed_c1q4_single_factor_shadow_replay"})
    stage_dir = OUT / "c1q_5"
    stage_dir.mkdir(parents=True, exist_ok=True)
    shadow = json.loads(
        (OUT / "c1q_4" / "c1q_4_shadow.json").read_text(encoding="utf-8"))
    candidates = shadow["safe_primary_candidates"]
    if candidates:
        status = "blocked_c1q5_safe_candidate_requires_explicit_code_repair"
        entered = True
    else:
        status = "not_entered_c1q5_no_safe_primary_candidate"
        entered = False
    checks = {
        "production_repair_gate_applied": True,
        "production_sources_unchanged": all(
            _sha256(path) == previous["source_sha256"][key]
            for path, key in (
                (AH_SOLVER, "ah_solver"), (AH_DISPERSION, "ah_dispersion"),
                (AH_PROVIDER, "ah_provider"), (CONFIG, "config"))),
        "no_empirical_fix": not shadow["documented_interpolation"]["production_eligible"],
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
    }
    figure_path = stage_dir / "c1q_5_production_gate.png"
    figure, axis = plt.subplots(figsize=(10, 4.8))
    axis.bar(["safe primary\ncandidate", "production\nmodified"],
             [int(bool(candidates)), 0], color=["#d28e00", "#187b3f"])
    axis.set_ylim(0, 1.1)
    axis.set_yticks((0, 1), ("NO", "YES"))
    axis.set_title("A4-C1Q-5 production repair gate")
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    checks["visualization"] = _image(figure_path)["passed"]
    report = stage_dir / "c1q_5_production_gate.json"
    _write_json(report, {
        "schema_version": "a4-c1q-5-production-gate-v1",
        "status": status, "checks": checks, "entered_repair": entered,
        "safe_primary_candidates": candidates,
        "repair_acceptance": {
            "high_trajectory_rms_km": 0.5,
            "high_max_and_landing_km": 1.0,
            "low_trajectory_rms_km": 0.15,
            "low_landing_km": 0.2,
            "group_and_geom_relative": 0.001,
        },
    })
    sources = {**previous["source_sha256"],
               "c1q4_manifest": _sha256(
                   OUT / "c1q_4" / "c1q_4_freeze_manifest.json")}
    _manifest(5, status, sources, [report, figure_path], checks)
    if not all(checks.values()):
        raise RuntimeError("C1Q-5 production gate integrity failed")
    return status


def run_c1q6() -> str:
    previous = _require(5, {
        "not_entered_c1q5_no_safe_primary_candidate",
        "blocked_c1q5_safe_candidate_requires_explicit_code_repair",
    })
    stage_dir = OUT / "c1q_6"
    stage_dir.mkdir(parents=True, exist_ok=True)
    medium = json.loads(
        (OUT / "c1q_2" / "c1q_2_state_medium.json").read_text(encoding="utf-8"))
    canonical = json.loads(
        (OUT / "c1q_3" / "c1q_3_canonical.json").read_text(encoding="utf-8"))
    shadow = json.loads(
        (OUT / "c1q_4" / "c1q_4_shadow.json").read_text(encoding="utf-8"))
    gate = json.loads(
        (OUT / "c1q_5" / "c1q_5_production_gate.json").read_text(encoding="utf-8"))
    if gate["safe_primary_candidates"]:
        final = "identified_but_not_safely_repairable"
    elif shadow["pharlap_effective_field_primary"]:
        final = "identified_but_not_safely_repairable"
    elif (medium["attribution"]
          == "pharlap_internal_field_or_output_semantics_candidate"
          and canonical[
              "effective_pharlap_n2_closes_momentum_better_than_common_provider"]):
        final = "multiple_unresolved_contributors"
    else:
        final = "pharlap_effective_dynamics_not_recoverable"
    checks = {
        "final_state_allowed": final in FINAL_STATES,
        "no_automatic_c2_admission": True,
        "formal_c1_failure_preserved": json.loads(
            (C1 / "final" / "a4_c1_validation.json").read_text(
                encoding="utf-8"))["status"] == FINAL_C1_STATUS,
        "production_sources_unchanged_without_accepted_repair": all(
            _sha256(path) == previous["source_sha256"][key]
            for path, key in (
                (AH_SOLVER, "ah_solver"), (AH_DISPERSION, "ah_dispersion"),
                (AH_PROVIDER, "ah_provider"), (CONFIG, "config"))),
    }
    figure_path = stage_dir / "c1q_6_final_attribution.png"
    labels = list(FINAL_STATES)
    values = [int(label == final) for label in labels]
    figure, axis = plt.subplots(figsize=(12, 5.5))
    axis.barh(np.arange(len(labels)), values,
              color=["#b22222" if value else "#c8c8c8" for value in values])
    axis.set_yticks(np.arange(len(labels)), [label.replace("_", " ") for label in labels])
    axis.set_xticks((0, 1), ("NO", "SELECTED"))
    axis.set_xlim(0, 1.05)
    axis.set_title("A4-C1Q-6 final attribution and manual decision gate")
    figure.tight_layout()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    checks["visualization"] = _image(figure_path)["passed"]
    status = f"completed_c1q6_{final}"
    report = stage_dir / "c1q_6_final_attribution.json"
    _write_json(report, {
        "schema_version": "a4-c1q-6-final-v1",
        "status": status, "final_attribution": final, "checks": checks,
        "state_medium_attribution": medium["attribution"],
        "canonical_effective_field_evidence": canonical[
            "effective_pharlap_n2_closes_momentum_better_than_common_provider"],
        "shadow_attribution": shadow["attribution"],
        "safe_primary_candidates": shadow["safe_primary_candidates"],
        "next_gate": (
            "manual_decision_required_before_A4_C2; "
            "recommend_third_auditable_Hamiltonian_solver_if_no_safe_repair"),
    })
    sources = {**previous["source_sha256"],
               "c1q5_manifest": _sha256(
                   OUT / "c1q_5" / "c1q_5_freeze_manifest.json")}
    _manifest(6, status, sources, [report, figure_path], checks)
    if not all(checks.values()):
        raise RuntimeError("C1Q-6 final gate integrity failed")
    return status


STAGES = {
    0: run_c1q0, 1: run_c1q1, 2: run_c1q2, 3: run_c1q3,
    4: run_c1q4, 5: run_c1q5, 6: run_c1q6,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=int, choices=STAGES)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and args.stage is None:
        parser.error("use --stage N or --all")
    selected = tuple(STAGES) if args.all else (args.stage,)
    statuses = {}
    for stage in selected:
        statuses[str(stage)] = STAGES[stage]()
        print(f"C1Q-{stage}: {statuses[str(stage)]}")
    print(json.dumps(statuses, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
