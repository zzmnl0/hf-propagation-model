"""Run A4-C1-0 admission and C1-1 fixed-same-angle comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_prep3 import _load_context, _model, _trace
from tests._run_a4_c1_prep4 import _turning_count
from utils import ecef_to_geocentric, enu_basis, geocentric_to_ecef


PREP = ROOT / "output2" / "a4_pharlap" / "c1_prep"
PREP1_REPORT = PREP / "prep1" / "a4_c1_prep1_validation.json"
PREP5_REPORT = PREP / "prep5" / "a4_c1_prep5_validation.json"
PREP5_MANIFEST = PREP / "prep5" / "a4_c1_prep5_freeze_manifest.json"
PREP6_REPORT = PREP / "prep6" / "a4_c1_prep6_validation.json"
PREP6_MANIFEST = PREP / "prep6" / "a4_c1_prep6_freeze_manifest.json"

BASE = ROOT / "output2" / "a4_pharlap" / "c1"
ADMISSION = BASE / "admission"
SAME = BASE / "same_angle"
RAW = SAME / "raw"
ADMISSION_REPORT = ADMISSION / "c1_0_admission.json"
ADMISSION_FIGURE = ADMISSION / "c1_0_source_status.png"
ADMISSION_MANIFEST = ADMISSION / "c1_0_freeze_manifest.json"
PH_SUMMARY = RAW / "c1_same_angle_pharlap_summary.json"
PH_LOG = RAW / "c1_same_angle_pharlap.log"
METRICS_CSV = SAME / "c1_same_angle_metrics.csv"
METRICS_JSON = SAME / "c1_same_angle_metrics.json"
TRAJECTORY_FIGURE = SAME / "c1_same_angle_trajectories.png"
RESIDUAL_FIGURE = SAME / "c1_same_angle_residuals.png"
METRICS_FIGURE = SAME / "c1_same_angle_metrics.png"
SAME_REPORT = SAME / "c1_1_validation.json"
SAME_MANIFEST = SAME / "c1_1_freeze_manifest.json"

MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SOURCE = ROOT / "tests" / "matlab" / "a4_c1_same_angle_pharlap.m"
MATLAB_SCRIPT = MATLAB_DIR / MATLAB_SOURCE.name
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
AH_SOLVER = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
AH_PROVIDER = ROOT / "models" / "ne_provider_geocentric.py"
AH_DRIVER = ROOT / "tests" / "_run_a4_c1_prep3.py"
CONFIG = ROOT / "config.py"

EXPECTED_INPUT_HASH = "0a978478b52f5363c93cbc77b0ee1868a4f50e488c82c829eb7794f0cf1235a7"
CASES = (
    ("low_pharlap_anchor", "low_ray", "PHaRLAP", 21.94290326504532),
    ("low_ah3d_anchor", "low_ray", "AH3D", 21.94562872092976),
    ("high_pharlap_anchor", "high_ray", "PHaRLAP", 39.54254824570964),
    ("high_ah3d_anchor", "high_ray", "AH3D", 39.58701138438118),
)
LIMITS = {
    "landing_absolute_error_km": 5.0,
    "landing_relative_error": 0.01,
    "apex_height_absolute_error_km": 2.0,
    "geometric_path_relative_error": 0.01,
    "group_path_relative_error": 0.01,
    "group_delay_relative_error": 0.01,
    "trajectory_rms_3d_km": 2.0,
    "trajectory_max_3d_km": 5.0,
}
REPEAT_SCALARS = (
    "ground_range_km", "along_track_residual_km", "cross_track_km", "apex_km",
    "geometric_path_km", "group_path_km", "tau_ms", "landing_lat_deg",
    "landing_lon_deg",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _verify_manifest(path: Path, expected_status: str) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != expected_status:
        raise RuntimeError(f"unexpected frozen status: {path}")
    for item in manifest.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"frozen artifact changed: {artifact}")
    return manifest


def _image_check(path: Path) -> dict:
    pixels = mpimg.imread(path)
    result = {
        "path": str(path), "sha256": _sha256(path),
        "width_px": int(pixels.shape[1]), "height_px": int(pixels.shape[0]),
        "finite_pixels": bool(np.isfinite(pixels).all()),
        "pixel_std": float(np.std(pixels)),
    }
    result["passed"] = bool(result["width_px"] >= 1200
                            and result["height_px"] >= 700
                            and result["finite_pixels"]
                            and result["pixel_std"] > 1.0e-4)
    return result


def _refined_input() -> tuple[Path, dict]:
    report = json.loads(PREP1_REPORT.read_text(encoding="utf-8"))
    record = next(item for item in report["inputs"] if item["variant"] == "refined_1km")
    return Path(record["path"]), record


def _current_sources(input_path: Path) -> dict:
    return {
        "runner": _sha256(Path(__file__)),
        "refined_input": _sha256(input_path),
        "matlab_source": _sha256(MATLAB_SOURCE),
        "matlab_deployed": _sha256(MATLAB_SCRIPT),
        "pharlap_mex": _sha256(MEX),
        "ah_solver": _sha256(AH_SOLVER),
        "ah_provider": _sha256(AH_PROVIDER),
        "ah_prep3_driver": _sha256(AH_DRIVER),
        "config": _sha256(CONFIG),
    }


def _plot_admission(checks: dict, source_tag: str) -> dict:
    labels = list(checks)
    values = [1 if checks[label] else 0 for label in labels]
    figure, axis = plt.subplots(figsize=(11, 5.5))
    colors = ["#187b3f" if value else "#b22222" for value in values]
    axis.barh(np.arange(len(labels)), values, color=colors)
    axis.set_yticks(np.arange(len(labels)), [label.replace("_", " ") for label in labels])
    axis.set_xlim(0, 1.05)
    axis.set_xticks((0, 1), ("FAIL", "PASS"))
    axis.grid(axis="x", alpha=0.2)
    axis.set_title(
        "A4-C1-0 frozen admission | Chapman F2 | one-way | 10 MHz | B=0 | collisions=off\n"
        f"hash/status evidence only | source={source_tag} | scattering/power=N/A")
    figure.tight_layout()
    figure.savefig(ADMISSION_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(ADMISSION_FIGURE)


def run_admission() -> int:
    ADMISSION.mkdir(parents=True, exist_ok=True)
    try:
        commands = []
        for script in ("_run_a4_c1_prep5.py", "_run_a4_c1_prep6.py"):
            result = subprocess.run(
                [sys.executable, str(ROOT / "tests" / script)], cwd=ROOT,
                text=True, capture_output=True, timeout=300, errors="replace")
            commands.append({"script": script, "returncode": result.returncode,
                             "stdout": result.stdout.strip(), "stderr": result.stderr.strip()})
        prep5_manifest = _verify_manifest(
            PREP5_MANIFEST, "passed_1km_grid_stability_and_signatures")
        prep6_manifest = _verify_manifest(PREP6_MANIFEST, "passed_c1_prep_admission")
        prep5 = json.loads(PREP5_REPORT.read_text(encoding="utf-8"))
        prep6 = json.loads(PREP6_REPORT.read_text(encoding="utf-8"))
        ph_environment = json.loads((PREP / "prep5" / "checkpoints"
                                     / "pharlap_refined_1km_summary.json").read_text(
                                         encoding="utf-8"))
        input_path, refined = _refined_input()
        sources = _current_sources(input_path)
        frozen_source = prep5["source_sha256"]
        signatures = prep5["frozen_signatures"]
        frozen_angles = {(item["solver"], item["branch_id"]): item["elevation_deg"]
                         for item in signatures}
        expected_angles = {
            ("PHaRLAP", "low_ray"): CASES[0][3],
            ("AH3D", "low_ray"): CASES[1][3],
            ("PHaRLAP", "high_ray"): CASES[2][3],
            ("AH3D", "high_ray"): CASES[3][3],
        }
        checks = {
            "prep5_strict_replay": commands[0]["returncode"] == 0,
            "prep6_strict_replay": commands[1]["returncode"] == 0,
            "prep5_manifest_frozen": prep5_manifest["status"].startswith("passed"),
            "prep6_manifest_frozen": prep6_manifest["status"].startswith("passed"),
            "refined_1km_input_hash": sources["refined_input"] == EXPECTED_INPUT_HASH
                == refined["sha256"],
            "four_frozen_roots_exact": len(signatures) == 4 and all(
                abs(float(frozen_angles[key]) - value) <= 1.0e-12
                for key, value in expected_angles.items()),
            "four_ground_single_turning_signatures": all(
                item["termination"] == "ground" and item["turning_count"] == 1
                and len(item["trajectory_ecef_sha256"]) == 64 for item in signatures),
            "ah3d_source_hashes_frozen": all(
                sources[key] == frozen_source[key]
                for key in ("ah_solver", "ah_provider", "ah_prep3_driver", "config")),
            "pharlap_mex_hash_frozen": sources["pharlap_mex"] == frozen_source["pharlap_mex"],
            "c1_matlab_source_deployed_exact": sources["matlab_source"]
                == sources["matlab_deployed"],
            "matlab_r2024a_available": MATLAB.is_file()
                and ph_environment.get("matlab_release") == "2024a",
            "pharlap_4p7p4_mex_selected": PHARLAP.is_dir() and MEX.is_file()
                and ph_environment.get("pharlap_version") == "4.7.4"
                and Path(ph_environment["raytrace_3d_sp_path"]).resolve() == MEX.resolve(),
            "prep6_admission_status": prep6["status"] == "passed_c1_prep_admission",
        }
        source_tag = hashlib.sha256(json.dumps(
            sources, sort_keys=True).encode("ascii")).hexdigest()[:12]
        figure = _plot_admission(checks, source_tag)
        checks["visualization"] = figure["passed"]
        status = "passed_c1_0_frozen_admission" if all(checks.values()) \
            else "blocked_input_or_source_drift"
        report = {
            "schema_version": "a4-c1-0-admission-v1", "status": status,
            "case_semantics": "frozen_source_and_environment_admission_not_cross_solver_validation",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "checks": checks, "strict_replay": commands,
            "source_sha256": sources, "source_set_short_hash": source_tag,
            "frozen_angles_deg": [dict(case_id=item[0], branch_id=item[1],
                                        anchor_solver=item[2], elevation_deg=item[3])
                                  for item in CASES],
            "prep5_manifest_sha256": _sha256(PREP5_MANIFEST),
            "prep6_manifest_sha256": _sha256(PREP6_MANIFEST),
            "visualization": figure,
            "limitations": [
                "C1-0 checks frozen provenance only; no cross-solver metric is accepted here",
                "The case is B=0 smooth Chapman F2 and does not validate O/X or full-wave physics",
            ],
        }
        _atomic_json(ADMISSION_REPORT, report)
        artifacts = [ADMISSION_REPORT, ADMISSION_FIGURE]
        _atomic_json(ADMISSION_MANIFEST, {
            "schema_version": "a4-c1-0-freeze-v1", "status": status,
            "source_sha256": sources,
            "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                                 for path in artifacts],
        })
        print(f"A4-C1-0: {'OK' if status.startswith('passed') else 'BLOCKED'} {status}")
        return 0 if status.startswith("passed") else 2
    except Exception as error:
        _atomic_json(ADMISSION_REPORT, {
            "status": "blocked_input_or_source_drift",
            "error": f"{type(error).__name__}: {error}"})
        print(f"A4-C1-0: BLOCKED {type(error).__name__}: {error}")
        return 2


def _run_matlab(input_path: Path) -> subprocess.CompletedProcess:
    quote = lambda value: str(value).replace("'", "''")
    angles = "[" + ",".join(f"{item[3]:.15g}" for item in CASES) + "]"
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); a4_c1_same_angle_pharlap("
        f"'{quote(input_path)}','{quote(RAW)}','{quote(PH_SUMMARY)}',{angles})")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run([str(MATLAB), "-batch", expression], env=env,
                            text=True, capture_output=True, timeout=1800,
                            errors="replace")
    PH_LOG.write_text((result.stdout or "") + "\n" + (result.stderr or ""),
                      encoding="utf-8")
    return result


def _ray_arrays(ray: dict) -> dict:
    keys = (
        "trajectory_ecef_km", "trajectory_group_tangent_ecef",
        "wavevector_ecef_km^-1", "trajectory_geocentric_lat_lon_height",
        "trajectory_physical_height_km", "trajectory_geom_path_km",
    )
    return {key: np.asarray(ray[key], dtype=float) for key in keys}


def _ray_summary(record: dict, ray: dict) -> dict:
    arrays = set(_ray_arrays(ray))
    return {**record, "raw_diagnostics": {
        key: value for key, value in ray.items() if key not in arrays
    }}


def _ah_worker(payload):
    case_id, input_path, elevation = payload
    context = _load_context(input_path)
    runs = []
    arrays = {}
    for repeat in (1, 2):
        record, ray = _trace(_model(context), context, elevation)
        runs.append(_ray_summary(record, ray))
        arrays.update({f"repeat{repeat}_{key}": value
                       for key, value in _ray_arrays(ray).items()})
    return case_id, runs, arrays


def _ah_paths(case_id: str) -> tuple[Path, Path]:
    return RAW / f"ah3d_{case_id}_summary.json", RAW / f"ah3d_{case_id}_raw.npz"


def _load_ah_case(case_id: str, sources: dict):
    summary_path, npz_path = _ah_paths(case_id)
    if not summary_path.is_file() or not npz_path.is_file():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("source_sha256") != sources
            or summary.get("raw_npz_sha256") != _sha256(npz_path)):
        raise RuntimeError(f"AH3D checkpoint source drift: {case_id}")
    with np.load(npz_path, allow_pickle=False) as source:
        arrays = {key: np.asarray(source[key], dtype=float) for key in source.files}
    return summary["runs"], arrays


def _save_ah_case(case_id: str, elevation: float, runs: list[dict], arrays: dict,
                  sources: dict) -> None:
    summary_path, npz_path = _ah_paths(case_id)
    np.savez_compressed(npz_path, **arrays)
    _atomic_json(summary_path, {
        "schema_version": "a4-c1-1-ah3d-fixed-angle-v1",
        "status": "completed_ah3d_fixed_angle_forward_runs",
        "case_semantics": "same_launch_angle_forward_trace_no_root_solving",
        "case_id": case_id, "elevation_deg": elevation, "runs": runs,
        "source_sha256": sources, "raw_npz_sha256": _sha256(npz_path),
    })


def _load_pharlap_case(path: Path, radius_km: float) -> list[dict]:
    paths = []
    with h5py.File(path, "r") as source:
        references = source["run_raw/ray_path_data"][()].reshape(-1)
        for reference in references:
            group = source[reference]
            lat = np.asarray(group["lat"], dtype=float).reshape(-1)
            lon = np.asarray(group["lon"], dtype=float).reshape(-1)
            height = np.asarray(group["height"], dtype=float).reshape(-1)
            geom = np.asarray(group["geometric_distance"], dtype=float).reshape(-1)
            geom -= geom[0]
            paths.append({
                "trajectory_ecef_km": geocentric_to_ecef(lat, lon, height, radius_km),
                "trajectory_geocentric_lat_lon_height": np.column_stack((lat, lon, height)),
                "trajectory_physical_height_km": height,
                "trajectory_geom_path_km": geom,
            })
    return paths


def _repeatability(runs: list[dict], paths: list[dict]) -> dict:
    scalar_ok = all(np.isclose(float(runs[0][key]), float(runs[1][key]),
                               rtol=1.0e-12, atol=1.0e-10)
                    for key in REPEAT_SCALARS)
    first = np.asarray(paths[0]["trajectory_ecef_km"], dtype=float)
    second = np.asarray(paths[1]["trajectory_ecef_km"], dtype=float)
    checks = {
        "termination": runs[0]["status"] == runs[1]["status"],
        "turning_count": int(runs[0]["turning_count"]) == int(runs[1]["turning_count"]),
        "point_count": int(runs[0]["path_point_count"]) == int(runs[1]["path_point_count"]),
        "critical_scalars": scalar_ok,
        "trajectory": first.shape == second.shape and np.allclose(
            first, second, rtol=1.0e-12, atol=1.0e-10),
    }
    return {"passed": all(checks.values()), "checks": checks}


def _resample(path: np.ndarray, values: np.ndarray, sample: np.ndarray) -> np.ndarray:
    distance = np.asarray(path, dtype=float).reshape(-1)
    distance -= distance[0]
    unique, indices = np.unique(distance, return_index=True)
    if unique.size < 2 or unique[-1] <= 0.0:
        raise ValueError("trajectory geometric path is not strictly usable")
    value = np.asarray(values, dtype=float)[indices]
    return np.column_stack([np.interp(sample, unique, value[:, axis])
                            for axis in range(value.shape[1])])


def _trajectory_metrics(ph: dict, ah: dict, radius_km: float) -> tuple[dict, dict]:
    ph_path = np.asarray(ph["trajectory_geom_path_km"], dtype=float)
    ah_path = np.asarray(ah["trajectory_geom_path_km"], dtype=float)
    u = np.linspace(0.0, 1.0, 201)
    ph_ecef = _resample(ph_path, ph["trajectory_ecef_km"], u * (ph_path[-1] - ph_path[0]))
    ah_ecef = _resample(ah_path, ah["trajectory_ecef_km"], u * (ah_path[-1] - ah_path[0]))
    delta = ah_ecef - ph_ecef
    norm = np.linalg.norm(delta, axis=1)
    geocentric = ecef_to_geocentric(ph_ecef, radius_km)
    enu = np.array([delta[index] @ enu_basis(row[0], row[1], radius_km)
                    for index, row in enumerate(geocentric)])
    horizontal = np.hypot(enu[:, 0], enu[:, 1])
    common = min(ph_path[-1] - ph_path[0], ah_path[-1] - ah_path[0])
    local_s = np.linspace(0.0, common, 201)
    local_delta = (_resample(ah_path, ah["trajectory_ecef_km"], local_s)
                   - _resample(ph_path, ph["trajectory_ecef_km"], local_s))
    metrics = {
        "trajectory_rms_3d_km": float(np.sqrt(np.mean(norm ** 2))),
        "trajectory_max_3d_km": float(np.max(norm)),
        "trajectory_east_rms_km": float(np.sqrt(np.mean(enu[:, 0] ** 2))),
        "trajectory_north_rms_km": float(np.sqrt(np.mean(enu[:, 1] ** 2))),
        "trajectory_up_rms_km": float(np.sqrt(np.mean(enu[:, 2] ** 2))),
        "trajectory_horizontal_rms_km": float(np.sqrt(np.mean(horizontal ** 2))),
        "common_geom_path_km": float(common),
        "common_path_diagnostic_rms_km": float(np.sqrt(np.mean(
            np.linalg.norm(local_delta, axis=1) ** 2))),
        "common_path_diagnostic_max_km": float(np.max(np.linalg.norm(local_delta, axis=1))),
    }
    residual = {"u": u, "total": norm, "east": enu[:, 0], "north": enu[:, 1],
                "up": enu[:, 2], "horizontal": horizontal}
    return metrics, residual


def _surface_distance(a: np.ndarray, b: np.ndarray, radius_km: float) -> float:
    a_rad = np.radians(a[:2]); b_rad = np.radians(b[:2])
    dot = (np.sin(a_rad[0]) * np.sin(b_rad[0])
           + np.cos(a_rad[0]) * np.cos(b_rad[0]) * np.cos(b_rad[1] - a_rad[1]))
    return float(radius_km * np.arccos(np.clip(dot, -1.0, 1.0)))


def _comparison(case, ph_runs, ph_paths, ah_runs, ah_paths,
                radius_km: float, target_range_km: float) -> tuple[dict, dict]:
    case_id, branch, anchor_solver, angle = case
    ph = ph_runs[0]; ah = ah_runs[0]
    ph_path = ph_paths[0]; ah_path = ah_paths[0]
    trajectory, residual = _trajectory_metrics(ph_path, ah_path, radius_km)
    landing = _surface_distance(
        np.asarray(ph_path["trajectory_geocentric_lat_lon_height"])[-1],
        np.asarray(ah_path["trajectory_geocentric_lat_lon_height"])[-1], radius_km)
    metric = {
        "case_id": case_id, "branch_id": branch, "anchor_solver": anchor_solver,
        "elevation_deg": angle,
        "pharlap_termination": ph["status"], "ah3d_termination": ah["status"],
        "pharlap_turning_count": int(ph["turning_count"]),
        "ah3d_turning_count": int(ah["turning_count"]),
        "landing_absolute_error_km": landing,
        "landing_relative_error": landing / target_range_km,
        "apex_height_absolute_error_km": abs(float(ah["apex_km"]) - float(ph["apex_km"])),
        "geometric_path_relative_error": abs(float(ah["geometric_path_km"])
            - float(ph["geometric_path_km"])) / abs(float(ph["geometric_path_km"])),
        "group_path_relative_error": abs(float(ah["group_path_km"])
            - float(ph["group_path_km"])) / abs(float(ph["group_path_km"])),
        "group_delay_relative_error": abs(float(ah["tau_ms"])
            - float(ph["tau_ms"])) / abs(float(ph["tau_ms"])),
        "pharlap_group_delay_closure_ms": abs(float(ph["tau_ms"])
            - 1000.0 * float(ph["group_path_km"]) / 299792.458),
        "ah3d_group_delay_closure_ms": abs(float(ah["tau_ms"])
            - 1000.0 * float(ah["group_path_km"]) / 299792.458),
        **trajectory,
    }
    gates = {
        "ground_single_turning": ph["status"] == ah["status"] == "ground"
            and int(ph["turning_count"]) == int(ah["turning_count"]) == 1,
        **{key: metric[key] <= threshold for key, threshold in LIMITS.items()},
        "group_delay_closure": metric["pharlap_group_delay_closure_ms"] <= 1.0e-10
            and metric["ah3d_group_delay_closure_ms"] <= 1.0e-10,
    }
    metric["gates"] = gates
    metric["passed"] = all(gates.values())
    return metric, residual


def _range_height(path: dict, tx: np.ndarray, radius_km: float) -> tuple[np.ndarray, np.ndarray]:
    geocentric = np.asarray(path["trajectory_geocentric_lat_lon_height"], dtype=float)
    ranges = np.array([_surface_distance(tx, point, radius_km) for point in geocentric])
    return ranges, geocentric[:, 2]


def _plot_same_angle(metrics: list[dict], case_data: dict, residuals: dict,
                     input_data: dict, source_tag: str) -> list[dict]:
    radius = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0
    tx = np.array([float(np.asarray(input_data[key]).item()) for key in
                   ("origin_lat_deg", "origin_lon_deg", "origin_height_km")])
    target = float(np.asarray(input_data["target_range_km"]).item())
    heights = np.asarray(input_data["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(input_data["iono_en_grid_m3"], dtype=float)
    profile = ne[ne.shape[0] // 2, ne.shape[1] // 2, :]

    figure, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for axis, metric in zip(axes.flat, metrics):
        image = axis.imshow(np.repeat(profile[:, None], 2, axis=1), origin="lower",
                            extent=(0, target * 1.06, heights[0], heights[-1]),
                            cmap="jet", alpha=0.5, aspect="auto")
        data = case_data[metric["case_id"]]
        for solver, style, color in (("PHaRLAP", "-", "#1368aa"),
                                     ("AH3D", "--", "#c43b3b")):
            path = data[solver.lower() + "_paths"][0]
            distance, height = _range_height(path, tx, radius)
            axis.plot(distance, height, style, color=color, linewidth=1.7, label=solver)
        axis.axvline(target, color="black", linestyle=":", linewidth=1)
        axis.set_title(f"{metric['branch_id']} | {metric['anchor_solver']} anchor | "
                       f"{metric['elevation_deg']:.8f} deg")
        axis.grid(alpha=0.2); axis.legend(fontsize=8)
    axes[1, 0].set_xlabel("Great-circle range from TX (km)")
    axes[1, 1].set_xlabel("Great-circle range from TX (km)")
    axes[0, 0].set_ylabel("Physical height (km)"); axes[1, 0].set_ylabel("Physical height (km)")
    figure.suptitle(
        "A4-C1-1 same-angle forward trajectories | one-way | 10 MHz | B=0 | collisions=off\n"
        f"spherical geocentric input / ECEF state | no root solving | source={source_tag} | scattering/power=N/A")
    figure.subplots_adjust(left=0.07, right=0.86, bottom=0.07, top=0.88,
                           hspace=0.20, wspace=0.20)
    colorbar_axis = figure.add_axes((0.89, 0.16, 0.015, 0.65))
    figure.colorbar(image, cax=colorbar_axis, label="Electron density Ne (m^-3)")
    figure.savefig(TRAJECTORY_FIGURE, dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
    for axis, metric in zip(axes.flat, metrics):
        values = residuals[metric["case_id"]]
        axis.plot(values["u"], values["total"], color="black", label="3-D norm")
        axis.plot(values["u"], values["horizontal"], color="#1368aa", label="horizontal")
        axis.plot(values["u"], np.abs(values["up"]), color="#c43b3b", label="|up|")
        axis.axhline(LIMITS["trajectory_max_3d_km"], color="black", linestyle="--",
                     linewidth=1, label="5 km max gate")
        axis.set_title(f"{metric['branch_id']} | {metric['anchor_solver']} anchor")
        axis.grid(alpha=0.2); axis.legend(fontsize=7)
    axes[1, 0].set_xlabel("Normalized cumulative geometric path u")
    axes[1, 1].set_xlabel("Normalized cumulative geometric path u")
    axes[0, 0].set_ylabel("AH3D - PHaRLAP residual (km)")
    axes[1, 0].set_ylabel("AH3D - PHaRLAP residual (km)")
    figure.suptitle(
        "A4-C1-1 normalized-path ECEF/ENU residuals | 201 samples | formal trajectory gate\n"
        f"one-way | 10 MHz | B=0 | collisions=off | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    figure.savefig(RESIDUAL_FIGURE, dpi=180)
    plt.close(figure)

    keys = ("landing_absolute_error_km", "apex_height_absolute_error_km",
            "geometric_path_relative_error", "group_path_relative_error",
            "trajectory_rms_3d_km", "trajectory_max_3d_km")
    labels = ("landing", "apex", "geom path", "group path", "trajectory RMS", "trajectory max")
    figure, axis = plt.subplots(figsize=(12, 6.5))
    x = np.arange(len(keys)); width = 0.19
    for index, metric in enumerate(metrics):
        ratios = [metric[key] / LIMITS[key] for key in keys]
        axis.bar(x + (index - 1.5) * width, ratios, width,
                 label=f"{metric['branch_id']} / {metric['anchor_solver']}",
                 color="#187b3f" if metric["passed"] else "#b22222", alpha=0.55 + 0.1 * index)
    axis.axhline(1.0, color="black", linestyle="--", label="acceptance limit")
    axis.set_xticks(x, labels)
    axis.set_ylabel("Metric / hard limit")
    axis.set_yscale("log"); axis.grid(axis="y", alpha=0.25); axis.legend(fontsize=8)
    axis.set_title(
        "A4-C1-1 same-angle metric gates | each branch/anchor evaluated independently\n"
        f"one-way | 10 MHz | B=0 | collisions=off | source={source_tag} | no averaging")
    figure.tight_layout()
    figure.savefig(METRICS_FIGURE, dpi=180)
    plt.close(figure)
    return [_image_check(path) for path in
            (TRAJECTORY_FIGURE, RESIDUAL_FIGURE, METRICS_FIGURE)]


def _write_metrics_csv(metrics: list[dict]) -> None:
    rows = []
    for metric in metrics:
        row = {key: value for key, value in metric.items()
               if not isinstance(value, (dict, list))}
        row.update({f"gate_{key}": value for key, value in metric["gates"].items()})
        rows.append(row)
    with METRICS_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def run_same_angle(workers: int = 2) -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    try:
        admission = _verify_manifest(ADMISSION_MANIFEST, "passed_c1_0_frozen_admission")
        input_path, _ = _refined_input()
        sources = _current_sources(input_path)
        if admission.get("source_sha256") != sources:
            raise RuntimeError("C1-0 source set changed before C1-1")
        process = _run_matlab(input_path)
        ph_summary = json.loads(PH_SUMMARY.read_text(encoding="utf-8")) \
            if PH_SUMMARY.is_file() else {}
        if (process.returncode != 0
                or ph_summary.get("status") != "passed_pharlap_fixed_angle_forward_runs"):
            raise RuntimeError(f"PHaRLAP fixed-angle run failed, returncode={process.returncode}")

        ah_sources = {key: sources[key] for key in (
            "refined_input", "ah_solver", "ah_provider", "ah_prep3_driver", "config")}
        ah_results = {}
        missing = []
        for case_id, _, _, elevation in CASES:
            loaded = _load_ah_case(case_id, ah_sources)
            if loaded is None:
                missing.append((case_id, str(input_path), elevation))
            else:
                ah_results[case_id] = loaded
        if missing:
            with ProcessPoolExecutor(max_workers=min(max(1, workers), len(missing))) as executor:
                futures = {executor.submit(_ah_worker, payload): payload for payload in missing}
                for future in as_completed(futures):
                    case_id, runs, arrays = future.result()
                    elevation = next(item[3] for item in CASES if item[0] == case_id)
                    _save_ah_case(case_id, elevation, runs, arrays, ah_sources)
                    ah_results[case_id] = (runs, arrays)
                    print(f"A4-C1-1 AH3D completed: {case_id}")

        input_data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
        radius = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0
        target = float(np.asarray(input_data["target_range_km"]).item())
        ph_cases = {item["case_id"]: item for item in ph_summary["cases"]}
        case_data = {}; metrics = []; residuals = {}; repeatability = {}
        for case in CASES:
            case_id = case[0]
            ph_record = ph_cases[case_id]
            ph_raw = Path(ph_record["raw_mat"])
            ph_paths = _load_pharlap_case(ph_raw, radius)
            ah_runs, ah_arrays = ah_results[case_id]
            ah_paths = [{key: ah_arrays[f"repeat{repeat}_{key}"] for key in (
                "trajectory_ecef_km", "trajectory_geocentric_lat_lon_height",
                "trajectory_physical_height_km", "trajectory_geom_path_km")}
                for repeat in (1, 2)]
            repeatability[case_id] = {
                "PHaRLAP": _repeatability(ph_record["runs"], ph_paths),
                "AH3D": _repeatability(ah_runs, ah_paths),
            }
            metric, residual = _comparison(
                case, ph_record["runs"], ph_paths, ah_runs, ah_paths, radius, target)
            metrics.append(metric); residuals[case_id] = residual
            case_data[case_id] = {"pharlap_paths": ph_paths, "ah3d_paths": ah_paths}
            _atomic_json(RAW / f"{case_id}_comparison.json", {
                "case_semantics": "same_angle_pair_no_root_solving",
                "metric": metric, "repeatability": repeatability[case_id],
                "pharlap_raw_mat": str(ph_raw),
                "ah3d_raw_npz": str(_ah_paths(case_id)[1]),
            })

        repeat_ok = all(value[solver]["passed"] for value in repeatability.values()
                        for solver in ("PHaRLAP", "AH3D"))
        _write_metrics_csv(metrics)
        _atomic_json(METRICS_JSON, {"limits": LIMITS, "metrics": metrics})
        source_hashes = {
            **sources, "admission_manifest": _sha256(ADMISSION_MANIFEST),
            "pharlap_summary": _sha256(PH_SUMMARY),
            **{f"pharlap_{case_id}_raw": _sha256(Path(ph_cases[case_id]["raw_mat"]))
               for case_id, *_ in CASES},
            **{f"ah3d_{case_id}_raw": _sha256(_ah_paths(case_id)[1])
               for case_id, *_ in CASES},
        }
        source_tag = hashlib.sha256(json.dumps(
            source_hashes, sort_keys=True).encode("ascii")).hexdigest()[:12]
        figures = _plot_same_angle(metrics, case_data, residuals, input_data, source_tag)
        checks = {
            "c1_0_admission_frozen": True,
            "four_symmetric_anchor_pairs": len(metrics) == 4,
            "eight_trajectories_plus_repeats_per_solver": all(
                len(ph_cases[item[0]]["runs"]) == 2 and len(ah_results[item[0]][0]) == 2
                for item in CASES),
            "no_root_solving": ph_summary["checks"]["no_root_solving"],
            "all_repeatability": repeat_ok,
            "all_semantic_ground_single_turning": all(
                metric["gates"]["ground_single_turning"] for metric in metrics),
            "all_numeric_hard_gates": all(metric["passed"] for metric in metrics),
            "all_visualizations": all(item["passed"] for item in figures),
            "raw_outputs_preserved": all(
                Path(ph_cases[item[0]]["raw_mat"]).is_file()
                and _ah_paths(item[0])[0].is_file() and _ah_paths(item[0])[1].is_file()
                for item in CASES),
        }
        if not repeat_ok or not checks["all_semantic_ground_single_turning"]:
            status = "failed_c1_1_semantic_or_repeatability"
        elif all(metric["passed"] for metric in metrics) and checks["all_visualizations"]:
            status = "passed_c1_1_same_angle_forward_comparison"
        else:
            status = "failed_c1_1_numeric_thresholds"
        report = {
            "schema_version": "a4-c1-1-validation-v1", "status": status,
            "case_semantics": (
                "same_launch_angle_forward_cross_solver_comparison_no_root_solving_"
                "not_independent_root_acceptance"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "checks": checks, "limits": LIMITS, "metrics": metrics,
            "repeatability": repeatability, "source_sha256": source_hashes,
            "source_set_short_hash": source_tag, "visualizations": figures,
            "artifacts": {"metrics_csv": str(METRICS_CSV),
                          "metrics_json": str(METRICS_JSON), "pharlap_log": str(PH_LOG)},
            "limitations": [
                "C1-1 tests fixed same-angle forward propagation, not P2P root agreement",
                "B=0 smooth Chapman F2 does not validate O/X, irregularities, scattering, PE, full-wave, or power",
                "Low/high labels are topology-consistent and do not directly establish variational Hessian type",
            ],
            "next_step": ("Proceed to C1-2 independent frozen-root comparison"
                          if status.startswith("passed") else
                          "Keep C1-2 diagnostic-only and prepare C1R attribution before solver changes"),
        }
        _atomic_json(SAME_REPORT, report)
        artifacts = [PH_SUMMARY, PH_LOG, METRICS_CSV, METRICS_JSON,
                     TRAJECTORY_FIGURE, RESIDUAL_FIGURE, METRICS_FIGURE, SAME_REPORT]
        artifacts.extend(Path(ph_cases[item[0]]["raw_mat"]) for item in CASES)
        artifacts.extend(path for item in CASES for path in _ah_paths(item[0]))
        artifacts.extend(RAW / f"{item[0]}_comparison.json" for item in CASES)
        _atomic_json(SAME_MANIFEST, {
            "schema_version": "a4-c1-1-freeze-v1", "status": status,
            "source_sha256": source_hashes,
            "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                                 for path in artifacts],
        })
        print(f"A4-C1-1: {'OK' if status.startswith('passed') else 'FAIL'} {status}")
        return 0 if status.startswith("passed") else 1
    except Exception as error:
        _atomic_json(SAME_REPORT, {"status": "blocked_c1_1", "error":
                                  f"{type(error).__name__}: {error}"})
        print(f"A4-C1-1: BLOCKED {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    arguments = sys.argv[1:]
    workers = 2
    stages = {"c1-0", "c1-1"}
    for argument in arguments:
        if argument.startswith("--workers="):
            workers = int(argument.split("=", 1)[1])
        elif argument.startswith("--stage="):
            stages = {argument.split("=", 1)[1]}
        else:
            raise SystemExit("usage: _run_a4_c1_validation.py [--stage=c1-0|c1-1] [--workers=N]")
    if not stages <= {"c1-0", "c1-1"}:
        raise SystemExit("stage must be c1-0 or c1-1")
    code = run_admission() if "c1-0" in stages else 0
    if code == 0 and "c1-1" in stages:
        code = run_same_angle(workers=workers)
    raise SystemExit(code)
