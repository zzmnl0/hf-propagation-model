"""A4-C1M staged PHaRLAP/AH3D attribution experiments through P2."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import copy
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.io import loadmat, savemat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import K_FP, UPDATE2  # noqa: E402
from models.dispersion_ah import C_KM_S  # noqa: E402
from models.ne_provider_geocentric import GeocentricGridNeProvider3D  # noqa: E402
from models.ray_tracer_ah_3d_ecef import (  # noqa: E402
    RefractiveIndexAH3DECEF,
    integrate_fixed_path_ah3d_ecef,
)
from utils import geocentric_to_ecef, launch_direction_ecef  # noqa: E402

BASE = ROOT / "output2" / "a4_pharlap"
OUT = BASE / "c1" / "c1m"
MATLAB_SOURCE = ROOT / "tests" / "matlab" / "a4_c1m_uniform_export.m"
MATLAB_EXE = Path(r"D:\Matlab2024a\bin\matlab.exe")
MEX = Path(r"D:\pharlap\pharlap_4.7.4\mex\raytrace_3d_sp.mexw64")
R_KM = 6371.0
S_KM = 100.0
DENSITIES = np.array([0.0, 1e10, 1e11, 5e11], dtype=np.float64)
FREQUENCIES = np.array([8.0, 10.0, 15.0], dtype=np.float64)
ELEVATIONS = np.array([22.0, 39.5], dtype=np.float64)
TOL = np.array([[1e-8, 0.01, 10.0], [1e-10, 0.001, 5.0],
                [1e-12, 0.001, 1.0]], dtype=np.float64)
SAMPLE = np.linspace(0.0, S_KM, 201)
LIMIT = {"rms_m": 10.0, "max_m": 50.0, "group_rel": 1e-5}
# A few ulps at terrestrial ECEF scale, converted from km to metres.
ROUND_FLOOR = {"rms_m": 1e-7, "max_m": 1e-7, "group_rel": 1e-13}
AH_PROFILES = (
    {"rtol": 1e-8, "atol_r_km": 1e-6, "atol_k": 1e-10,
     "ds_min_km": 0.01, "ds_max_km": 2.0},
    {"rtol": 1e-9, "atol_r_km": 1e-7, "atol_k": 1e-11,
     "ds_min_km": 0.005, "ds_max_km": 1.0},
    {"rtol": 1e-10, "atol_r_km": 1e-8, "atol_k": 1e-12,
     "ds_min_km": 0.002, "ds_max_km": 0.5},
    {"rtol": 1e-11, "atol_r_km": 1e-9, "atol_k": 1e-13,
     "ds_min_km": 0.001, "ds_max_km": 0.25},
    {"rtol": 1e-12, "atol_r_km": 1e-10, "atol_k": 1e-14,
     "ds_min_km": 0.0005, "ds_max_km": 0.125},
)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha(value: np.ndarray) -> str:
    value = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(value.dtype).encode("ascii"))
    digest.update(str(value.shape).encode("ascii"))
    digest.update(value.tobytes())
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    def scalar(item):
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, Path):
            return str(item)
        raise TypeError(f"unsupported JSON value: {type(item).__name__}")
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False, default=scalar) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_files() -> dict[str, Path]:
    c1 = BASE / "c1"
    prep1 = BASE / "c1_prep" / "prep1" / "a4_c1_prep1_validation.json"
    prep = read_json(prep1)
    refined = Path(next(row["path"] for row in prep["inputs"]
                        if row["variant"] == "refined_1km"))
    return {
        "runner": Path(__file__), "matlab_source": MATLAB_SOURCE,
        "plan": ROOT / "word" / "a4_c1m_implementation_plan.txt",
        "verifier": ROOT / "tests" / "_verify_a4_c1m_attribution.py",
        "matlab_executable": MATLAB_EXE,
        "ah_solver": ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        "ah_dispersion": ROOT / "models" / "dispersion_ah.py",
        "ah_provider": ROOT / "models" / "ne_provider_geocentric.py",
        "config": ROOT / "config.py", "pharlap_mex": MEX,
        "c1_prep1_report": prep1, "c1_refined_input": refined,
        "c1_admission_manifest": c1 / "admission" / "c1_0_freeze_manifest.json",
        "c1_same_angle_manifest": c1 / "same_angle" / "c1_1_freeze_manifest.json",
        "c1w_manifest": c1 / "c1w" / "a4_c1w_manifest.json",
        "c1w_input": c1 / "c1w" / "raw" / "a4_c1w_wgs84_ox_input.mat",
    }


def input_arrays() -> dict[str, np.ndarray | float]:
    lat = np.arange(28.0, 34.0 + 0.25, 0.5)
    lon = np.arange(118.0, 122.0 + 0.25, 0.5)
    height = np.arange(0.0, 400.0 + 2.5, 5.0)
    shape = (lat.size, lon.size, height.size)
    origin = np.asarray(geocentric_to_ecef(30.0, 120.0, 100.0, R_KM)) * 1000.0
    return {
        "latitude_deg": lat, "longitude_deg": lon, "height_km": height,
        "Ne_m3": np.broadcast_to(DENSITIES, shape + (4,)).copy(),
        "Bx_T": np.zeros(shape), "By_T": np.zeros(shape),
        "Bz_T": np.zeros(shape), "collision_freq_Hz": np.zeros(shape),
        "origin_lat_deg": 30.0, "origin_lon_deg": 120.0,
        "origin_height_km": 100.0, "earth_radius_m": R_KM * 1000.0,
        "densities_m3": DENSITIES, "frequencies_MHz": FREQUENCIES,
        "elevations_deg": ELEVATIONS, "tol_matrix": TOL,
        "origin_ecef_m": origin,
        "launch_unit_ecef": np.vstack([
            launch_direction_ecef(30.0, 120.0, float(e), 0.0, R_KM)
            for e in ELEVATIONS
        ]),
    }


def stage_a0(run: Path) -> None:
    files = source_files()
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing frozen source(s): " + ", ".join(missing))
    run.mkdir(parents=True, exist_ok=False)
    arrays = input_arrays()
    input_mat = run / "input.mat"
    savemat(input_mat, arrays, do_compression=True, oned_as="row")
    loaded = loadmat(input_mat, simplify_cells=True)
    for key, value in arrays.items():
        if not np.array_equal(np.asarray(value), np.asarray(loaded[key])):
            raise AssertionError(f"input MAT roundtrip mismatch: {key}")
    c1_arrays = {}
    for label in ("c1_refined_input", "c1w_input"):
        old = loadmat(files[label], simplify_cells=True)
        c1_arrays[label] = {key: array_sha(np.asarray(value)) for key, value in old.items()
                            if not key.startswith("__") and np.asarray(value).dtype.kind
                            in "biufc"}
    version_result = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         f"(Get-Item -LiteralPath '{MATLAB_EXE}').VersionInfo.ProductVersion"],
        capture_output=True, text=True, errors="replace", check=True,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    manifest = {
        "schema": "a4-c1m-a0-v1", "stage": "a0", "run_id": run.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "frozen_a0", "source_sha256": {k: sha(p) for k, p in files.items()},
        "matlab_executable_product_version": version_result.stdout.strip(),
        "source_path": {k: str(p) for k, p in files.items()},
        "c1_input_array_sha256": c1_arrays,
        "input_array_sha256": {k: array_sha(np.asarray(v)) for k, v in arrays.items()},
        "input_mat_sha256": sha(input_mat),
        "constants": {"common_absolute_arc_km": S_KM, "B_T": 0,
                      "collision_freq_Hz": 0, "earth_radius_km": R_KM,
                      "density_m3": DENSITIES.tolist(),
                      "frequency_MHz": FREQUENCIES.tolist(),
                      "elevation_deg": ELEVATIONS.tolist(),
                      "tol_matrix": TOL.tolist(), "ah_profiles": AH_PROFILES,
                      "limits": LIMIT, "self_convergence_fraction": 0.2,
                      "uncertainty_factor": 3,
                      "machine_roundoff_floor": ROUND_FLOOR},
    }
    write_json(run / "a0_manifest.json", manifest)
    (run / "a0_report.txt").write_text(
        "A0 冻结完成：24 个均匀介质案例，B=0、碰撞=0、球形地球 R=6371 km。\n"
        "输入 MAT 已回读核对；C1 基线、WGS84 输入与源代码记录 SHA-256。\n"
        "A1 需运行 PHaRLAP 与 AH3D 后才能判断数值收敛及物理误差。\n",
        encoding="utf-8")
    print(f"A0 frozen: {run}")
    print("PHaRLAP launch command: " + subprocess.list2cmdline(matlab_batch(run)))


def require_a0(run: Path) -> dict:
    manifest = read_json(run / "a0_manifest.json")
    if manifest.get("schema") != "a4-c1m-a0-v1" or manifest.get("status") != "frozen_a0":
        raise RuntimeError("A0 manifest is invalid")
    files = source_files()
    for key, original_hash in manifest["source_sha256"].items():
        if sha(files[key]) != original_hash:
            raise RuntimeError(f"frozen source changed since A0: {key}")
    if sha(run / "input.mat") != manifest["input_mat_sha256"]:
        raise RuntimeError("frozen input MAT changed since A0")
    return manifest


def matlab_batch(run: Path) -> list[str]:
    def q(path: Path) -> str:
        return str(path).replace("'", "''")

    expression = (f"addpath('{q(MATLAB_SOURCE.parent)}'); "
                  f"a4_c1m_uniform_export('{q(run / 'input.mat')}',"
                  f"'{q(run / 'pharlap_raw')}','{q(run / 'pharlap_summary.json')}')")
    return [str(MATLAB_EXE), "-batch", expression]


def stage_pharlap(run: Path) -> None:
    require_a0(run)
    summary = run / "pharlap_summary.json"
    raw_dir = run / "pharlap_raw"
    if summary.exists() or raw_dir.exists():
        raise FileExistsError("PHaRLAP raw/summary already exists; use a new run-id")
    if not MATLAB_EXE.is_file():
        raise FileNotFoundError(MATLAB_EXE)
    raw_dir.mkdir()
    command = matlab_batch(run)
    print("MATLAB command: " + subprocess.list2cmdline(command), flush=True)
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                            errors="replace", check=False,
                            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    (run / "pharlap_matlab.log").write_text(
        result.stdout + "\n[STDERR]\n" + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"MATLAB exited {result.returncode}; see pharlap_matlab.log")
    if read_json(summary).get("status") != "completed_uniform_export":
        raise RuntimeError("PHaRLAP summary does not report completion")
    print(f"PHaRLAP raw complete: {raw_dir}")


def mat(path: Path) -> dict:
    return loadmat(path, simplify_cells=True)


def vector(obj: dict, key: str) -> np.ndarray:
    value = np.asarray(obj[key], dtype=np.float64).reshape(-1)
    if value.size < 2 or not np.all(np.isfinite(value)):
        raise ValueError(f"invalid MATLAB vector {key}")
    return value


def scalar(obj: dict, key: str) -> float:
    value = float(np.asarray(obj[key], dtype=np.float64).reshape(-1)[0])
    if not np.isfinite(value):
        raise ValueError(f"invalid MATLAB scalar {key}")
    return value


def resample(x: np.ndarray, y: np.ndarray, sample: np.ndarray = SAMPLE,
             *, kind: str = "pchip") -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float)
    if y.shape[0] != x.size or np.any(~np.isfinite(y)):
        raise ValueError("sample coordinate/value mismatch")
    if x[0] > 1e-8 or x[-1] < sample[-1] - 1e-8:
        raise ValueError("ray does not cover common absolute 100 km arc")
    if np.any(np.diff(x) < -1e-9):
        raise ValueError("raw geometric path decreases")
    keep = np.r_[True, np.diff(x) > 1e-11]
    x, y = x[keep], y[keep]
    if np.any(np.diff(x) <= 0):
        raise ValueError("nonmonotone geometric path")
    if kind == "pchip":
        return np.asarray(PchipInterpolator(x, y, axis=0)(sample))
    if y.ndim == 1:
        return np.interp(sample, x, y)
    return np.column_stack([np.interp(sample, x, y[:, j])
                            for j in range(y.shape[1])])


def read_ph(path: Path, origin_ecef_m: np.ndarray) -> dict:
    raw = mat(path)
    state, raypath = raw["ray_state_vec"], raw["ray_path_data"]
    if not isinstance(state, dict) or not isinstance(raypath, dict):
        raise ValueError("PHaRLAP MAT must contain scalar state/path structs")
    position = np.column_stack([vector(state, key) for key in
                                ("pos_x", "pos_y", "pos_z")]) / 1000.0
    geometry = vector(state, "geometrical_path") / 1000.0
    group = vector(state, "group_path") / 1000.0
    phase_vector = np.column_stack([vector(state, key) for key in
                                    ("dir_x", "dir_y", "dir_z")])
    if geometry.size != position.shape[0] or group.size != geometry.size:
        raise ValueError("PHaRLAP state length mismatch")
    if np.linalg.norm(position[0] - origin_ecef_m / 1000.0) > 1e-6:
        raise ValueError("PHaRLAP initial ECEF position differs from frozen input")
    geometry -= geometry[0]
    group -= group[0]
    if np.any(np.diff(group) < -1e-8):
        raise ValueError("PHaRLAP group path decreases")
    lat, lon, height = (vector(raypath, key) for key in ("lat", "lon", "height"))
    if lat.size != geometry.size:
        raise ValueError("PHaRLAP path/state sample count mismatch")
    path_position = np.asarray(geocentric_to_ecef(lat, lon, height, R_KM))
    if path_position.shape != position.shape:
        raise ValueError("PHaRLAP path ECEF shape mismatch")
    path_state_max_m = float(np.max(np.linalg.norm(path_position - position, axis=1)) * 1000)
    path_group_max_m = float(np.max(np.abs(vector(raypath, "group_range")
                                           - group)) * 1000)
    ne_raw = vector(raypath, "electron_density") * 1e6
    n_raw = vector(raypath, "refractive_index")
    ng_raw = vector(raypath, "group_refractive_index")
    p2_raw = np.sum(phase_vector ** 2, axis=1)
    magnetic_max = max(float(np.max(np.abs(vector(raypath, key)))) for key in
                       ("geomag_x", "geomag_y", "geomag_z"))
    collision_max = float(np.max(np.abs(vector(raypath, "collision_frequency"))))
    xyz = resample(geometry, position)
    group_sample = resample(geometry, group)
    xyz_linear = resample(geometry, position, kind="linear")
    group_linear = resample(geometry, group, kind="linear")
    return {"s_raw_km": geometry, "xyz_raw_km": position,
            "group_raw_km": group, "xyz_km": xyz,
            "group_km": group_sample, "xyz_linear_km": xyz_linear,
            "group_linear_km": group_linear,
            "ne_m3": resample(geometry, ne_raw),
            "n": resample(geometry, n_raw),
            "ng": resample(geometry, ng_raw),
            "p2": resample(geometry, p2_raw),
            "phase_index_vector_first": phase_vector[0],
            "natural_end_label": float(np.asarray(raw["ray_data"]["ray_label"]).reshape(-1)[-1]),
            "height_100km_km": float(resample(geometry, height)[-1]),
            "height_initial_km": float(height[0]),
            "path_group_max_m": path_group_max_m,
            "magnetic_max": magnetic_max, "collision_max": collision_max,
            "path_state_max_m": path_state_max_m,
            "point_count": int(geometry.size),
            "raw_sha256": sha(path)}


def calibration(run: Path, summary: dict) -> dict:
    raw_path = run / "pharlap_raw" / "calibration_raw.mat"
    raw = mat(raw_path)
    raypath = raw["ray_path_data"]
    density = vector(raypath, "electron_density") * 1e6
    phase = vector(raypath, "refractive_index")
    inside = density > 0.9e11
    if (np.count_nonzero(inside) < 10 or
        np.max(np.abs(phase[inside] - 1.0)) > 1e-8):
        raise ValueError("vacuum-to-uniform calibration control changed")
    if summary.get("coefficient_identifiable") is not False:
        raise ValueError("uniform PH coefficient must be marked unidentifiable")
    if summary.get("calibration_k2") is not None:
        raise ValueError("uniform PH calibration_k2 must remain null")
    launch_k2 = float(summary["launch_k2"])
    if not np.isclose(launch_k2, K_FP ** 2, rtol=0, atol=1e-12):
        raise ValueError("PH launch K2 differs from controlled AH K2")
    return {"coefficient_identifiable": False, "pharlap_intrinsic_K2": None,
            "controlled_launch_K2": launch_k2, "ah_K2": K_FP ** 2,
            "uniform_calibration_probe_index": float(np.median(phase[inside])),
            "uniform_calibration_probe_density_m3": float(np.median(density[inside])),
            "sample_count": int(np.count_nonzero(inside)),
            "raw_sha256": sha(raw_path),
            "interpretation": "PHaRLAP preserves supplied initial index in uniform B=0 medium; no intrinsic coefficient inferred"}


def trace_ah(provider: GeocentricGridNeProvider3D, frequency: float,
             elevation: float, level: int) -> dict:
    name = f"c1m_level_{level}"
    UPDATE2["AH3D_INTEGRATION"][name] = dict(AH_PROFILES[level - 1])
    model = RefractiveIndexAH3DECEF(provider, earth_radius_km=R_KM,
                                     freq_MHz=frequency, wave_mode="O",
                                     gyro_MHz=0.0, collision_s=0.0,
                                     integration_profile=name)
    result = integrate_fixed_path_ah3d_ecef(
        model, tx_lat_deg=30.0, tx_lon_deg=120.0, tx_height_km=100.0,
        elevation_deg=elevation, bearing_deg=0.0, path_length_km=S_KM,
        profile=name, stop_on_physical_events=False)
    if result["termination_status"] != "fixed_geometric_path_completed":
        raise RuntimeError(f"AH3D did not finish fixed arc: {result['termination_status']}")
    xyz = resample(result["trajectory_geom_path_km"],
                   result["trajectory_ecef_km"])
    density = float(provider.Ne_m3.flat[0])
    n2 = 1 - K_FP ** 2 * density / (frequency * 1e6) ** 2
    if n2 <= 0:
        raise ValueError("nonpropagating AH3D uniform case")
    analytic_group = SAMPLE / np.sqrt(n2)
    native_group = float(result["group_path_km"])
    if abs(native_group - analytic_group[-1]) / analytic_group[-1] > 2e-7:
        raise RuntimeError("AH3D native group integral disagrees with uniform analytic check")
    xyz_linear = resample(result["trajectory_geom_path_km"],
                          result["trajectory_ecef_km"], kind="linear")
    return {"xyz_km": xyz, "xyz_linear_km": xyz_linear,
            "raw_xyz_km": result["trajectory_ecef_km"],
            "raw_wavevector_km_inv": result["wavevector_ecef_km^-1"],
            "raw_s_km": result["trajectory_geom_path_km"],
            "group_km": analytic_group,
            "group_native_km": native_group,
            "group_analytic_residual_rel":
            abs(native_group - analytic_group[-1]) / analytic_group[-1],
            "point_count": len(result["trajectory_geom_path_km"]),
            "model_provider_hash": provider.grid_hash}


def check_readback(run: Path) -> dict:
    frozen = mat(run / "input.mat")
    readback_path = run / "pharlap_raw" / "input_readback.mat"
    readback = mat(readback_path)
    keys = [key for key in frozen if not key.startswith("__")]
    for key in keys:
        if key not in readback or not np.array_equal(np.asarray(frozen[key]),
                                                      np.asarray(readback[key])):
            raise ValueError(f"MATLAB input readback mismatch: {key}")
    return {"readback_sha256": sha(readback_path), "matched_fields": keys}


def check_initialization_controls(run: Path, summary: dict) -> list[dict]:
    records = summary.get("initialization_records", [])
    if isinstance(records, dict):
        records = [records]
    if len(records) != 9:
        raise ValueError("expected nine explicit-state PHaRLAP initialization controls")
    expected = [(ki, mode) for ki in range(1, 4) for mode in (0, 1, -1)]
    origin = np.asarray(mat(run / "input.mat")["origin_ecef_m"], dtype=float).reshape(3)
    result = []
    for record, (ki, mode) in zip(records, expected):
        name = f"init_n{ki}_mode{mode}"
        path = run / "pharlap_raw" / f"{name}.mat"
        if (record["case_id"] != name or int(record["OX_mode"]) != mode
            or Path(record["raw_mat"]).resolve() != path.resolve()):
            raise ValueError(f"initialization control identity mismatch: {name}")
        raw = mat(path)
        first = np.array([vector(raw["ray_state_vec"], key)[0] for key in
                          ("dir_x", "dir_y", "dir_z")])
        initial = raw["initial"]
        supplied = np.array([scalar(initial, key) for key in
                             ("dir_x", "dir_y", "dir_z")])
        first_pos = np.array([vector(raw["ray_state_vec"], key)[0] for key in
                              ("pos_x", "pos_y", "pos_z")])
        if np.linalg.norm(first_pos - origin) > 1e-3:
            raise ValueError(f"initialization position mismatch: {name}")
        supplied_norm = float(np.linalg.norm(supplied))
        first_norm = float(np.linalg.norm(first))
        if (abs(supplied_norm - float(record["initial_index"])) > 1e-12
            or abs(first_norm - supplied_norm) > 1e-8
            or np.linalg.norm(first - supplied) > 1e-8):
            raise ValueError(f"initialization norm mismatch: {name}")
        result.append({"case_id": name, "OX_mode": mode,
                       "supplied_index": supplied_norm,
                       "exported_initial_index": first_norm,
                       "raw_sha256": sha(path)})
    return result


def case_keys() -> list[tuple[int, int, int]]:
    return [(ni, fi, ei) for ni in range(1, 5)
            for fi in range(1, 4) for ei in range(1, 3)]


def norm_metrics(delta_km: np.ndarray) -> dict:
    dist_m = np.linalg.norm(delta_km, axis=1) * 1000.0
    return {"rms_m": float(np.sqrt(np.mean(dist_m ** 2))),
            "max_m": float(np.max(dist_m))}


def per_case_metrics(ah: dict, ph: dict) -> dict:
    result = norm_metrics(ah["xyz_km"] - ph["xyz_km"])
    result["group_rel"] = abs(ah["group_native_km"] - ph["group_km"][-1]) / ph["group_km"][-1]
    result["ah_group_native_km"] = ah["group_native_km"]
    result["ph_group_100km_km"] = float(ph["group_km"][-1])
    result["ph_path_state_max_m"] = ph["path_state_max_m"]
    result["ph_path_group_max_m"] = ph["path_group_max_m"]
    result["ah_group_analytic_residual_rel"] = ah["group_analytic_residual_rel"]
    return result


def ph_uniform_checks(ph: dict, density: float, frequency: float,
                      launch_unit: np.ndarray) -> dict:
    n2 = 1.0 - K_FP ** 2 * density / (frequency * 1e6) ** 2
    ne_error = float(np.max(np.abs(ph["ne_m3"] - density)))
    ne_rel = ne_error / max(density, 1.0)
    n2_error = float(np.max(np.abs(ph["n"] ** 2 - n2)))
    p2_error = float(np.max(np.abs(ph["p2"] - n2)))
    ng_error = float(np.max(np.abs(ph["ng"] - 1 / np.sqrt(n2))))
    launch_error = float(np.linalg.norm(ph["phase_index_vector_first"]
                                        - np.sqrt(n2) * launch_unit))
    checks = {"ne": ne_rel <= 1e-10, "n2": n2_error <= 1e-8,
              "p2": p2_error <= 1e-8, "ng": ng_error <= 1e-8,
              "launch": launch_error <= 1e-8,
              "path_state": ph["path_state_max_m"] <= 1e-4,
              "path_group": ph["path_group_max_m"] <= 1e-4,
              "B_zero": ph["magnetic_max"] == 0.0,
              "collision_zero": ph["collision_max"] == 0.0,
              "upward_100km": ph["height_100km_km"] > ph["height_initial_km"]}
    return {"ne_relative_max": ne_rel, "n2_absolute_max": n2_error,
            "p2_absolute_max": p2_error, "ng_absolute_max": ng_error,
            "launch_vector_absolute": launch_error,
            "path_state_max_m": ph["path_state_max_m"],
            "path_group_max_m": ph["path_group_max_m"],
            "B_max": ph["magnetic_max"], "collision_max": ph["collision_max"],
            "natural_end_label": ph["natural_end_label"],
            "checks": checks, "pass": all(checks.values()),
            "interpretation": "controlled initial state consistency, not identification of PH intrinsic K2"}


def uncertainty(fine: dict, previous: dict) -> dict:
    ah, ph = fine["ah"], fine["ph"]
    ah_prev, ph_prev = previous["ah"], previous["ph"]
    ah_delta = norm_metrics(ah["xyz_km"] - ah_prev["xyz_km"])
    ph_delta = norm_metrics(ph["xyz_km"] - ph_prev["xyz_km"])
    method = norm_metrics((ah["xyz_km"] - ah["xyz_linear_km"])
                          - (ph["xyz_km"] - ph["xyz_linear_km"]))
    half_grid = np.linspace(0.0, S_KM, 401)
    ah_half = resample(ah["raw_s_km"], ah["raw_xyz_km"], half_grid)
    ph_half = resample(ph["s_raw_km"], ph["xyz_raw_km"], half_grid)
    coarse_error = norm_metrics(ah["xyz_km"] - ph["xyz_km"])
    half_error = norm_metrics(ah_half - ph_half)
    grid_delta = {key: abs(half_error[key] - coarse_error[key])
                  for key in ("rms_m", "max_m")}
    resampling = {key: max(method[key], grid_delta[key])
                  for key in ("rms_m", "max_m")}
    ah_delta["group_rel"] = (abs(ah["group_native_km"] - ah_prev["group_native_km"])
                              / ph["group_km"][-1])
    ph_delta["group_rel"] = (abs(ph["group_km"][-1] - ph_prev["group_km"][-1])
                              / ph["group_km"][-1])
    resampling["group_rel"] = (abs(ph["group_km"][-1]
                                    - ph["group_linear_km"][-1])
                               / ph["group_km"][-1])
    total = {key: ah_delta[key] + ph_delta[key] + resampling[key] for key in LIMIT}
    inflated = {key: max(ROUND_FLOOR[key], 3.0 * total[key]) for key in LIMIT}
    self_pass = {key: total[key] <= 0.2 * LIMIT[key] for key in LIMIT}
    return {"ah_final_increment": ah_delta, "ph_final_increment": ph_delta,
            "resampling_uncertainty": resampling, "self_convergence_sum": total,
            "half_interval_check": {"sample_km": 0.25,
                                    "half_grid_error": half_error,
                                    "coarse_grid_error": coarse_error,
                                    "metric_change": grid_delta,
                                    "interpolation_method_change": method},
            "self_convergence_pass": self_pass,
            "inflated_uncertainty": inflated,
            "roundoff_floor_applied": {key: inflated[key] == ROUND_FLOOR[key]
                                        for key in LIMIT}}


def write_samples(path: Path, cases: dict, levels: int) -> None:
    headers = ["case_id", "level", "s_km", "ah_x_km", "ah_y_km", "ah_z_km",
               "ph_x_km", "ph_y_km", "ph_z_km", "dx_m", "dy_m", "dz_m",
               "xyz_error_m", "ah_group_analytic_km", "ph_group_km",
               "ph_ne_m3", "ph_n", "ph_ng", "ph_p2"]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        for case_id, runs in cases.items():
            for level in range(1, levels + 1):
                if level not in runs:
                    continue
                ah, ph = runs[level]["ah"], runs[level]["ph"]
                difference_m = (ah["xyz_km"] - ph["xyz_km"]) * 1000.0
                magnitude_m = np.linalg.norm(difference_m, axis=1)
                for j, s in enumerate(SAMPLE):
                    writer.writerow([case_id, level, float(s),
                                     *ah["xyz_km"][j], *ph["xyz_km"][j],
                                     *difference_m[j], float(magnitude_m[j]),
                                     float(ah["group_km"][j]),
                                     float(ph["group_km"][j]),
                                     float(ph["ne_m3"][j]), float(ph["n"][j]),
                                     float(ph["ng"][j]), float(ph["p2"][j])])


def plot_difference(path: Path, cases: dict, levels: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for case_id, runs in cases.items():
        if levels not in runs:
            continue
        ah, ph = runs[levels]["ah"], runs[levels]["ph"]
        error_m = np.linalg.norm(ah["xyz_km"] - ph["xyz_km"], axis=1) * 1000.0
        axes[0].plot(SAMPLE, error_m, lw=0.8, alpha=0.6)
        ni = int(case_id.split("_")[0][1:])
        axes[1].scatter(DENSITIES[ni - 1] / 1e11,
                        (ah["group_native_km"] - ph["group_km"][-1]) * 1000,
                        color=f"C{ni - 1}", s=16)
    axes[0].axhline(LIMIT["rms_m"], color="black", ls="--", lw=1)
    axes[0].set(xlabel="Absolute geometric arc (km)", ylabel="ECEF difference (m)")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(xlabel="Ne (1e11 m^-3)", ylabel="AH3D - PH group path (m)")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def stage_a1(run: Path) -> None:
    frozen = require_a0(run)
    analysis = run / "a1"
    if analysis.exists():
        raise FileExistsError("A1 analysis exists; no overwrite")
    summary_path = run / "pharlap_summary.json"
    summary = read_json(summary_path)
    if summary.get("status") != "completed_uniform_export":
        raise RuntimeError("PHaRLAP export incomplete; see summary and MATLAB log")
    readback = check_readback(run)
    calibrated = calibration(run, summary)
    initialization = check_initialization_controls(run, summary)
    records = summary["records"]
    if isinstance(records, dict):
        records = [records]
    levels = len(records) // 24
    if len(records) != 24 * levels or not 3 <= levels <= 5:
        raise ValueError("expected 24 cases times 3-5 numerical levels")
    record_by_id = {record["case_id"]: record for record in records}
    if len(record_by_id) != len(records):
        raise ValueError("duplicate PHaRLAP case ids")
    raw_dir = run / "pharlap_raw"
    input_data = mat(run / "input.mat")
    origin = np.asarray(input_data["origin_ecef_m"], dtype=float).reshape(3)
    axes = [np.asarray(input_data[key], dtype=float).reshape(-1)
            for key in ("latitude_deg", "longitude_deg", "height_km")]
    densities = np.asarray(input_data["densities_m3"], dtype=float).reshape(-1)
    data: dict[str, dict[int, dict]] = {}
    case_reports = []
    raw_hashes = {}
    analysis.mkdir()
    for ni, fi, ei in case_keys():
        density = float(densities[ni - 1])
        frequency = float(FREQUENCIES[fi - 1])
        elevation = float(ELEVATIONS[ei - 1])
        case_id = f"n{ni}_f{fi}_e{ei}"
        provider = GeocentricGridNeProvider3D(
            np.asarray(input_data["Ne_m3"][:, :, :, ni - 1], dtype=float),
            *axes, earth_radius_km=R_KM, interpolation="linear")
        data[case_id] = {}
        case_failure = None
        for level in range(1, levels + 1):
            name = f"case_{case_id}_l{level}"
            record = record_by_id.get(name)
            if record is None:
                raise ValueError(f"missing PHaRLAP record: {name}")
            if (float(record["density_m3"]) != density or
                float(record["frequency_MHz"]) != frequency or
                float(record["elevation_deg"]) != elevation or
                int(record["level"]) != level or
                Path(record["raw_mat"]).resolve() != (raw_dir / f"{name}.mat").resolve()):
                raise ValueError(f"PHaRLAP record identity mismatch: {name}")
            try:
                ph = read_ph(raw_dir / f"{name}.mat", origin)
                ah = trace_ah(provider, frequency, elevation, level)
            except (ValueError, RuntimeError) as exc:
                raw = mat(raw_dir / f"{name}.mat")
                end_km = float(vector(raw["ray_state_vec"],
                                      "geometrical_path")[-1] / 1000.0)
                case_failure = {"level": level, "reason": str(exc),
                                "natural_end_label": record["ray_label"],
                                "natural_end_arc_km": end_km,
                                "covers_100km": end_km >= S_KM}
                break
            data[case_id][level] = {"ah": ah, "ph": ph}
            raw_hashes[name] = ph["raw_sha256"]
            np.savez_compressed(analysis / f"{name}_samples.npz", s_km=SAMPLE,
                                ah_xyz_km=ah["xyz_km"], ph_xyz_km=ph["xyz_km"],
                                ah_raw_s_km=ah["raw_s_km"],
                                ah_raw_xyz_km=ah["raw_xyz_km"],
                                ah_raw_wavevector_km_inv=ah["raw_wavevector_km_inv"],
                                ah_group_analytic_km=ah["group_km"],
                                ah_group_native_100km_km=ah["group_native_km"],
                                ph_group_km=ph["group_km"],
                                ph_ne_m3=ph["ne_m3"], ph_n=ph["n"],
                                ph_ng=ph["ng"], ph_p2=ph["p2"],
                                ph_raw_s_km=ph["s_raw_km"],
                                ph_raw_xyz_km=ph["xyz_raw_km"],
                                ph_raw_group_km=ph["group_raw_km"])
        if case_failure is not None:
            case_reports.append({"case_id": case_id, "density_m3": density,
                                 "frequency_MHz": frequency,
                                 "elevation_deg": elevation,
                                 "failure": case_failure, "pass": False})
            print(f"A1 {len(case_reports)}/24 {case_id}: failed "
                  f"{case_failure['reason']}", flush=True)
            continue
        finest = data[case_id][levels]
        metrics = per_case_metrics(finest["ah"], finest["ph"])
        direction = np.asarray(input_data["launch_unit_ecef"][ei - 1], dtype=float)
        physics = ph_uniform_checks(finest["ph"], density, frequency, direction)
        straight_xyz = origin[None, :] / 1000.0 + SAMPLE[:, None] * direction[None, :]
        n2 = 1 - K_FP ** 2 * density / (frequency * 1e6) ** 2
        straight_group = S_KM / np.sqrt(n2)
        analytic_control = {
            "ah_straight_ray": norm_metrics(finest["ah"]["xyz_km"] - straight_xyz),
            "ph_straight_ray": norm_metrics(finest["ph"]["xyz_km"] - straight_xyz),
            "ah_group_rel": abs(finest["ah"]["group_native_km"] - straight_group)
                            / straight_group,
            "ph_group_rel": abs(finest["ph"]["group_km"][-1] - straight_group)
                            / straight_group,
        }
        numeric = uncertainty(finest, data[case_id][levels - 1])
        budget = {key: {"physical_limit": metrics[key] <= LIMIT[key],
                        "numerical_explanation":
                        metrics[key] <= numeric["inflated_uncertainty"][key]}
                  for key in LIMIT}
        case_reports.append({"case_id": case_id, "density_m3": density,
                             "frequency_MHz": frequency, "elevation_deg": elevation,
                             "analytic_uniform_control": analytic_control,
                             "ph_uniform_controls": physics,
                             "levels": levels, "metrics_by_level":
                             {str(level): per_case_metrics(**data[case_id][level])
                              for level in range(1, levels + 1)},
                             "numeric": numeric, "physical_budget_pass": budget,
                             "pass": physics["pass"] and
                             all(all(gate.values()) for gate in budget.values()) and
                             all(numeric["self_convergence_pass"].values())})
        print(f"A1 {len(case_reports)}/24 {case_id}: "
              f"rms={metrics['rms_m']:.6g}m max={metrics['max_m']:.6g}m "
              f"group_rel={metrics['group_rel']:.6g}", flush=True)
    cache_records = summary["cache_records"]
    if isinstance(cache_records, dict):
        cache_records = [cache_records]
    if len(cache_records) != 4:
        raise ValueError("expected four alternating cache repeats")
    cache_results = []
    for index, ni in enumerate((2, 4, 2, 4), start=1):
        rec = cache_records[index - 1]
        name = f"cache_{index}"
        if rec["case_id"] != name or Path(rec["raw_mat"]).resolve() != (raw_dir / f"{name}.mat").resolve():
            raise ValueError("cache repeat record identity mismatch")
        ph = read_ph(raw_dir / f"{name}.mat", origin)
        raw_hashes[name] = ph["raw_sha256"]
        if levels not in data[f"n{ni}_f2_e1"]:
            cache_results.append({"case_id": name, "pass": False,
                                  "reason": "reference case did not cover common arc"})
            continue
        reference = data[f"n{ni}_f2_e1"][levels]["ph"]
        delta = norm_metrics(ph["xyz_km"] - reference["xyz_km"])
        group_rel = abs(ph["group_km"][-1] - reference["group_km"][-1]) / reference["group_km"][-1]
        cache_results.append({"case_id": name, **delta, "group_rel": group_rel,
                              "semantics": "full_input_reset_repeat_only",
                              "pass": delta["max_m"] <= 1e-5 and group_rel <= 1e-10})
    write_samples(analysis / "per_point.csv", data, levels)
    plot_difference(analysis / "difference.png", data, levels)
    needs_more = any("numeric" in row and not all(
        row["numeric"]["self_convergence_pass"].values()) for row in case_reports)
    has_early_failure = any("failure" in row for row in case_reports)
    passed = (all(row["pass"] for row in case_reports)
              and all(row["pass"] for row in cache_results))
    status = ("completed_a1_uniform_initialization_identified_constant_unidentifiable"
              if passed else ("needs_extra_levels" if needs_more and levels < 5
                              and not has_early_failure else "failed"))
    report = {"schema": "a4-c1m-a1-v1", "run_id": run.name,
              "status": status, "common_absolute_arc_km": S_KM,
        "comparison": "same controlled launch K2, first 100 km absolute geometric arc, no landing/P2P",
              "calibration": calibrated, "input_readback": readback,
              "initialization_controls": initialization,
              "limits": LIMIT, "machine_roundoff_floor": ROUND_FLOOR,
              "cases": case_reports, "cache_repeats": cache_results,
              "cache_scope": "full-grid reset/repeat; cached-grid shortcut not exercised",
              "pharlap_matlab_version": summary.get("matlab_version"),
              "pharlap_mex_path": summary.get("mex_path"),
              "pharlap_summary_sha256": sha(summary_path),
              "pharlap_raw_sha256": raw_hashes,
              "a0_manifest_sha256": sha(run / "a0_manifest.json")}
    write_json(analysis / "a1_report.json", report)
    text_lines = [
        f"A1 均匀介质归因：{status}",
        f"案例数：{len(case_reports)}；数值等级：{levels}；固定绝对几何弧长：100 km。",
        f"共同受控发射 K²={calibrated['controlled_launch_K2']:.12g}；PHaRLAP 内禀系数不可由均匀介质辨识。",
        "门槛：三维位置 RMS≤10 m、最大≤50 m、群程相对误差≤1e-5；",
        "自收敛预算为各门槛 20%，外推不确定度为三倍末级增量与重采样之和，",
        "机器舍入下限另行列示。判定同时要求误差低于物理门槛，且能由三倍数值增量解释。",
        f"未通过案例：{sum(not row['pass'] for row in case_reports)}；缓存重复异常："
        f"{sum(not row['pass'] for row in cache_results)}。",
    ]
    if status == "needs_extra_levels":
        text_lines.append("前三档不足以确认收敛；需追加合法的第 4/5 档原始追踪，当前不得判通过。")
    (analysis / "a1_report.txt").write_text("\n".join(text_lines) + "\n", encoding="utf-8")
    artifacts = sorted(path for path in analysis.iterdir() if path.is_file())
    write_json(analysis / "a1_manifest.json", {
        "schema": "a4-c1m-a1-freeze-v1", "status": status,
        "a0_source_sha256": frozen["source_sha256"],
        "artifacts": {path.name: sha(path) for path in artifacts}})
    print(f"A1 {status}: {analysis}")


from models.ne_provider_geocentric import GeocentricNeProvider3D


class _A2RadialProvider(GeocentricNeProvider3D):
    """Exact horizontal constant grid interpolation, with a vacuum floor."""
    name = "a2_radial_grid_ne_provider"

    def __init__(self, height: np.ndarray, profile: np.ndarray,
                 *, vacuum_below: bool, grid_hash: str,
                 interpolation: str = "linear"):
        self.height = np.asarray(height, dtype=float)
        self.profile = np.asarray(profile, dtype=float)
        self.earth_radius_km = R_KM
        self.floor = 0.0 if vacuum_below else float(self.height[0])
        self.vacuum_below = vacuum_below
        self.grid_hash = grid_hash
        self.interpolation = interpolation
        if self.height.ndim != 1 or self.profile.shape != self.height.shape:
            raise ValueError("invalid radial grid")
        if interpolation == "hermite":
            from scipy.interpolate import CubicHermiteSpline
            slopes = np.gradient(self.profile, self.height, edge_order=2)
            self._spline = CubicHermiteSpline(self.height, self.profile, slopes)
        elif interpolation != "linear":
            raise ValueError("unsupported radial interpolation")

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        if not np.isfinite(lat_deg + lon_deg + height_km):
            raise ValueError("nonfinite Ne query")
        if height_km < self.height[0] or height_km > self.height[-1]:
            return 0.0
        if self.interpolation == "hermite":
            return max(0.0, float(self._spline(height_km)))
        return float(np.interp(height_km, self.height, self.profile))

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                  height_km: float) -> float:
        return float(min(lat_deg - 28.0, 65.0 - lat_deg,
                         lon_deg - 118.0, 122.0 - lon_deg,
                         height_km - self.floor, self.height[-1] - height_km))

    def get_metadata(self) -> dict:
        return {**super().get_metadata(), "earth_radius_km": R_KM,
                "grid_hash": self.grid_hash,
                "vacuum_below_grid": self.vacuum_below,
                "interpolation": self.interpolation}


def _a2_profile(kind: str, height: np.ndarray, *, hmin: float = 0.0,
                gradient: float = 0.0) -> np.ndarray:
    if kind == "affine":
        return 1e11 + gradient * (height - 100.0)
    peak = 607635874822.0494  # frozen C1 input, foF2=7 MHz at hmF2=350 km
    u = (height - 350.0) / 100.0
    chapman = peak * np.exp(0.5 * (1.0 - u - np.exp(-u)))
    if kind == "boundary":
        return np.where(height >= hmin, chapman, 0.0)
    if kind == "smooth":
        x = np.clip((height - 20.0) / 60.0, 0.0, 1.0)
        return chapman * (x * x * x * (10.0 + x * (-15.0 + 6.0 * x)))
    if kind == "internal":
        return chapman
    raise ValueError(kind)


def _a2_job(case: str, stage: str, grid_path: Path, grid_hash: str,
            frequency: float, elevation: float, height: float,
            level: int, *, explicit: bool, launch_k2: float = 80.6404) -> dict:
    job_id = (f"{case}_f{frequency:g}_e{elevation:g}_t{level}"
              .replace(".", "p").replace("-", "m"))
    job = {
        "job_id": job_id, "case_id": case, "stage": stage,
        "grid_mat": str(grid_path.resolve()), "grid_sha256": grid_hash,
        "frequency_MHz": float(frequency), "elevation_deg": float(elevation),
        "origin_height_km": float(height), "explicit_state": bool(explicit),
        "launch_k2": float(launch_k2), "level": int(level),
        "tol": TOL[level - 1].copy(),
    }
    encoded = {key: value.tolist() if isinstance(value, np.ndarray) else value
               for key, value in job.items()}
    job["job_sha256"] = hashlib.sha256(
        json.dumps(encoded, sort_keys=True, separators=(",", ":"),
                   allow_nan=False).encode("utf-8")).hexdigest()
    return job


def stage_a2_prepare(run: Path) -> None:
    """Freeze PHaRLAP grids and independent, resumable A2 jobs."""
    a2 = run / "a2"
    if a2.exists():
        raise FileExistsError(f"A2 input already exists: {a2}")
    source = BASE / "c1_prep" / "prep1" / "inputs" / "a4_c1_chapman_refined_1km.mat"
    if not source.is_file():
        raise FileNotFoundError(source)
    original = mat(source)
    if not np.isclose(float(original["iono_en_grid_m3"][0, 0, 290]),
                      607635874822.0494, rtol=0.0, atol=1e-3):
        raise ValueError("frozen C1 Chapman peak changed")
    a2.mkdir(parents=True)
    grids = a2 / "grids"
    grids.mkdir()
    lat = np.arange(28.0, 65.0 + 0.25, 0.5)
    lon = np.arange(118.0, 122.0 + 0.25, 0.5)
    bh = np.arange(0.0, 460.0 + 5.0, 10.0)
    bshape = (lat.size, lon.size, bh.size)
    jobs: list[dict] = []
    grid_manifest = {}

    def add_grid(case: str, stage: str, height: np.ndarray,
                 profile: np.ndarray, origin: float,
                 frequencies: tuple[float, ...], explicit: bool,
                 elevations: tuple[float, ...] = (22.0, 39.5)) -> None:
        if np.any(profile < 0.0) or not np.all(np.isfinite(profile)):
            raise ValueError(f"invalid Ne profile: {case}")
        grid_path = grids / f"{case}.mat"
        shape = (lat.size, lon.size, height.size)
        grid = {
            "latitude_deg": lat, "longitude_deg": lon, "height_km": height,
            "Ne_m3": np.broadcast_to(profile, shape).copy(),
            "ne_profile_m3": profile,
            "geomag_latitude_deg": lat, "geomag_longitude_deg": lon,
            "geomag_height_km": bh,
            "Bx_T": np.zeros(bshape), "By_T": np.zeros(bshape),
            "Bz_T": np.zeros(bshape),
            "collision_freq_Hz": np.zeros(shape),
            "earth_radius_m": R_KM * 1000.0,
        }
        savemat(grid_path, grid, do_compression=True, oned_as="row")
        grid_hash = sha(grid_path)
        grid_manifest[case] = {
            "stage": stage, "grid_mat": str(grid_path.resolve()),
            "grid_sha256": grid_hash, "height_start_km": float(height[0]),
            "height_stop_km": float(height[-1]),
            "height_step_km": float(height[1] - height[0]),
            "Ne_at_start_m3": float(profile[0]),
            "Ne_at_origin_m3": float(np.interp(origin, height, profile)),
        }
        for frequency in frequencies:
            for elevation in elevations:
                for level in range(1, 4):
                    jobs.append(_a2_job(case, stage, grid_path, grid_hash,
                                        frequency, elevation, origin, level,
                                        explicit=explicit))

    frozen_angles = (21.94290326504532, 21.94562872092976,
                     39.54254824570964, 39.58701138438118)
    for hmin in (20.0, 40.0, 60.0, 80.0):
        h = np.arange(hmin, 460.0 + 1.0, 2.0)
        add_grid(f"boundary_h{int(hmin):03d}", "boundary", h,
                 _a2_profile("boundary", h, hmin=hmin), 0.0, (10.0,),
                 False, frozen_angles)
    h = np.arange(0.0, 460.0 + 1.0, 2.0)
    add_grid("smooth_s5", "smooth", h, _a2_profile("smooth", h),
             0.0, (10.0,), False, frozen_angles)
    add_grid("internal_common", "internal", h, _a2_profile("internal", h),
             100.0, (10.0,), True)
    for gradient_name, gradient in (("p1", 1e8), ("m1", -1e8),
                                    ("p2", 2e8), ("m2", -2e8)):
        for step in (2.0, 1.0, 0.5, 0.25):
            h = np.linspace(80.0, 180.0, int(round(100.0 / step)) + 1)
            step_name = str(step).replace(".", "p")
            case = f"affine_{gradient_name}_d{step_name}"
            add_grid(case, "affine", h,
                     _a2_profile("affine", h, gradient=gradient),
                     100.0, (10.0, 8.0, 15.0) if step == 0.25 else (10.0,), True)
            grid_manifest[case]["gradient_m3_per_km"] = gradient
    payload = np.empty((1, len(jobs)), dtype=object)
    payload[0, :] = jobs
    input_mat = a2 / "input.mat"
    savemat(input_mat, {"a2_jobs": payload}, do_compression=True,
            oned_as="row")
    manifest = {
        "schema": "a4-c1m-a2-input-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_c1_mat": str(source.resolve()), "source_c1_sha256": sha(source),
        "runner_sha256_at_prepare": sha(Path(__file__)),
        "input_mat": str(input_mat.resolve()), "input_mat_sha256": sha(input_mat),
        "grid": grid_manifest, "jobs": [
            {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in job.items()}
            for job in jobs],
    }
    write_json(a2 / "input_manifest.json", manifest)
    print(f"A2 prepared {len(jobs)} jobs: {input_mat}", flush=True)


def _a2_coefficient_from_ph(raw_path: Path, grid_path: Path,
                            frequency: float) -> dict:
    raw = mat(raw_path)
    state = raw["ray_state_vec"]
    path = raw["ray_path_data"]
    grid = mat(grid_path)
    height = np.asarray(path["height"], dtype=float).reshape(-1)
    axis = np.asarray(grid["height_km"], dtype=float).reshape(-1)
    profile = np.asarray(grid["ne_profile_m3"], dtype=float).reshape(-1)
    ne_input = np.interp(height, axis, profile, left=np.nan, right=np.nan)
    ne_ph = np.asarray(path["electron_density"], dtype=float).reshape(-1) * 1e6
    p2 = sum(np.asarray(state[key], dtype=float).reshape(-1) ** 2
             for key in ("dir_x", "dir_y", "dir_z"))
    if height.size != p2.size or not np.all(np.isfinite(p2)):
        raise ValueError(f"invalid PH state: {raw_path.name}")
    inside = np.isfinite(ne_input)
    readback = float(np.max(np.abs(ne_ph[inside] - ne_input[inside])))
    # Independent prescribed Ne field is the denominator; PH Ne is readback only.
    delta_ne = ne_input - ne_input[0]
    valid = inside & (np.abs(delta_ne) >= 1e8)
    if np.count_nonzero(valid) < 10:
        raise ValueError(f"insufficient gradient samples: {raw_path.name}")
    ceff = -(frequency * 1e6) ** 2 * (p2[valid] - p2[0]) / delta_ne[valid]
    return {
        "raw_mat": str(raw_path.resolve()), "raw_sha256": sha(raw_path),
        "grid_sha256": sha(grid_path), "sample_count": int(np.count_nonzero(valid)),
        "input_Ne_readback_max_m3": readback,
        "C_eff_median": float(np.median(ceff)),
        "C_eff_range": [float(np.min(ceff)), float(np.max(ceff))],
    }


def stage_a2_identified_prepare(run: Path) -> None:
    a2 = run / "a2"
    manifest = read_json(a2 / "input_manifest.json")
    summary = read_json(a2 / "pharlap_summary.json")
    if summary.get("status") != "completed_a2_export":
        raise ValueError("A2 PH export is incomplete")
    identified = a2 / "identified"
    if identified.exists():
        raise FileExistsError(identified)
    training = [job for job in manifest["jobs"]
                if job["stage"] == "affine" and job["frequency_MHz"] == 10.0
                and job["case_id"].endswith("d0p25") and job["level"] == 3]
    if len(training) != 8:
        raise ValueError("expected eight finest-grid training rays")
    evidence = []
    for job in training:
        raw_path = a2 / "pharlap_raw" / f"{job['job_id']}.mat"
        raw = mat(raw_path)
        if str(raw["job_sha256"]) != job["job_sha256"]:
            raise ValueError(f"PH job SHA mismatch: {job['job_id']}")
        evidence.append(_a2_coefficient_from_ph(
            raw_path, Path(job["grid_mat"]), job["frequency_MHz"]))
    identified_k2 = float(np.median([item["C_eff_median"] for item in evidence]))
    identified.mkdir()
    freeze = {
        "schema": "a4-c1m-a2-identified-k2-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection": "median of eight 10 MHz finest-grid level-3 affine-ray medians",
        "delta_Ne_min_m3": 1e8,
        "identified_K2": identified_k2,
        "training": evidence,
        "holdout_frequencies_MHz": [8.0, 15.0],
        "input_manifest_sha256": sha(a2 / "input_manifest.json"),
    }
    write_json(identified / "coefficient_freeze.json", freeze)
    grid_path = Path(manifest["grid"]["internal_common"]["grid_mat"])
    grid_hash = sha(grid_path)
    jobs = [_a2_job("internal_identified", "internal", grid_path, grid_hash,
                    10.0, float(elevation), 100.0, level, explicit=True,
                    launch_k2=identified_k2)
            for elevation in ELEVATIONS for level in range(1, 4)]
    payload = np.empty((1, len(jobs)), dtype=object)
    payload[0, :] = jobs
    input_mat = identified / "input.mat"
    savemat(input_mat, {"a2_jobs": payload}, do_compression=True, oned_as="row")
    write_json(identified / "input_manifest.json", {
        "schema": "a4-c1m-a2-identified-input-v1", "identified_K2": identified_k2,
        "coefficient_freeze_sha256": sha(identified / "coefficient_freeze.json"),
        "input_mat_sha256": sha(input_mat), "jobs": [
            {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in job.items()}
            for job in jobs],
    })
    print(f"A2 identified K2={identified_k2:.15g}; six PH jobs: {input_mat}", flush=True)


def stage_a2_identified_holdout_prepare(run: Path) -> None:
    a2 = run / "a2"
    base = read_json(a2 / "input_manifest.json")
    freeze_path = a2 / "identified" / "coefficient_freeze.json"
    freeze = read_json(freeze_path)
    identified_k2 = float(freeze["identified_K2"])
    output = a2 / "identified_holdout"
    if output.exists():
        raise FileExistsError(output)
    selected = [job for job in base["jobs"]
                if job["stage"] == "affine" and
                job["frequency_MHz"] in (8.0, 15.0) and
                job["case_id"].endswith("d0p25")]
    if len(selected) != 48:
        raise ValueError("expected 48 held-out affine jobs")
    jobs = [_a2_job(job["case_id"], "affine", Path(job["grid_mat"]),
                    job["grid_sha256"], job["frequency_MHz"],
                    job["elevation_deg"], 100.0, job["level"],
                    explicit=True, launch_k2=identified_k2)
            for job in selected]
    output.mkdir()
    payload = np.empty((1, len(jobs)), dtype=object)
    payload[0, :] = jobs
    input_mat = output / "input.mat"
    savemat(input_mat, {"a2_jobs": payload}, do_compression=True, oned_as="row")
    write_json(output / "input_manifest.json", {
        "schema": "a4-c1m-a2-identified-holdout-input-v1",
        "identified_K2": identified_k2,
        "coefficient_freeze_sha256": sha(freeze_path),
        "base_input_manifest_sha256": sha(a2 / "input_manifest.json"),
        "input_mat_sha256": sha(input_mat),
        "jobs": [{k: v.tolist() if isinstance(v, np.ndarray) else v
                  for k, v in job.items()} for job in jobs],
    })
    print(f"A2 identified holdout: {len(jobs)} PH jobs, {input_mat}", flush=True)


def stage_a2_extra_prepare(run: Path, *, dataset: str, job_id: str | None) -> None:
    raise ValueError("PHaRLAP MEX requires max_step >= 1 km; level 3 already uses the legal lower bound")


def _a2_ah_source_hashes() -> dict[str, str]:
    files = {
        "runner": Path(__file__), "config": ROOT / "config.py",
        "dispersion": ROOT / "models" / "dispersion_ah.py",
        "tracer": ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        "provider": ROOT / "models" / "ne_provider_geocentric.py",
    }
    return {key: sha(path) for key, path in files.items()}


def _a2_snapshot_ah_sources(a2: Path, source_hashes: dict[str, str]) -> Path:
    import shutil
    files = {
        "runner": Path(__file__), "config": ROOT / "config.py",
        "dispersion": ROOT / "models" / "dispersion_ah.py",
        "tracer": ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        "provider": ROOT / "models" / "ne_provider_geocentric.py",
    }
    snapshot = a2 / "sources" / "ah" / source_hashes["runner"][:16]
    snapshot.mkdir(parents=True, exist_ok=True)
    inventory = {}
    for key, source in files.items():
        target = snapshot / f"{key}{source.suffix}"
        if not target.exists():
            shutil.copyfile(source, target)
        if sha(target) != source_hashes[key]:
            raise ValueError(f"AH source snapshot mismatch: {key}")
        inventory[key] = {"source": str(source.resolve()),
                          "snapshot": str(target.resolve()),
                          "sha256": source_hashes[key]}
    manifest_path = snapshot / "source_manifest.json"
    if not manifest_path.exists():
        write_json(manifest_path, {"schema": "a4-c1m-a2-ah-source-snapshot-v1",
                                   "files": inventory})
    elif read_json(manifest_path)["files"] != inventory:
        raise ValueError("AH source snapshot manifest changed")
    return manifest_path


def _a2_ah_one(job: dict, out_dir: str, coefficient: float, interpolation: str,
               fd_base_km: float,
               source_hashes: dict[str, str]) -> dict:
    """Run one production AH3D integrator call in its own process."""
    start = time.perf_counter()
    job_id = job["job_id"]
    directory = Path(out_dir)
    result_path = directory / f"{job_id}.json"
    raw_path = directory / f"{job_id}.npz"
    result = {
        "schema": "a4-c1m-a2-ah-raw-v1", "job_id": job_id,
        "job_sha256": job["job_sha256"], "grid_sha256": job["grid_sha256"],
        "source_sha256": source_hashes, "coefficient_K2": coefficient,
        "interpolation": interpolation,
        "spatial_fd_base_km": fd_base_km,
    }
    try:
        grid_path = Path(job["grid_mat"])
        if sha(grid_path) != job["grid_sha256"]:
            raise ValueError("grid MAT SHA mismatch")
        grid = mat(grid_path)
        height = np.asarray(grid["height_km"], dtype=float).reshape(-1)
        profile = np.asarray(grid["ne_profile_m3"], dtype=float).reshape(-1)
        ne_grid = np.asarray(grid["Ne_m3"], dtype=float)
        if ne_grid.shape != (len(grid["latitude_deg"]), len(grid["longitude_deg"]),
                             len(height)):
            raise ValueError("grid shape mismatch")
        if (not np.array_equal(ne_grid[0, 0], profile)
                or not np.array_equal(ne_grid[-1, -1], profile)):
            raise ValueError("radial profile differs grid")
        provider = _A2RadialProvider(height, profile,
                                     vacuum_below=job["stage"] == "boundary",
                                     grid_hash=job["grid_sha256"],
                                     interpolation=interpolation)
        if job["stage"] == "boundary":
            if provider.get_Ne_geocentric(30.0, 120.0, height[0] - 1e-6) != 0.0:
                raise ValueError("boundary vacuum floor failed")
        from models import dispersion_ah
        from models.ray_tracer_ah_3d_ecef import trace_ray_ah3d_ecef
        dispersion_ah.K_FP = float(np.sqrt(coefficient))
        level = int(job["level"])
        name = f"a2_level_{level}"
        UPDATE2["AH3D_INTEGRATION"][name] = dict(AH_PROFILES[level - 1])
        model = RefractiveIndexAH3DECEF(
            provider, earth_radius_km=R_KM,
            freq_MHz=float(job["frequency_MHz"]), wave_mode="O",
            gyro_MHz=0.0, collision_s=0.0, integration_profile=name,
            spatial_fd_base_km=fd_base_km)
        common = dict(tx_lat_deg=float(job.get("origin_lat_deg", 30.0)),
                      tx_lon_deg=float(job.get("origin_lon_deg", 120.0)),
                      tx_height_km=float(job["origin_height_km"]),
                      elevation_deg=float(job["elevation_deg"]),
                      bearing_deg=0.0, profile=name)
        if job["stage"] in ("boundary", "smooth"):
            traced = trace_ray_ah3d_ecef(model, **common,
                                         h_stop_km=460.0, max_path_km=5000.0)
        else:
            traced = integrate_fixed_path_ah3d_ecef(
                model, **common, path_length_km=100.0,
                stop_on_physical_events=False,
                apply_dispersion_projection=bool(job.get("apply_projection", True)))
        xyz = np.asarray(traced["trajectory_ecef_km"], dtype=float)
        s = np.asarray(traced["trajectory_geom_path_km"], dtype=float)
        wavevector = np.asarray(traced["wavevector_ecef_km^-1"], dtype=float)
        coords = np.asarray(traced["trajectory_lat_lon_height"], dtype=float)
        tangent = np.asarray(traced["trajectory_group_tangent_ecef"], dtype=float)
        ne_native = np.array([provider.get_Ne_geocentric(*row) for row in coords])
        n2 = 1.0 - coefficient * ne_native / (float(job["frequency_MHz"]) * 1e6) ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            group_check = float(np.trapezoid(1.0 / np.sqrt(n2), s))
        tmp_path = raw_path.with_suffix(".tmp")
        with tmp_path.open("wb") as stream:
            np.savez_compressed(stream, xyz_km=xyz, geom_path_km=s,
                                wavevector_km_inv=wavevector,
                                group_tangent_ecef=tangent, lat_lon_height=coords,
                                Ne_m3_native=ne_native, n2_native=n2)
        tmp_path.replace(raw_path)
        result.update({
            "status": ("completed" if traced["termination_status"] in
                       ("ground", "upper_boundary", "grid_exit",
                        "fixed_geometric_path_completed") else
                       "physics_domain_rejected" if
                       traced["termination_status"] == "physics_domain_rejected" else
                       "incomplete_max_path"),
            "termination_status": traced["termination_status"],
            "failure_reason": traced.get("failure_reason"),
            "raw_npz": str(raw_path.resolve()), "raw_npz_sha256": sha(raw_path),
            "point_count": int(s.size), "geom_path_km": float(traced["geom_path_km"]),
            "group_path_km_native": float(traced["group_path_km"]),
            "group_path_km_trapezoid_diagnostic": group_check,
            "group_trapezoid_relative_difference":
                float((group_check - traced["group_path_km"]) /
                      traced["group_path_km"]),
            "event_history": traced["event_history"],
            "event_localization_error_max_km": traced["event_localization_error_max_km"],
            "initial_n2_prescribed": float(n2[0]),
        })
    except Exception as exc:
        result.update({"status": "failed", "error_type": type(exc).__name__,
                       "error": str(exc)})
    result["elapsed_s"] = time.perf_counter() - start
    write_json(result_path, result)
    return result


def stage_a2_ah(run: Path, *, dataset: str, subset: str,
                workers: int, requested_k2: float | None, interpolation: str,
                job_id: str | None, tag: str, extra_id: str | None,
                fd_base_km: float, level: int | None) -> None:
    a2 = run / "a2"
    root = (a2 / "extra" / extra_id if dataset == "extra" and extra_id else
            a2 / dataset if dataset != "base" else a2)
    if dataset == "extra" and not extra_id:
        raise ValueError("extra dataset requires --a2-extra-id")
    manifest = read_json(root / "input_manifest.json")
    input_mat = root / "input.mat"
    if sha(input_mat) != manifest["input_mat_sha256"]:
        raise ValueError("A2 input MAT changed")
    coefficient = (float(requested_k2) if requested_k2 is not None else
                   float(manifest["identified_K2"]) if dataset != "base" else
                   K_FP ** 2)
    if not np.isfinite(coefficient) or coefficient <= 0.0:
        raise ValueError("invalid coefficient")
    label = "common" if np.isclose(coefficient, K_FP ** 2, rtol=0.0, atol=1e-12) else "identified"
    output = root / "ah" / tag / label / interpolation
    if fd_base_km != 0.5:
        output = output / f"fd{fd_base_km:g}"
    output.mkdir(parents=True, exist_ok=True)
    source_hashes = _a2_ah_source_hashes()
    source_snapshot = _a2_snapshot_ah_sources(a2, source_hashes)
    selected = [job for job in manifest["jobs"]
                if (subset == "all" or job["stage"] == subset)
                and (job_id is None or job["job_id"] == job_id)
                and (level is None or job["level"] == level)]
    if not selected:
        raise ValueError("no matching A2 jobs")
    pending = []
    for job in selected:
        result_path = output / f"{job['job_id']}.json"
        if result_path.is_file():
            prior = read_json(result_path)
            raw_path = output / f"{job['job_id']}.npz"
            if (prior.get("status") in ("completed", "physics_domain_rejected",
                                        "incomplete_max_path")
                    and prior.get("job_sha256") == job["job_sha256"]
                    and prior.get("grid_sha256") == job["grid_sha256"]
                    and prior.get("source_sha256") == source_hashes
                    and prior.get("coefficient_K2") == coefficient
                    and prior.get("interpolation") == interpolation
                    and prior.get("spatial_fd_base_km") == fd_base_km
                    and raw_path.is_file()
                    and prior.get("raw_npz_sha256") == sha(raw_path)):
                continue
            raise RuntimeError(f"existing AH result mismatches frozen input/source: {job['job_id']}")
        pending.append(job)
    print(f"A2 AH {dataset}/{subset}: {len(selected)} selected, {len(pending)} pending, K2={coefficient:.15g}",
          flush=True)
    print(f"A2 AH source snapshot: {source_snapshot}", flush=True)
    results = []
    if workers == 1:
        for job in pending:
            result = _a2_ah_one(job, str(output), coefficient,
                                interpolation, fd_base_km, source_hashes)
            results.append(result)
            print(f"{job['job_id']}: {result['status']} {result['elapsed_s']:.1f}s", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_a2_ah_one, job, str(output), coefficient,
                                   interpolation, fd_base_km,
                                   source_hashes) for job in pending]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                if result["status"] != "completed" or len(results) % 10 == 0:
                    print(f"A2 AH {len(results)}/{len(pending)}: {result['job_id']} {result['status']}",
                          flush=True)
    failures = [result for result in results if result["status"] != "completed"]
    print(f"A2 AH batch finished: {len(results)} run, {len(failures)} failed, output={output}",
          flush=True)


def _a2_ph_trace(raw_path: Path, job: dict) -> dict:
    raw = mat(raw_path)
    if (str(raw["job_sha256"]) != job["job_sha256"]
            or str(raw["grid_sha256"]) != job["grid_sha256"]):
        raise ValueError(f"PH raw hash metadata mismatch: {job['job_id']}")
    state = raw["ray_state_vec"]
    path = raw["ray_path_data"]
    xyz = np.column_stack([vector(state, key) for key in
                           ("pos_x", "pos_y", "pos_z")]) / 1000.0
    s = vector(state, "geometrical_path") / 1000.0
    group = vector(state, "group_path") / 1000.0
    s -= s[0]
    group -= group[0]
    if not np.all(np.diff(s) > 0.0):
        raise ValueError(f"PH geometric path is not increasing: {job['job_id']}")
    prescribed = np.asarray(geocentric_to_ecef(
        30.0, 120.0, float(job["origin_height_km"]), R_KM), dtype=float)
    origin_error_m = float(np.linalg.norm(xyz[0] - prescribed) * 1000.0)
    if origin_error_m > 1e-3:
        raise ValueError(f"PH origin mismatch: {job['job_id']} ({origin_error_m:g}m)")
    p2 = sum(vector(state, key) ** 2 for key in ("dir_x", "dir_y", "dir_z"))
    return {"xyz": xyz, "s": s, "group": group, "p2": p2,
            "height": vector(path, "height"),
            "ne_m3": vector(path, "electron_density") * 1e6,
            "ray_label": float(np.asarray(raw["ray_data"]["ray_label"]).reshape(-1)[-1]),
            "origin_error_m": origin_error_m, "raw_sha256": sha(raw_path)}


def _a2_compare(job: dict, dataset_dir: Path, ah_dir: Path,
                interpolation: str) -> dict:
    raw_path = dataset_dir / "pharlap_raw" / f"{job['job_id']}.mat"
    ah_path = ah_dir / f"{job['job_id']}.json"
    ph = _a2_ph_trace(raw_path, job)
    grid = mat(Path(job["grid_mat"]))
    height = np.asarray(grid["height_km"], dtype=float).reshape(-1)
    profile = np.asarray(grid["ne_profile_m3"], dtype=float).reshape(-1)
    provider = _A2RadialProvider(height, profile,
                                 vacuum_below=job["stage"] == "boundary",
                                 grid_hash=job["grid_sha256"],
                                 interpolation=interpolation)
    field_ne = np.array([provider.get_Ne_geocentric(30.0, 120.0, h)
                         for h in ph["height"]])
    interior = ((ph["height"] >= height[0] + (height[1] - height[0]))
                & (ph["height"] <= height[-1]))
    first_interval = ((ph["height"] >= height[0]) &
                      (ph["height"] < height[1]))
    below = ph["height"] < height[0]
    row = {
        "job_id": job["job_id"], "case_id": job["case_id"],
        "stage": job["stage"], "frequency_MHz": job["frequency_MHz"],
        "elevation_deg": job["elevation_deg"], "level": job["level"],
        "ph_raw_sha256": ph["raw_sha256"], "ph_ray_label": ph["ray_label"],
        "ph_point_count": int(ph["s"].size),
        "ph_geom_path_km": float(ph["s"][-1]),
        "ph_group_path_km": float(ph["group"][-1]),
        "ph_origin_error_m": ph["origin_error_m"],
        "ph_Ne_below_grid_max_m3": float(np.max(np.abs(ph["ne_m3"][below])))
            if np.any(below) else None,
        "ph_Ne_interior_rms_m3": float(np.sqrt(np.mean(
            (ph["ne_m3"][interior] - field_ne[interior]) ** 2)))
            if np.any(interior) else None,
        "ph_Ne_first_interval_rms_m3": float(np.sqrt(np.mean(
            (ph["ne_m3"][first_interval] - field_ne[first_interval]) ** 2)))
            if np.any(first_interval) else None,
    }
    if job["stage"] == "affine":
        row["coefficient_probe"] = _a2_coefficient_from_ph(
            raw_path, Path(job["grid_mat"]), float(job["frequency_MHz"]))
    if not ah_path.exists():
        row["ah_status"] = "missing"
        return row
    ah = read_json(ah_path)
    if (ah.get("job_sha256") != job["job_sha256"] or
            ah.get("grid_sha256") != job["grid_sha256"]):
        raise ValueError(f"AH job/grid hash mismatch: {job['job_id']}")
    row.update({"ah_status": ah["status"],
                "ah_termination_status": ah.get("termination_status"),
                "ah_failure_reason": ah.get("failure_reason", ah.get("error")),
                "ah_elapsed_s": ah["elapsed_s"],
                "ah_raw_npz_sha256": ah.get("raw_npz_sha256"),
                "ah_geom_path_km": ah.get("geom_path_km"),
                "ah_group_path_km_native": ah.get("group_path_km_native"),
                "ah_group_trapezoid_relative_difference":
                    ah.get("group_trapezoid_relative_difference"),
                "ph_initial_p2": float(ph["p2"][0]),
                "ah_initial_n2": ah.get("initial_n2_prescribed")})
    if ah.get("initial_n2_prescribed") is not None:
        row["initial_n2_difference"] = abs(
            ph["p2"][0] - ah["initial_n2_prescribed"])
        row["same_initial_state"] = row["initial_n2_difference"] <= 1e-10
    if ah["status"] != "completed":
        return row
    raw_npz = ah_dir / f"{job['job_id']}.npz"
    if sha(raw_npz) != ah["raw_npz_sha256"]:
        raise ValueError(f"AH raw NPZ SHA mismatch: {job['job_id']}")
    with np.load(raw_npz) as sample:
        ah_s = np.array(sample["geom_path_km"])
        ah_xyz = np.array(sample["xyz_km"])
    fixed = job["stage"] in ("internal", "affine")
    end_s = 100.0 if fixed else min(float(ph["s"][-1]), float(ah_s[-1]))
    if end_s > ph["s"][-1] + 1e-9 or end_s > ah_s[-1] + 1e-9:
        raise ValueError(f"missing 100-km comparison prefix: {job['job_id']}")
    sample_s = np.linspace(0.0, end_s, 201 if fixed else 501)
    delta_m = np.linalg.norm(
        resample(ph["s"], ph["xyz"], sample_s) -
        resample(ah_s, ah_xyz, sample_s), axis=1) * 1000.0
    ph_group = (float(resample(ph["s"], ph["group"], np.array([end_s]))[0])
                if fixed else float(ph["group"][-1]))
    ah_group = float(ah["group_path_km_native"])
    row.update({
        "comparison_arc_km": end_s, "comparison_fixed_100km": fixed,
        "position_rms_m": float(np.sqrt(np.mean(delta_m ** 2))),
        "position_max_m": float(np.max(delta_m)),
        "position_endpoint_m": float(np.linalg.norm(
            (ph["xyz"][-1] if not fixed else
             resample(ph["s"], ph["xyz"], np.array([100.0]))[0]) -
            ah_xyz[-1]) * 1000.0),
        "group_comparison_kind": "fixed_arc_100km" if fixed else "natural_terminal",
        "ph_group_for_gate_km": ph_group,
        "group_relative_error": float(abs(ah_group - ph_group) / ph_group),
    })
    return row


def _a2_numeric_gate(levels: list[dict], jobs_by_id: dict[str, dict],
                     dataset_dir: Path, ah_dir: Path) -> dict:
    """Independent level increments and interpolation budget for one case."""
    if len(levels) != 3 or any(row["ah_status"] != "completed" for row in levels):
        return {"status": "native_incomplete"}
    if any(not row.get("same_initial_state", False) for row in levels):
        return {"status": "initial_state_mismatch"}
    traces = []
    for row in levels:
        job = jobs_by_id[row["job_id"]]
        ph = _a2_ph_trace(dataset_dir / "pharlap_raw" / f"{row['job_id']}.mat", job)
        ah = read_json(ah_dir / f"{row['job_id']}.json")
        with np.load(ah_dir / f"{row['job_id']}.npz") as raw:
            traces.append({"ph": ph, "ah": ah,
                           "ah_s": np.array(raw["geom_path_km"]),
                           "ah_xyz": np.array(raw["xyz_km"])})
    fixed = levels[0]["stage"] in ("internal", "affine")

    def endpoint_ph(trace: dict) -> tuple[np.ndarray, float]:
        ph = trace["ph"]
        if fixed:
            return (resample(ph["s"], ph["xyz"], np.array([100.0]))[0],
                    float(resample(ph["s"], ph["group"], np.array([100.0]))[0]))
        return ph["xyz"][-1], float(ph["group"][-1])

    def increment(a: dict, b: dict) -> dict:
        result = {}
        for side in ("ph", "ah"):
            if side == "ph":
                s1, xyz1 = a["ph"]["s"], a["ph"]["xyz"]
                s2, xyz2 = b["ph"]["s"], b["ph"]["xyz"]
            else:
                s1, xyz1 = a["ah_s"], a["ah_xyz"]
                s2, xyz2 = b["ah_s"], b["ah_xyz"]
            end = 100.0 if fixed else min(float(s1[-1]), float(s2[-1]))
            arc = np.linspace(0.0, end, 201 if fixed else 501)
            diff = resample(s1, xyz1, arc) - resample(s2, xyz2, arc)
            result[side] = norm_metrics(diff)
        ph1_xyz, ph1_group = endpoint_ph(a)
        ph2_xyz, ph2_group = endpoint_ph(b)
        result["ph"]["landing_m"] = float(np.linalg.norm(ph1_xyz - ph2_xyz) * 1000.0)
        result["ph"]["group_rel"] = abs(ph1_group - ph2_group) / ph2_group
        result["ah"]["landing_m"] = float(np.linalg.norm(
            a["ah_xyz"][-1] - b["ah_xyz"][-1]) * 1000.0)
        result["ah"]["group_rel"] = abs(
            a["ah"]["group_path_km_native"] - b["ah"]["group_path_km_native"]
        ) / ph2_group
        return result

    d12 = increment(traces[0], traces[1])
    d23 = increment(traces[1], traces[2])
    fine = traces[2]
    ph, ah = fine["ph"], fine["ah"]
    end = 100.0 if fixed else min(float(ph["s"][-1]), float(fine["ah_s"][-1]))

    def positions(step: float, kind: str) -> tuple[np.ndarray, np.ndarray]:
        arc = np.linspace(0.0, end, int(np.ceil(end / step)) + 1)
        return (resample(ph["s"], ph["xyz"], arc, kind=kind),
                resample(fine["ah_s"], fine["ah_xyz"], arc, kind=kind))

    ph05, ah05 = positions(0.5, "pchip")
    ph025, ah025 = positions(0.25, "pchip")
    phlin, ahlin = positions(0.5, "linear")
    coarse = norm_metrics(ph05 - ah05)
    refined = norm_metrics(ph025 - ah025)
    method = norm_metrics((ph05 - ah05) - (phlin - ahlin))
    resampling = {
        "rms_m": max(abs(refined["rms_m"] - coarse["rms_m"]), method["rms_m"]),
        "max_m": max(abs(refined["max_m"] - coarse["max_m"]), method["max_m"]),
        "landing_m": 0.0,
        "group_rel": 0.0,
    }
    if fixed:
        ph_at_end = resample(ph["s"], ph["xyz"], np.array([100.0]))[0]
        ph_linear = resample(ph["s"], ph["xyz"], np.array([100.0]),
                             kind="linear")[0]
        resampling["landing_m"] = float(np.linalg.norm(ph_at_end - ph_linear) * 1000.0)
        ph_group = endpoint_ph(fine)[1]
        linear_group = float(resample(ph["s"], ph["group"],
                                      np.array([100.0]), kind="linear")[0])
        resampling["group_rel"] = abs(ph_group - linear_group) / ph_group
    limits = {"rms_m": LIMIT["rms_m"], "max_m": LIMIT["max_m"],
              "landing_m": LIMIT["max_m"], "group_rel": LIMIT["group_rel"]}
    actual = {"rms_m": levels[2]["position_rms_m"],
              "max_m": levels[2]["position_max_m"],
              "landing_m": levels[2]["position_endpoint_m"],
              "group_rel": levels[2]["group_relative_error"]}
    floor = {"rms_m": ROUND_FLOOR["rms_m"],
             "max_m": ROUND_FLOOR["max_m"],
             "landing_m": ROUND_FLOOR["max_m"],
             "group_rel": ROUND_FLOOR["group_rel"]}
    checks = {}
    for key in limits:
        last_delta = d23["ph"][key] + d23["ah"][key]
        prior_delta = d12["ph"][key] + d12["ah"][key]
        uncertainty = max(floor[key], 3.0 * (last_delta + resampling[key]))
        checks[key] = {
            "residual": actual[key], "physical_limit": limits[key],
            "delta12_ph": d12["ph"][key], "delta12_ah": d12["ah"][key],
            "delta23_ph": d23["ph"][key], "delta23_ah": d23["ah"][key],
            "resampling": resampling[key], "three_sigma_budget": uncertainty,
            "physical_pass": actual[key] <= limits[key],
            "explained_by_numeric": actual[key] <= uncertainty,
            "self_convergence_pass": last_delta <= 0.2 * limits[key],
            "monotonic_or_roundoff": (last_delta <= prior_delta + floor[key] or
                                      last_delta <= floor[key]),
        }
    passed = all(all(check[name] for name in
                     ("physical_pass", "explained_by_numeric",
                      "self_convergence_pass", "monotonic_or_roundoff"))
                 for check in checks.values())
    return {"status": "pass" if passed else "needs_extra_levels",
            "checks": checks, "sample_step_km": [0.5, 0.25]}


def stage_a2_summary(run: Path, *, dataset: str, interpolation: str,
                     requested_k2: float | None, tag: str,
                     subset: str, extra_id: str | None,
                     fd_base_km: float) -> None:
    a2 = run / "a2"
    root = (a2 / "extra" / extra_id if dataset == "extra" and extra_id else
            a2 / dataset if dataset != "base" else a2)
    if dataset == "extra" and not extra_id:
        raise ValueError("extra dataset requires --a2-extra-id")
    manifest = read_json(root / "input_manifest.json")
    coefficient = (float(requested_k2) if requested_k2 is not None else
                   float(manifest["identified_K2"]) if dataset != "base" else
                   K_FP ** 2)
    label = "common" if np.isclose(coefficient, K_FP ** 2, rtol=0.0, atol=1e-12) else "identified"
    ah_dir = root / "ah" / tag / label / interpolation
    if fd_base_km != 0.5:
        ah_dir = ah_dir / f"fd{fd_base_km:g}"
    jobs = [job for job in manifest["jobs"]
            if subset == "all" or job["stage"] == subset]
    rows = [_a2_compare(job, root, ah_dir, interpolation) for job in jobs]
    groups = {}
    for row in rows:
        key = f"{row['case_id']}|f{row['frequency_MHz']:g}|e{row['elevation_deg']:.14g}"
        groups.setdefault(key, []).append(row)
    jobs_by_id = {job["job_id"]: job for job in jobs}
    gates = []
    for key, levels in groups.items():
        levels.sort(key=lambda item: item["level"])
        third = levels[-1]
        numeric = _a2_numeric_gate(levels, jobs_by_id, root, ah_dir)
        gates.append({"case": key, "gate": numeric["status"],
                      "numeric": numeric,
                      "levels": [item["level"] for item in levels],
                      "rms_m": [item.get("position_rms_m") for item in levels],
                      "max_m": [item.get("position_max_m") for item in levels],
                      "group_rel": [item.get("group_relative_error") for item in levels]})
    analysis = root / "analysis"
    analysis.mkdir(exist_ok=True)
    fd_suffix = "" if fd_base_km == 0.5 else f"_fd{str(fd_base_km).replace('.', 'p')}"
    path = analysis / f"a2_{dataset}_{tag}_{label}_{interpolation}_{subset}{fd_suffix}.json"
    write_json(path, {
        "schema": "a4-c1m-a2-analysis-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": dataset, "subset": subset, "tag": tag,
        "coefficient_K2": coefficient, "interpolation": interpolation,
        "spatial_fd_base_km": fd_base_km,
        "input_manifest_sha256": sha(root / "input_manifest.json"),
        "coefficient_freeze_sha256": sha(a2 / "identified" / "coefficient_freeze.json")
            if (a2 / "identified" / "coefficient_freeze.json").exists() else None,
        "limits": LIMIT, "rows": rows, "gates": gates,
        "completed_ah_count": sum(row["ah_status"] == "completed" for row in rows),
        "missing_ah_count": sum(row["ah_status"] == "missing" for row in rows),
    })
    print(f"A2 summary: {path} ({len(rows)} jobs)", flush=True)


def stage_a2_artifacts(run: Path) -> None:
    a2 = run / "a2"
    analysis = a2 / "analysis"
    named = {
        "common_linear": analysis / "a2_base_formal_common_linear_all.json",
        "common_hermite": analysis / "a2_base_formal_common_hermite_smooth.json",
        "internal_common_hermite": analysis / "a2_base_formal_common_hermite_internal.json",
        "identified_linear": analysis / "a2_base_formal_identified_linear_smooth.json",
        "identified_hermite": analysis / "a2_base_formal_identified_hermite_all.json",
        "identified_hermite_fd025": analysis / "a2_base_formal_identified_hermite_smooth_fd0p25.json",
        "identified_hermite_fd0125": analysis / "a2_base_formal_identified_hermite_smooth_fd0p125.json",
        "identified_hermite_fd1": analysis / "a2_base_formal_identified_hermite_smooth_fd1p0.json",
        "internal_identified": a2 / "identified" / "analysis" /
            "a2_identified_formal_identified_hermite_internal.json",
        "holdout_identified": a2 / "identified_holdout" / "analysis" /
            "a2_identified_holdout_formal_identified_hermite_affine.json",
    }
    missing = [name for name, path in named.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing A2 summary: {missing}")
    summaries = {name: read_json(path) for name, path in named.items()}
    manifests = {}
    for dataset in ("base", "identified", "identified_holdout"):
        folder = a2 if dataset == "base" else a2 / dataset
        manifests[dataset] = {
            job["job_id"]: job for job in read_json(folder / "input_manifest.json")["jobs"]}
    case_csv = analysis / "a2_per_case.csv"
    point_csv = analysis / "a2_per_point.csv"
    case_fields = ["candidate", "dataset", "job_id", "case_id", "stage",
                   "frequency_MHz", "elevation_deg", "level", "gate",
                   "ah_status", "ah_failure_reason", "comparison_arc_km",
                   "group_comparison_kind", "position_rms_m", "position_max_m",
                   "position_endpoint_m", "group_relative_error",
                   "ph_Ne_interior_rms_m3", "ph_Ne_first_interval_rms_m3",
                   "ph_Ne_below_grid_max_m3", "initial_n2_difference",
                   "same_initial_state", "ah_geom_path_km", "ph_geom_path_km",
                   "ah_group_path_km_native", "ph_group_for_gate_km"]
    point_fields = ["candidate", "dataset", "job_id", "s_km",
                    "ph_x_km", "ph_y_km", "ph_z_km",
                    "ah_x_km", "ah_y_km", "ah_z_km", "position_error_m"]
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
    colors = {"common_linear": "#8c564b", "common_hermite": "#ff7f0e",
              "identified_linear": "#1f77b4", "identified_hermite": "#2ca02c",
              "identified_hermite_fd0125": "#9467bd",
              "identified_hermite_fd1": "#7f7f7f"}
    with case_csv.open("w", encoding="utf-8", newline="") as case_stream, \
            point_csv.open("w", encoding="utf-8", newline="") as point_stream:
        case_writer = csv.DictWriter(case_stream, fieldnames=case_fields)
        point_writer = csv.writer(point_stream)
        case_writer.writeheader()
        point_writer.writerow(point_fields)
        for candidate, summary in summaries.items():
            dataset = summary["dataset"]
            dataset_dir = a2 if dataset == "base" else a2 / dataset
            label = ("common" if np.isclose(summary["coefficient_K2"], K_FP ** 2,
                                            rtol=0.0, atol=1e-12) else "identified")
            ah_dir = (dataset_dir / "ah" / summary["tag"] / label /
                      summary["interpolation"])
            fd = float(summary.get("spatial_fd_base_km", 0.5))
            if fd != 0.5:
                ah_dir = ah_dir / f"fd{fd:g}"
            gate_by_case = {item["case"]: item["gate"] for item in summary["gates"]}
            for row in summary["rows"]:
                key = (f"{row['case_id']}|f{row['frequency_MHz']:g}|"
                       f"e{row['elevation_deg']:.14g}")
                case_writer.writerow({field: (candidate if field == "candidate" else
                                               dataset if field == "dataset" else
                                               gate_by_case.get(key) if field == "gate" else
                                               row.get(field)) for field in case_fields})
                if row["ah_status"] != "completed":
                    continue
                job = manifests[dataset][row["job_id"]]
                ph = _a2_ph_trace(dataset_dir / "pharlap_raw" /
                                  f"{row['job_id']}.mat", job)
                with np.load(ah_dir / f"{row['job_id']}.npz") as raw:
                    ah_s = np.array(raw["geom_path_km"])
                    ah_xyz = np.array(raw["xyz_km"])
                end = float(row["comparison_arc_km"])
                fixed = bool(row["comparison_fixed_100km"])
                samples = np.linspace(0.0, end, 201 if fixed else 501)
                ph_xyz = resample(ph["s"], ph["xyz"], samples)
                ah_sample = resample(ah_s, ah_xyz, samples)
                errors = np.linalg.norm(ph_xyz - ah_sample, axis=1) * 1000.0
                for s, p, q, error in zip(samples, ph_xyz, ah_sample, errors):
                    point_writer.writerow([candidate, dataset, row["job_id"], s,
                                           *p, *q, error])
                if (row["stage"] == "smooth" and row["level"] == 3
                        and candidate in colors):
                    angles = (21.94290326504532, 21.94562872092976,
                              39.54254824570964, 39.58701138438118)
                    index = int(np.argmin(np.abs(np.array(angles) - row["elevation_deg"])))
                    axes.flat[index].plot(samples, np.maximum(errors, 1e-7),
                                          label=candidate, color=colors[candidate],
                                          linewidth=1.15)
    for index, ax in enumerate(axes.flat):
        angles = (21.94290326504532, 21.94562872092976,
                  39.54254824570964, 39.58701138438118)
        ax.set_title(f"elevation {angles[index]:.6f} deg")
        ax.set_yscale("log")
        ax.set_xlabel("common geometric arc (km)")
        ax.set_ylabel("position difference (m)")
        ax.grid(True, alpha=0.25)
    axes.flat[0].legend(fontsize=7)
    figure.tight_layout()
    figure_path = analysis / "a2_smooth_difference.png"
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    smooth = {}
    for angle in (21.94290326504532, 21.94562872092976,
                  39.54254824570964, 39.58701138438118):
        smooth[f"{angle:.14g}"] = {}
        for candidate in ("common_linear", "common_hermite",
                          "identified_linear", "identified_hermite",
                          "identified_hermite_fd025", "identified_hermite_fd0125",
                          "identified_hermite_fd1"):
            rows = [row for row in summaries[candidate]["rows"]
                    if row["stage"] == "smooth" and row["level"] == 3
                    and row["elevation_deg"] == angle]
            if len(rows) != 1:
                raise ValueError(f"missing smooth level-3 row: {angle}, {candidate}")
            row = rows[0]
            smooth[f"{angle:.14g}"][candidate] = {
                "rms_m": row["position_rms_m"],
                "max_m": row["position_max_m"],
                "natural_landing_m": row["position_endpoint_m"],
                "natural_group_rel": row["group_relative_error"]}
    boundary_rows = [row for row in summaries["identified_hermite"]["rows"]
                     if row["stage"] == "boundary"]
    boundary = {}
    for hmin in (20, 40, 60, 80):
        rows = [row for row in boundary_rows
                if row["case_id"] == f"boundary_h{hmin:03d}"]
        boundary[str(hmin)] = {
            "job_count": len(rows),
            "native_completed": sum(row["ah_status"] == "completed" for row in rows),
            "native_rejected": sum(row["ah_status"] == "physics_domain_rejected"
                                   for row in rows),
            "max_PH_Ne_below_grid_m3": max(
                (row["ph_Ne_below_grid_max_m3"] or 0.0) for row in rows),
            "first_interval_PH_vs_Hermite_rms_m3": [
                row["ph_Ne_first_interval_rms_m3"] for row in rows
                if row["ph_Ne_first_interval_rms_m3"] is not None],
        }
    holdout = summaries["holdout_identified"]
    freeze = read_json(a2 / "identified" / "coefficient_freeze.json")
    import shutil
    runner_hash = sha(Path(__file__))
    runner_snapshot = a2 / "sources" / "analysis" / f"runner_{runner_hash[:16]}.py"
    runner_snapshot.parent.mkdir(parents=True, exist_ok=True)
    if not runner_snapshot.exists():
        shutil.copyfile(Path(__file__), runner_snapshot)
    if sha(runner_snapshot) != runner_hash:
        raise ValueError("analysis runner snapshot mismatch")
    matrix_path = analysis / "a2_cause_matrix.json"
    write_json(matrix_path, {
        "schema": "a4-c1m-a2-cause-matrix-v1",
        "analysis_runner_sha256": runner_hash,
        "analysis_runner_snapshot": str(runner_snapshot.resolve()),
        "identified_K2": freeze["identified_K2"],
        "coefficient_freeze_sha256": sha(a2 / "identified" / "coefficient_freeze.json"),
        "smooth_level3": smooth,
        "boundary": boundary,
        "matched_holdout": {
            "jobs": len(holdout["rows"]),
            "gate_counts": {status: sum(g["gate"] == status for g in holdout["gates"])
                            for status in sorted(set(g["gate"] for g in holdout["gates"]))},
            "max_level3_rms_m": max(row["position_rms_m"] for row in holdout["rows"]
                                    if row["level"] == 3),
            "max_level3_group_rel": max(row["group_relative_error"]
                                        for row in holdout["rows"] if row["level"] == 3),
        },
        "summary_sha256": {name: sha(path) for name, path in named.items()},
        "prepare_replay_sha256": sha(analysis / "a2_prepare_replay.json")
            if (analysis / "a2_prepare_replay.json").exists() else None,
        "artifacts": {"per_case_csv_sha256": sha(case_csv),
                      "per_point_csv_sha256": sha(point_csv),
                      "difference_png_sha256": sha(figure_path)},
    })
    print(f"A2 artifacts: {matrix_path}, {case_csv}, {point_csv}, {figure_path}", flush=True)


def numeric_decision_v2(actual: dict, d12: dict, d23: dict,
                        resampling: dict, *, eligible: bool = True) -> dict:
    """New-run gate: each solver must converge; method differences are not bounds.

    resampling contains independently established errors, not the difference
    between two interpolation choices. Historical A2 gates remain unchanged.
    """
    limits = {**LIMIT, "landing_m": LIMIT["max_m"]}
    checks = {}
    for key, limit in limits.items():
        values = [actual[key], resampling[key],
                  *(d[key] for level in (d12, d23) for d in level.values())]
        if any(not np.isfinite(v) or v < 0 for v in values):
            raise ValueError("non-finite or negative numerical evidence")
        sides = {side: {
            "monotonic": bool(d23[side][key] <= d12[side][key]),
            "last_change_small": bool(d23[side][key] <= 0.2 * limit),
        } for side in ("ph", "ah")}
        budget = 3.0 * (d23["ph"][key] + d23["ah"][key] + resampling[key])
        checks[key] = {"physical_pass": bool(actual[key] <= limit),
                       "explained_by_numeric": bool(actual[key] <= budget),
                       "three_times_budget": float(budget), "sides": sides}
    passed = eligible and all(
        row["physical_pass"] and row["explained_by_numeric"]
        and all(all(side.values()) for side in row["sides"].values())
        for row in checks.values())
    return {"passed": bool(passed), "eligible": bool(eligible), "checks": checks,
            "roundoff_waiver_applied": False,
            "interpolation_method_difference_in_budget": False}


def _p0_trace_contract(trace: dict, scope: str) -> list[str]:
    errors = []
    expected_n2 = float(trace["expected_initial_n2"])
    if (not np.isfinite(expected_n2) or expected_n2 <= 0 or
            not np.all(np.isfinite(trace["expected_origin_km"])) or
            not np.all(np.isfinite(trace["expected_initial_p"]))):
        errors.append("invalid_expected_state")
    if trace["units"] != {"xyz": "km", "s": "km", "group": "km", "p": "dimensionless"}:
        errors.append("unit_contract")
    arrays = {key: np.asarray(trace[key], dtype=float) for key in ("s", "xyz", "group", "p")}
    s, xyz, group, p = (arrays[key] for key in ("s", "xyz", "group", "p"))
    if (s.ndim != 1 or s.size < 2 or xyz.shape != (s.size, 3)
            or p.shape != xyz.shape or group.shape != s.shape
            or not all(np.all(np.isfinite(v)) for v in arrays.values())):
        return errors + ["invalid_raw_arrays"]
    if (abs(s[0]) > 1e-10 or abs(group[0]) > 1e-10 or group[-1] <= 0 or
            np.any(np.diff(s) <= 0) or np.any(np.diff(group) < 0)):
        errors.append("invalid_path_coordinate")
    if not np.allclose(xyz[0], trace["expected_origin_km"], rtol=0, atol=1e-8):
        errors.append("origin_mismatch")
    if not np.allclose(p[0], trace["expected_initial_p"], rtol=0, atol=1e-10):
        errors.append("initial_wavevector_mismatch")
    if abs(float(p[0] @ p[0]) - trace["expected_initial_n2"]) > 1e-10:
        errors.append("initial_dispersion_mismatch")
    chords = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    if np.any(chords > np.diff(s) * 1.001 + 1e-8):
        errors.append("position_arc_scale_mismatch")
    if trace["termination"] not in {"path_complete", "ground", "domain_exit", "penetrated"}:
        errors.append("invalid_termination")
    if scope == "local":
        if s[-1] < S_KM - 1e-9:
            errors.append("insufficient_local_coverage")
    elif scope == "full":
        if trace["termination"] != "ground":
            errors.append("full_path_not_ground")
        if abs(float(np.linalg.norm(xyz[-1])) - R_KM) > 5e-5:
            errors.append("ground_position_mismatch")
    else:
        raise ValueError("scope must be local or full")
    return errors


def _p0_pair_metrics(a: dict, b: dict, scope: str = "local", step: float = 0.25) -> dict:
    end = S_KM if scope == "local" else min(a["s"][-1], b["s"][-1])
    samples = np.linspace(0.0, end, int(np.ceil(end / step)) + 1)
    result = norm_metrics(resample(a["s"], a["xyz"], samples)
                          - resample(b["s"], b["xyz"], samples))
    endpoints = [resample(t["s"], t["xyz"], np.array([end]))[0]
                 if scope == "local" else np.asarray(t["xyz"])[-1] for t in (a, b)]
    groups = [float(resample(t["s"], t["group"], np.array([end]))[0])
              if scope == "local" else float(t["group"][-1]) for t in (a, b)]
    return {"rms_m": result["rms_m"], "max_m": result["max_m"],
            "landing_m": float(np.linalg.norm(endpoints[0] - endpoints[1]) * 1000),
            "group_rel": abs(groups[0] - groups[1]) / max(abs(groups[0]), 1e-15)}


def _p0_fixture(curved: bool = False, nodes: int = 41, offset_m: float = 0.0) -> dict:
    s = np.linspace(0.0, S_KM, nodes)
    if curved:
        xyz = np.column_stack((R_KM * np.sin(s / R_KM), np.zeros_like(s),
                               R_KM * np.cos(s / R_KM)))
        p = np.column_stack((np.cos(s / R_KM), np.zeros_like(s), -np.sin(s / R_KM)))
    else:
        xyz = np.column_stack((s, np.zeros_like(s), np.full_like(s, R_KM)))
        p = np.tile([1.0, 0.0, 0.0], (nodes, 1))
    xyz[:, 1] += offset_m / 1000.0 * (s / S_KM) ** 2
    return {"s": s, "xyz": xyz, "group": s.copy(), "p": p,
            "units": {"xyz": "km", "s": "km", "group": "km", "p": "dimensionless"},
            "expected_origin_km": [0.0, 0.0, R_KM], "expected_initial_p": [1.0, 0.0, 0.0],
            "expected_initial_n2": 1.0, "termination": "path_complete"}


def p0_checks() -> dict:
    """Deterministic analytic fixtures; this does not execute either ray tracer."""
    checks, evidence = {}, {}
    zero = {key: 0.0 for key in (*LIMIT, "landing_m")}
    for name, curved in (("straight", False), ("curved", True)):
        trace = _p0_fixture(curved)
        metrics = _p0_pair_metrics(trace, trace)
        decision = numeric_decision_v2(metrics, {"ph": zero, "ah": zero},
                                       {"ph": zero, "ah": zero}, zero)
        checks[name + "_identity_passes"] = not _p0_trace_contract(trace, "local") and decision["passed"]
        evidence[name] = {"metrics": metrics, "decision": decision}
    variants = {}
    wrong = _p0_fixture(); wrong["group"] -= 1000.0
    variants["negative_group_origin_rejected"] = (wrong, "local", "invalid_path_coordinate")
    wrong = _p0_fixture(); wrong["expected_initial_n2"] = float("nan")
    variants["invalid_expected_dispersion_rejected"] = (wrong, "local", "invalid_expected_state")
    wrong = _p0_fixture(); wrong["units"]["xyz"] = "m"
    variants["unit_contract_rejected"] = (wrong, "local", "unit_contract")
    wrong = _p0_fixture(); wrong["xyz"] *= 1000.0
    variants["mislabeled_position_scale_rejected"] = (wrong, "local", "origin_mismatch")
    wrong = _p0_fixture(); wrong["p"][0] *= -1.0
    variants["reversed_initial_direction_rejected"] = (wrong, "local", "initial_wavevector_mismatch")
    wrong = _p0_fixture(); wrong["p"][0] *= 1000.0
    variants["wavevector_scale_rejected"] = (wrong, "local", "initial_dispersion_mismatch")
    wrong = _p0_fixture()
    for key in ("s", "xyz", "group", "p"):
        wrong[key] = wrong[key][:-1]
    variants["early_local_exit_rejected"] = (wrong, "local", "insufficient_local_coverage")
    wrong = _p0_fixture(True); wrong["termination"] = "domain_exit"
    variants["ph_full_early_exit_rejected"] = (wrong, "full", "full_path_not_ground")
    wrong = _p0_fixture(); wrong["termination"] = "ground"
    variants["false_ground_position_rejected"] = (wrong, "full", "ground_position_mismatch")
    for name, (trace, scope, error) in variants.items():
        errors = _p0_trace_contract(trace, scope)
        checks[name] = error in errors
        evidence[name] = {"errors": errors, "expected_error": error}
    ground = _p0_fixture(True); ground["termination"] = "ground"
    checks["valid_full_ground_passes"] = not _p0_trace_contract(ground, "full")
    # A PH ray can exit its domain after the comparison window without invalidating that prefix.
    later = _p0_fixture(); later["termination"] = "domain_exit"
    checks["valid_local_domain_endpoint_passes"] = not _p0_trace_contract(later, "local")
    ph = [_p0_fixture(offset_m=x) for x in (1.05, 0.05, 0.0)]
    ah = [_p0_fixture(offset_m=x) for x in (0.25, 0.30, 0.50)]
    d12 = {side: _p0_pair_metrics(traces[0], traces[1]) for side, traces in (("ph", ph), ("ah", ah))}
    d23 = {side: _p0_pair_metrics(traces[1], traces[2]) for side, traces in (("ph", ph), ("ah", ah))}
    result = numeric_decision_v2(_p0_pair_metrics(ph[-1], ah[-1]), d12, d23, zero)
    checks["single_solver_nonmonotonic_rejected"] = not result["passed"] and not result["checks"]["rms_m"]["sides"]["ah"]["monotonic"]
    evidence["single_solver_nonmonotonic"] = {"d12": d12, "d23": d23, "decision": result}
    # A large difference between interpolation methods is NOT a measured error bound.
    coarse, fine = _p0_fixture(nodes=3), _p0_fixture(nodes=5)
    for trace in (coarse, fine):
        trace["xyz"][:, 1] = 0.01 * np.sin(np.pi * trace["s"] / S_KM)
    residual = _p0_pair_metrics(coarse, fine)
    result = numeric_decision_v2(residual, {"ph": zero, "ah": zero}, {"ph": zero, "ah": zero}, zero)
    checks["method_difference_cannot_rescue_plateau"] = not result["passed"]
    evidence["method_difference_plateau"] = {"metrics": residual, "decision": result}
    # Check actual interpolation accuracy against the analytic circle, with nested raw nodes.
    errors = []
    sample = np.linspace(0, S_KM, 1601)
    exact = np.column_stack((R_KM*np.sin(sample/R_KM), np.zeros_like(sample), R_KM*np.cos(sample/R_KM)))
    for nodes in (11, 21, 41):
        trace = _p0_fixture(True, nodes)
        errors.append(norm_metrics(resample(trace["s"], trace["xyz"], sample) - exact))
    checks["curved_raw_refinement_converges"] = all(
        errors[i + 1][key] < errors[i][key] for i in (0, 1) for key in ("rms_m", "max_m"))
    evidence["curved_analytic_refinement"] = errors
    full_a, full_b = copy.deepcopy(ground), copy.deepcopy(ground)
    full_b["group"] = np.asarray(full_b["group"]) * 1.0001
    checks["native_total_group_difference_detected"] = _p0_pair_metrics(full_a, full_b, "full")["group_rel"] > LIMIT["group_rel"]
    return {"schema": "c1m-p0-checks-v2", "checks": {k: bool(v) for k, v in checks.items()},
            "passed": bool(all(checks.values())), "evidence": evidence}


def p0_a2_audit(parent: Path) -> dict:
    """Read-only audit of the previously accepted smooth/holdout candidates."""
    datasets = [
        (parent, "a2_base_formal_identified_hermite_smooth_fd0p125.json",
         parent / "ah/formal/identified/hermite/fd0.125", True),
        (parent / "identified_holdout", "a2_identified_holdout_formal_identified_hermite_affine.json",
         parent / "identified_holdout/ah/formal/identified/hermite", False)]
    violations, raw_checks, sources = [], [], {}
    group_count = 0
    for folder, name, ah_folder, full in datasets:
        summary_path = folder / "analysis" / name
        summary = read_json(summary_path)
        sources[str(summary_path)] = sha(summary_path)
        for gate in summary["gates"]:
            group_count += 1
            for metric, values in gate["numeric"]["checks"].items():
                for side in ("ph", "ah"):
                    first, last = values["delta12_" + side], values["delta23_" + side]
                    if last > first:
                        violations.append({"case": gate["case"], "side": side, "metric": metric,
                                           "delta12": first, "delta23": last, "full_path": full})
        for row in summary["rows"]:
            job_id = row["job_id"]
            ph_path = folder / "pharlap_raw" / (job_id + ".mat")
            ah_path = ah_folder / (job_id + ".npz")
            ah_record_path = ah_folder / (job_id + ".json")
            sources[str(ah_record_path)] = sha(ah_record_path)
            if sha(ph_path) != row["ph_raw_sha256"] or sha(ah_path) != row["ah_raw_npz_sha256"]:
                raise ValueError("historical raw hash mismatch")
            sources[str(ph_path)] = row["ph_raw_sha256"]
            sources[str(ah_path)] = row["ah_raw_npz_sha256"]
            raw = mat(ph_path)
            state = raw["ray_state_vec"]
            job = raw["job"]
            grid_path = Path(job["grid_mat"])
            if sha(grid_path) != job["grid_sha256"]:
                raise ValueError("historical initial-state grid hash mismatch")
            sources[str(grid_path)] = sha(grid_path)
            grid = mat(grid_path)
            ne0 = float(np.interp(float(job["origin_height_km"]),
                np.asarray(grid["height_km"]).reshape(-1),
                np.asarray(grid["ne_profile_m3"]).reshape(-1)))
            expected_n2 = 1 - float(job["launch_k2"]) * ne0 / (float(job["frequency_MHz"]) * 1e6)**2
            expected_p = np.sqrt(expected_n2) * launch_direction_ecef(
                30.0, 120.0, float(job["elevation_deg"]), 0.0, R_KM)
            expected_xyz = np.asarray(geocentric_to_ecef(
                30.0, 120.0, float(job["origin_height_km"]), R_KM))
            pp = np.array([vector(state, key)[0] for key in ("dir_x", "dir_y", "dir_z")])
            ph_xyz = np.column_stack([vector(state, key) for key in ("pos_x", "pos_y", "pos_z")]) / 1000
            ph_s = vector(state, "geometrical_path") / 1000
            with np.load(ah_path) as ah:
                ah_origin = np.asarray(ah["xyz_km"][0])
                pa = ah["wavevector_km_inv"][0] * C_KM_S / (2*np.pi*row["frequency_MHz"]*1e6)
                origin_error = float(np.linalg.norm(ah["xyz_km"][0] - ph_xyz[0]) * 1000)
                dispersion_error = float(abs(pa @ pa - ah["n2_native"][0]))
                ah_s = ah["geom_path_km"]
                covered = bool(ph_s[-1] >= 100 - 1e-9 and ah_s[-1] >= 100 - 1e-9)
                ph_end_h = float(np.linalg.norm(ph_xyz[-1]) - R_KM)
                ah_end_h = float(np.linalg.norm(ah["xyz_km"][-1]) - R_KM)
            label = int(scalar(raw["ray_data"], "ray_label"))
            path = raw["ray_path_data"]
            ph_event_h = float(vector(path, "height")[-1])
            ph_event_xyz = np.asarray(geocentric_to_ecef(
                float(vector(path, "lat")[-1]), float(vector(path, "lon")[-1]),
                ph_event_h, R_KM), dtype=float)
            event_state_distance_m = float(np.linalg.norm(ph_event_xyz - ph_xyz[-1]) * 1000)
            terminal_ok = (label == 1 and row["ah_termination_status"] == "ground"
                           and abs(ph_end_h) < 5e-5 and abs(ah_end_h) < 5e-5) if full else (
                               covered and label in (1, -2, -3)
                               and row["ah_termination_status"] == "fixed_geometric_path_completed")
            vector_error = float(np.linalg.norm(pp - pa))
            launch_error = float(max(np.linalg.norm(pp - expected_p), np.linalg.norm(pa - expected_p)))
            position_error = float(max(np.linalg.norm(ph_xyz[0] - expected_xyz),
                np.linalg.norm(ah_origin - expected_xyz)) * 1000)
            raw_checks.append({"job_id": job_id, "full_path": full,
                               "initial_p_vector_error": vector_error,
                "prescribed_launch_p_error": launch_error,
                "prescribed_origin_error_m": position_error,
                               "initial_position_error_m": origin_error,
                               "ah_initial_dispersion_error": dispersion_error,
                "ph_ray_label": label, "terminal_or_local_coverage_ok": bool(terminal_ok),
                "ah_termination_status": row["ah_termination_status"],
                "ph_state_end_height_m": ph_end_h * 1000,
                "ph_event_end_height_m": ph_event_h * 1000,
                "ah_end_height_m": ah_end_h * 1000,
                "ph_event_state_distance_m": event_state_distance_m,
                "reported_ground_event_ok": bool(label == 1 and abs(ph_event_h) < 5e-5),
                "state_endpoint_ground_ok": bool(abs(ph_end_h) < 5e-5),
                               "passed": bool(vector_error <= 1e-10 and origin_error <= 1e-5
                                             and dispersion_error <= 1e-10 and launch_error <= 1e-10
                                             and position_error <= 1e-5 and terminal_ok)})
    return {"schema": "c1m-p0-parent-audit-v1", "source_sha256": sources,
            "groups_checked": group_count, "single_solver_metrics_checked": group_count * 8,
            "strict_nonmonotonic": violations,
            "strict_nonmonotonic_groups": len({v["case"] for v in violations}),
            "raw_contract_checks": raw_checks, "raw_contract_passed": bool(all(r["passed"] for r in raw_checks)),
            "historical_records_modified": False, "new_ray_computation": False,
            "interpretation": "No automatic roundoff waiver. Nanometre-scale nonmonotonic changes need separate platform evidence; this is not a new physical-error attribution."}


def stage_p0(run: Path) -> None:
    """Freeze first, then create and check synthetic inputs. Never overwrites a run."""
    parent = OUT / "a2_20260928_v1" / "a2"
    frozen = read_json(parent / "sha256_manifest.json")
    assets = source_files()
    assets["utils"] = ROOT / "utils.py"
    assets["a2_report"] = ROOT / "word" / "a4_c1m_a2_attribution_report.txt"
    assets["review_plan"] = ROOT / "word" / "a4_c1m_next_attribution_and_repair_plan_20261007.txt"
    assets["a2_completion"] = parent / "completion.json"
    assets["a2_readonly"] = parent / "diagnostics" / "a2_readonly_verification.json"
    for record in frozen["records"]:
        if record["path"] in {str(assets["a2_completion"]), str(assets["a2_readonly"])}:
            if sha(Path(record["path"])) != record["sha256"]:
                raise ValueError("parent A2 evidence changed")
    if any(not path.is_file() for path in assets.values()):
        raise FileNotFoundError("P0 source or parent evidence missing")
    run.mkdir(parents=True, exist_ok=False)
    sources = run / "sources"; sources.mkdir()
    records = {}
    for key, path in assets.items():
        digest = sha(path)
        snapshot = None
        if path.suffix.lower() not in {".exe", ".mexw64", ".mat"}:
            target = sources / (key + path.suffix)
            shutil.copyfile(path, target)
            if sha(target) != digest:
                raise ValueError("source changed during freezing")
            snapshot = str(target.relative_to(run))
        records[key] = {"path": str(path), "sha256": digest, "snapshot": snapshot}
    version = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         f"(Get-Item -LiteralPath '{MATLAB_EXE}').VersionInfo.ProductVersion"],
        capture_output=True, text=True, check=True,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)).stdout.strip()
    # This manifest is saved BEFORE any fixture generation or comparison.
    write_json(run / "p0_manifest.json", {
        "schema": "c1m-p0-freeze-v1", "run_id": run.name,
        "created_utc": datetime.now(timezone.utc).isoformat(), "sources": records,
        "parent_a2": str(parent), "preparation_snapshot_before_inputs": True,
        "python_version": sys.version, "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "matlab_executable_product_version": version,
        "ph_tolerance_levels": TOL, "ah_integration_levels": AH_PROFILES,
        "production_modified": False, "pharlap_executed": False})
    fixtures = {"straight": _p0_fixture(), "curved": _p0_fixture(True)}
    write_json(run / "p0_inputs.json", fixtures)
    result = p0_checks()
    result.update({"input_sha256": sha(run / "p0_inputs.json"),
                   "evidence_integrity_passed": True,
                   "synthetic_numeric_checks_passed": result["passed"],
                   "physical_rays_validated_in_p0": False,
                   "parent_a2_global_passed": read_json(assets["a2_readonly"])["global_passed"]})
    write_json(run / "p0_checks.json", result)
    write_json(run / "p0_parent_audit.json", p0_a2_audit(parent))
    print(f"P0 synthetic checks: {sum(result['checks'].values())}/{len(result['checks'])}")
    if not result["passed"]:
        raise RuntimeError("P0 fixture check failed; frozen failed run retained")


P1_C = 80.616384943129


def _p1_job_hash(job: dict) -> dict:
    job = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in job.items()
           if k != "job_sha256"}
    job["job_sha256"] = hashlib.sha256(json.dumps(job, sort_keys=True,
        separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
    return job


def _p1_save_jobs(folder: Path, jobs: list[dict]) -> None:
    payload = np.empty((1, len(jobs)), dtype=object)
    for i, job in enumerate(jobs):
        payload[0, i] = job
    savemat(folder / "input.mat", {"a2_jobs": payload}, do_compression=True)
    write_json(folder / "jobs.json", {"jobs": jobs, "input_sha256": sha(folder / "input.mat")})


def stage_p1_prepare(run: Path) -> None:
    """Freeze P1 source before constructing any new medium; no production edits."""
    run.mkdir(parents=True, exist_ok=False)
    sources = run / "sources"; sources.mkdir()
    assets = source_files()
    assets.update({"utils": ROOT / "utils.py",
        "p0_report": ROOT / "word/a4_c1m_p0_implementation_report_20261007.txt",
        "p0_verification": OUT / "p0_20261007_v2/verification.json",
        "p0_manifest": OUT / "p0_20261007_v2/p0_manifest.json",
        "p0_delivery_sha256": OUT / "p0_20261007_v2/delivery_sha256.json",
                   "a2_report": ROOT / "word/a4_c1m_a2_attribution_report.txt",
                   "pharlap_api_doc": Path(r"D:\pharlap\pharlap_4.7.4\src\matlab\raytrace_3d_sp.m")})
    records = {}
    for key, path in assets.items():
        digest = sha(path)
        snapshot = None
        if path.suffix.lower() not in {".exe", ".mexw64", ".mat"}:
            target = sources / (key + path.suffix); shutil.copyfile(path, target)
            assert sha(target) == digest
            snapshot = str(target.relative_to(run))
        records[key] = {"path": str(path), "sha256": digest, "snapshot": snapshot}
    write_json(run / "p1_manifest.json", {"schema": "c1m-p1-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(), "sources": records,
        "python_executable": sys.executable, "python_version": sys.version,
        "numpy_version": np.__version__, "scipy_version": __import__("scipy").__version__,
        "coefficient": P1_C, "coefficient_source": "A2 independent affine identification",
        "ph_profiles": TOL, "ah_profiles": AH_PROFILES,
        "production_modified": False, "parent": "p0_20261007_v2",
        "scope": "unmagnetized radial entrance and output-state attribution"})
    grids = run / "grids"; grids.mkdir()
    jobs, grid_records = [], []
    lat = np.arange(28., 65.25, .5); lon = np.arange(118., 122.25, .5)
    bh = np.arange(0., 470., 10.); bshape = (len(lat), len(lon), len(bh))

    def add(case, kind, hb, nb, gradient=0., width=0., dh=2., frequencies=(10.,),
            full=False):
        lower = 0. if full else hb
        upper = 460. if full else lower + 100.
        h = np.linspace(lower, upper, int(round((upper-lower)/dh))+1)
        assert len(h) <= 401
        if full:
            h = np.arange(60., 461., 2.)
            ne = _a2_profile("internal", h)
        elif kind == "smooth":
            u = np.clip((h-hb)/width, 0., 1.)
            ne = nb * u**3 * (10.-15.*u+6.*u*u)
        else:
            ne = nb + gradient*(h-hb)
        shape = (len(lat), len(lon), len(h))
        path = grids / (case + ".mat")
        grid = {"latitude_deg": lat, "longitude_deg": lon, "height_km": h,
            "ne_profile_m3": ne, "Ne_m3": np.broadcast_to(ne, shape).copy(),
            "geomag_latitude_deg": lat, "geomag_longitude_deg": lon,
            "geomag_height_km": bh, "Bx_T": np.zeros(bshape),
            "By_T": np.zeros(bshape), "Bz_T": np.zeros(bshape),
            "collision_freq_Hz": np.zeros(shape), "earth_radius_m": R_KM*1000.}
        savemat(path, grid, do_compression=True)
        digest = sha(path)
        grid_records.append({"case_id": case, "path": str(path), "sha256": digest,
            "kind": kind, "h_boundary_km": hb, "Nb_m3": nb, "g_m3_per_km": gradient,
            "width_km": width, "dh_km": dh,
            "array_sha256": {k: array_sha(np.asarray(v)) for k, v in grid.items()}})
        angles = (21.94290326504532, 21.94562872092976, 39.54254824570964, 39.58701138438118) if full else (22., 39.5)
        for f in frequencies:
            for angle in angles:
                for level in (1, 2, 3):
                    job = _a2_job(case, "boundary" if full else "p1_local", path, digest,
                        f, angle, 0. if full else hb-10., level, explicit=False, launch_k2=P1_C)
                    job.update({"p1_case": kind, "h_boundary_km": hb, "Nb_m3": nb,
                        "g_m3_per_km": gradient, "width_km": width, "dh_km": dh,
                        "scope": "full" if full else "local", "role": "training" if f==10. else "holdout",
                        "exporter_sha256": records["matlab_source"]["sha256"]})
                    jobs.append(_p1_job_hash(job))

    for nb in (0., 1e8, 5e8, 1e9):
        add(f"uniform_n{nb:g}_h60_d2", "uniform", 60., nb)
    add("gradient_h60_d2", "gradient", 60., 5e8, gradient=1e8)
    for hb in (20., 80.):
        add(f"uniform_n5e8_h{hb:g}_d2", "uniform", hb, 5e8)
    for width in (2., 4., 8.):
        add(f"smooth_w{width:g}_d2", "smooth", 80.-width, 5e8, width=width)
    for dh in (1., .5, .25):
        fs = (10., 8., 15.) if dh == .25 else (10.,)
        add(f"uniform_n5e8_h60_d{dh:g}", "uniform", 60., 5e8, dh=dh, frequencies=fs)
        add(f"gradient_h60_d{dh:g}", "gradient", 60., 5e8, gradient=1e8, dh=dh, frequencies=fs)
        for width in (2., 4., 8.):
            add(f"smooth_w{width:g}_d{dh:g}", "smooth", 80.-width, 5e8,
                width=width, dh=dh, frequencies=fs)
    # Density/position holdouts reuse identical grids; no fitted parameter changes.
    for original in list(jobs):
        if original["dh_km"] == 2 and original["p1_case"] == "uniform":
            for f in (8., 15.):
                job = dict(original); job["frequency_MHz"] = f; job["role"] = "holdout"
                job["job_id"] = job["job_id"].replace("_f10_", f"_f{f:g}_")
                jobs.append(_p1_job_hash(job))
    add("c1_h60_d2", "chapman_hard", 60., float(_a2_profile("internal", np.array([60.]))[0]), full=True)
    # A/B/A/B explicitly cleared full-grid calls test deterministic order independence.
    for index, nb in enumerate((1e8, 1e9, 1e8, 1e9)):
        source = next(j for j in jobs if j["p1_case"]=="uniform" and j["Nb_m3"]==nb
                      and j["frequency_MHz"]==10. and j["elevation_deg"]==22. and j["level"]==3)
        job = dict(source); job["job_id"] = f"cache_{index+1}_"+job["job_id"]
        job["role"] = "cache_control"; jobs.append(_p1_job_hash(job))
    _p1_save_jobs(run, jobs)
    write_json(run / "grids.json", {"grids": grid_records})
    print(f"P1 prepared {len(jobs)} jobs, {len(grid_records)} grids", flush=True)


def _p1_ah_one(job: dict, output: str, hashes: dict, fd: float) -> dict:
    import models.ray_tracer_ah_3d_ecef as solver
    original = solver.launch_state_ah3d_ecef
    if "initial_p" in job:
        def restart(model, **kwargs):
            return np.r_[np.asarray(job["initial_xyz_km"]), model.k0*np.asarray(job["initial_p"]), 0.]
        solver.launch_state_ah3d_ecef = restart
    try:
        return _a2_ah_one(job, output, P1_C, "hermite", fd, hashes)
    finally:
        solver.launch_state_ah3d_ecef = original


def stage_p1_ah(run: Path, workers: int, subset: str, fd: float) -> None:
    if workers < 1 or not np.isfinite(fd) or fd <= 0:
        raise ValueError("P1 requires positive workers and finite positive FD scale")
    manifest = read_json(run / "p1_manifest.json")
    for key in ("runner", "ah_solver", "ah_dispersion", "ah_provider", "config"):
        rec = manifest["sources"][key]
        if sha(Path(rec["path"])) != rec["sha256"]:
            raise ValueError(f"P1 execution source changed: {key}")
    jobs = read_json(run / "jobs.json")["jobs"]
    if subset != "all":
        jobs = [j for j in jobs if j["p1_case"] == subset or j["role"] == subset]
    if not jobs:
        raise ValueError("empty P1 selection")
    output = run / "ah" / f"fd{fd:g}"; output.mkdir(parents=True, exist_ok=True)
    hashes = _a2_ah_source_hashes()
    pending = []
    for job in jobs:
        path = output / (job["job_id"] + ".json")
        if path.exists():
            old = read_json(path)
            if old["job_sha256"] != job["job_sha256"] or old["source_sha256"] != hashes:
                raise ValueError("P1 AH resume identity mismatch")
            raw = output / (job["job_id"] + ".npz")
            if raw.exists() and sha(raw) != old.get("raw_npz_sha256"):
                raise ValueError("P1 AH resume raw mismatch")
            continue
        pending.append(job)
    print(f"P1 AH selected={len(jobs)} pending={len(pending)} fd={fd}", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_p1_ah_one, j, str(output), hashes, fd) for j in pending]
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            print(f"P1 AH {i}/{len(pending)} {result['job_id']} {result['status']} {result['elapsed_s']:.2f}s", flush=True)


def stage_p1_restart_prepare(run: Path) -> None:
    """Use one interior native state per frequency/angle, shared across three levels."""
    folder = run / "restart"; folder.mkdir(exist_ok=False)
    shutil.copyfile(run / "p1_manifest.json", folder / "p1_manifest.json")
    jobs = []
    for parent in read_json(run / "jobs.json")["jobs"]:
        if not (parent["p1_case"] == "gradient" and parent["dh_km"] == .25 and parent["level"] == 3):
            continue
        rawpath = run / "pharlap_raw" / (parent["job_id"] + ".mat")
        raw = mat(rawpath); state = raw["ray_state_vec"]
        xyz = np.column_stack([vector(state, k) for k in ("pos_x", "pos_y", "pos_z")])/1000.
        heights = np.linalg.norm(xyz, axis=1)-R_KM
        index = int(np.flatnonzero(heights >= 80.)[0])
        position = xyz[index]; h = float(heights[index]); unit = position/np.linalg.norm(position)
        lat = float(np.rad2deg(np.arcsin(unit[2]))); lon = float(np.rad2deg(np.arctan2(unit[1], unit[0])))
        p = np.array([vector(state, k)[index] for k in ("dir_x", "dir_y", "dir_z")])
        elevation = float(np.rad2deg(np.arcsin(np.dot(p/np.linalg.norm(p), unit))))
        ne = parent["Nb_m3"] + parent["g_m3_per_km"]*(h-parent["h_boundary_km"])
        n2 = 1.-P1_C*ne/(parent["frequency_MHz"]*1e6)**2
        for kind in ("correct", "raw"):
            for level in (1, 2, 3):
                job = dict(parent)
                job.update({"job_id": f"restart_{kind}_"+parent["job_id"].rsplit("_t",1)[0]+f"_t{level}",
                    "case_id": "restart_"+kind, "p1_case": "restart_"+kind,
                    "level": level, "tol": TOL[level-1], "origin_lat_deg": lat,
                    "origin_lon_deg": lon, "origin_height_km": h, "elevation_deg": elevation,
                    "initial_xyz_km": position, "initial_p": p if kind=="raw" else p/np.linalg.norm(p)*np.sqrt(n2),
                    "apply_projection": kind != "raw", "parent_raw": str(rawpath),
                    "parent_raw_sha256": sha(rawpath), "parent_state_index": index})
                jobs.append(_p1_job_hash(job))
    _p1_save_jobs(folder, jobs)
    print(f"P1 restart prepared {len(jobs)}", flush=True)


def _p1_stats(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float).reshape(-1)
    finite = values[np.isfinite(values)]
    return {"count": int(values.size), "finite_count": int(finite.size),
            "median": float(np.median(finite)) if finite.size else None,
            "min": float(np.min(finite)) if finite.size else None,
            "max": float(np.max(finite)) if finite.size else None,
            "rms": float(np.sqrt(np.mean(finite**2))) if finite.size else None}


def _p1_medium(grid: dict, heights: np.ndarray) -> np.ndarray:
    """Frozen input medium, including its explicit vacuum outside the grid."""
    from scipy.interpolate import CubicHermiteSpline
    h = np.asarray(grid["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(grid["ne_profile_m3"], dtype=float).reshape(-1)
    values = np.asarray(heights, dtype=float)
    spline = CubicHermiteSpline(h, ne, np.gradient(ne, h, edge_order=2))
    return np.where((values >= h[0]) & (values <= h[-1]),
                    np.maximum(0., spline(values)), 0.)


def _p1_segment(heights: np.ndarray, grid: dict) -> np.ndarray:
    h = np.asarray(grid["height_km"], dtype=float).reshape(-1)
    values = np.asarray(heights, dtype=float)
    return np.select([values < h[0], values < h[1], values <= h[-1]],
                     ["below", "first_cell", "interior"], default="above")


def _p1_uniform_exact(job: dict, samples: np.ndarray) -> dict:
    """Independent straight segments and spherical Snell interface; no fit."""
    if job["p1_case"] != "uniform" or "initial_p" in job:
        raise ValueError("analytic layer solution requires the prescribed vacuum launch")
    xyz0 = np.asarray(geocentric_to_ecef(
        job.get("origin_lat_deg", 30.), job.get("origin_lon_deg", 120.),
        job["origin_height_km"], R_KM), dtype=float)
    direction = np.asarray(launch_direction_ecef(
        job.get("origin_lat_deg", 30.), job.get("origin_lon_deg", 120.),
        job["elevation_deg"], 0., R_KM), dtype=float)
    direction /= np.linalg.norm(direction)
    radius = R_KM + job["h_boundary_km"]
    projection = float(xyz0 @ direction)
    discriminant = projection**2 + radius**2 - float(xyz0 @ xyz0)
    if discriminant < 0:
        raise ValueError("vacuum ray does not meet the prescribed spherical interface")
    crossing = -projection + np.sqrt(discriminant)
    position = xyz0 + crossing * direction
    normal = position / np.linalg.norm(position)
    tangential_p = direction - float(direction @ normal) * normal
    n2 = 1. - P1_C * job["Nb_m3"] / (job["frequency_MHz"] * 1e6)**2
    normal_p2 = n2 - float(tangential_p @ tangential_p)
    if crossing < 0 or n2 <= 0 or normal_p2 <= 0:
        raise ValueError("analytic transmission unavailable (cutoff or total reflection)")
    outgoing_p = tangential_p + np.sqrt(normal_p2) * normal
    outgoing = outgoing_p / np.sqrt(n2)
    samples = np.asarray(samples, dtype=float).reshape(-1)
    entered = samples >= crossing
    xyz = xyz0 + samples[:, None] * direction
    xyz[entered] = position + (samples[entered] - crossing)[:, None] * outgoing
    p = np.tile(direction, (samples.size, 1))
    p[entered] = outgoing_p
    group = np.minimum(samples, crossing) + np.maximum(samples - crossing, 0.) / np.sqrt(n2)
    return {"s": samples, "xyz": xyz, "p": p, "group": group,
            "crossing_s_km": float(crossing), "crossing_xyz_km": position,
            "incoming_p": direction, "outgoing_p": outgoing_p,
            "normal": normal, "tangential_p": tangential_p, "layer_n2": float(n2),
            "solution": "exact_vacuum_and_uniform_straight_segments_spherical_Snell"}


def p1_checks() -> dict:
    """Small analytic checks, independent of MATLAB and the AH integrator."""
    job = {"p1_case": "uniform", "origin_height_km": 50.,
           "elevation_deg": 22., "h_boundary_km": 60.,
           "Nb_m3": 5e8, "frequency_MHz": 10.}
    exact = _p1_uniform_exact(job, np.linspace(0., 100., 101))
    radial = exact["normal"]
    outgoing = exact["outgoing_p"]
    vacuum = _p1_uniform_exact({**job, "Nb_m3": 0.}, exact["s"])
    origin = np.asarray(geocentric_to_ecef(30., 120., 50., R_KM))
    straight = origin + vacuum["s"][:, None] * vacuum["incoming_p"]
    angular = np.cross(exact["xyz"], exact["p"])
    grid = {"height_km": np.array([60., 62., 64.]),
            "ne_profile_m3": np.full(3, 5e8)}
    checks = {
        "vacuum_straight_position": bool(np.max(np.abs(vacuum["xyz"] - straight)) < 2e-12),
        "vacuum_group_equals_geometric": bool(np.array_equal(vacuum["group"], vacuum["s"])),
        "interface_tangential_p_conserved": bool(np.linalg.norm(
            outgoing - (outgoing @ radial)*radial - exact["tangential_p"]) < 5e-15),
        "interface_outgoing_p2_equals_n2": bool(abs(outgoing @ outgoing - exact["layer_n2"]) < 5e-15),
        "spherical_angular_momentum_conserved": bool(np.max(np.linalg.norm(
            angular - angular[0], axis=1)) < 5e-12),
        "known_vacuum_below_and_above_grid": bool(np.array_equal(
            _p1_medium(grid, np.array([59., 60., 61., 65.])), [0., 5e8, 5e8, 0.])),
        "medium_segment_boundaries": bool(np.array_equal(_p1_segment(
            np.array([59., 60., 61., 62., 65.]), grid),
            ["below", "first_cell", "first_cell", "interior", "above"])),
    }
    return {"schema": "c1m-p1-analytic-checks-v1", "checks": checks,
            "passed": bool(all(checks.values())), "coefficient_C": P1_C,
            "crossing_s_km": exact["crossing_s_km"], "fitted_parameters": False}


def _p1_ph_read(path: Path, job: dict) -> dict:
    raw = mat(path)
    for key in ("job_sha256", "grid_sha256", "exporter_sha256"):
        if key not in raw or str(raw[key]) != str(job[key]):
            raise ValueError("PH frozen identity mismatch: " + key)
    state, output = raw["ray_state_vec"], raw["ray_path_data"]
    xyz = np.column_stack([vector(state, k) for k in ("pos_x", "pos_y", "pos_z")]) / 1000.
    p = np.column_stack([vector(state, k) for k in ("dir_x", "dir_y", "dir_z")])
    s_absolute = vector(state, "geometrical_path") / 1000.
    group_absolute = vector(state, "group_path") / 1000.
    lat, lon, height = (vector(output, k) for k in ("lat", "lon", "height"))
    path_xyz = np.asarray(geocentric_to_ecef(lat, lon, height, R_KM))
    s_path = vector(output, "geometric_distance")
    group_path = vector(output, "group_range")
    ne = vector(output, "electron_density") * 1e6
    n, ng = (vector(output, k) for k in ("refractive_index", "group_refractive_index"))
    count = s_absolute.size
    arrays = (xyz, p, group_absolute, lat, lon, height, path_xyz,
              s_path, group_path, ne, n, ng)
    if any(v.shape[0] != count for v in arrays):
        raise ValueError("PH path/state point counts differ")
    if (np.any(np.diff(s_absolute) <= 0.) or np.any(np.diff(group_absolute) <= 0.)
            or np.any(np.diff(s_path) <= 0.) or np.any(np.diff(group_path) <= 0.)):
        raise ValueError("PH native coordinates do not strictly increase")
    for key in ("geomag_x", "geomag_y", "geomag_z", "collision_frequency"):
        if np.any(vector(output, key) != 0.):
            raise ValueError("P1 B0/collisionless contract violated")
    return {"s": s_absolute - s_absolute[0], "xyz": xyz, "p": p,
            "group": group_absolute - group_absolute[0],
            "s_absolute": s_absolute, "group_absolute": group_absolute,
            "path_xyz": path_xyz, "path_height": height,
            "path_lat": lat, "path_lon": lon,
            "path_s_absolute": s_path, "path_group_absolute": group_path,
            "ne_out": ne, "n_out": n, "ng_out": ng,
            "height": np.linalg.norm(xyz, axis=1) - R_KM,
            "termination": "ground" if scalar(raw["ray_data"], "ray_label") == 1. else "other",
            "ray_label": scalar(raw["ray_data"], "ray_label"),
            "reported_event_group_km": scalar(raw["ray_data"], "group_range"),
            "reported_event_geometric_km": scalar(raw["ray_data"], "geometric_path_length"),
            "reported_event_lat_deg": scalar(raw["ray_data"], "lat"),
            "reported_event_lon_deg": scalar(raw["ray_data"], "lon"),
            "mex_sha256": str(raw.get("mex_sha256", "")), "raw_sha256": sha(path)}


def _p1_ah_read(folder: Path, job: dict) -> tuple[dict | None, dict | None]:
    path = folder / (job["job_id"] + ".json")
    if not path.exists():
        return None, None
    meta = read_json(path)
    if (meta.get("job_sha256") != job["job_sha256"]
            or meta.get("grid_sha256") != job["grid_sha256"]
            or meta.get("coefficient_K2") != P1_C or meta.get("interpolation") != "hermite"):
        raise ValueError("AH frozen identity/medium mismatch")
    raw_path = folder / (job["job_id"] + ".npz")
    if not raw_path.exists():
        return meta, None
    if sha(raw_path) != meta.get("raw_npz_sha256"):
        raise ValueError("AH raw hash mismatch")
    with np.load(raw_path) as raw:
        s = np.asarray(raw["geom_path_km"], dtype=float)
        xyz = np.asarray(raw["xyz_km"], dtype=float)
        p = np.asarray(raw["wavevector_km_inv"], dtype=float) * C_KM_S / (2*np.pi*job["frequency_MHz"]*1e6)
        ne = np.asarray(raw["Ne_m3_native"], dtype=float)
        n2 = np.asarray(raw["n2_native"], dtype=float)
    if (s.ndim != 1 or xyz.shape != (s.size, 3) or p.shape != xyz.shape
            or ne.shape != s.shape or n2.shape != s.shape
            or any(np.any(~np.isfinite(v)) for v in (s, xyz, p, ne, n2))
            or np.any(np.diff(s) <= 0.)):
        raise ValueError("invalid AH raw native arrays")
    return meta, {"s": s, "xyz": xyz, "p": p, "ne_out": ne,
                  "n2_native": n2, "height": np.linalg.norm(xyz, axis=1) - R_KM,
                  "termination": meta.get("termination_status"),
                  "group_total_km": meta.get("group_path_km_native")}


def _p1_native_diagnostics(trace: dict, job: dict, grid: dict, *, ph: bool) -> tuple[dict, dict]:
    xyz, p, s = trace["xyz"], trace["p"], trace["s"]
    height = trace["height"]
    p2 = np.sum(p*p, axis=1)
    input_ne = _p1_medium(grid, height)
    coefficient = P1_C / (job["frequency_MHz"]*1e6)**2
    input_n2 = 1. - coefficient * input_ne
    angular = np.cross(xyz, p)
    angular_norm = np.linalg.norm(angular, axis=1)
    angular_delta = np.linalg.norm(angular - angular[0], axis=1)
    segment = _p1_segment(height, grid)
    selected = np.ones(s.size, dtype=bool) if job["scope"] == "full" else s <= S_KM
    values = {"p2": p2, "input_ne_state_m3": input_ne, "input_n2_state": input_n2,
              "R_input_state": p2 - input_n2,
              "angular_momentum_norm_km": angular_norm,
              "angular_momentum_vector_delta_km": angular_delta,
              "angular_momentum_relative_delta": angular_delta / max(angular_norm[0], 1e-30)}
    summary = {"sample_scope": "full_native" if job["scope"] == "full" else "native_s_le_100km",
               "point_count": int(s.size), "diagnostic_point_count": int(np.sum(selected)),
               "first_geometric_km": float(s[0]), "last_geometric_km": float(s[-1]),
               "below_grid_points": int(np.sum(selected & (segment == "below")))}
    if ph:
        path_ne = _p1_medium(grid, trace["path_height"])
        n2_out = trace["n_out"]**2
        density_n2 = 1. - coefficient * trace["ne_out"]
        values.update({"input_ne_path_m3": path_ne, "input_n2_path": 1.-coefficient*path_ne,
                       "R_output_n": p2-n2_out, "R_output_density": p2-density_n2,
                       "output_n2_minus_density_n2": n2_out-density_n2,
                       "R_input_path": p2-(1.-coefficient*path_ne),
                       "Ne_out_minus_input_state_m3": trace["ne_out"]-input_ne,
                       "Ne_out_minus_input_path_m3": trace["ne_out"]-path_ne,
                       "path_state_position_m": np.linalg.norm(trace["path_xyz"]-xyz, axis=1)*1000.,
                       "path_state_geometric_m": (trace["path_s_absolute"]-trace["s_absolute"])*1000.,
                       "path_state_group_m": (trace["path_group_absolute"]-trace["group_absolute"])*1000.})
        summary["path_height_segments"] = {name: int(np.sum(selected & (_p1_segment(
            trace["path_height"], grid) == name))) for name in ("below", "first_cell", "interior", "above")}
        # Native interval diagnostics; discard the terminal interval and every
        # interval touching an entrance/cell crossing. No gradient is inferred.
        stable = (segment[:-1] == segment[1:]) & selected[:-1] & selected[1:]
        crossing = np.flatnonzero(segment[:-1] != segment[1:])
        for index in crossing:
            stable[max(0, index-1):min(stable.size, index+2)] = False
        stable[-1:] = False
        native_ng = np.diff(trace["group"])/np.diff(s)
        native_path_ng = np.diff(trace["path_group_absolute"])/np.diff(trace["path_s_absolute"])
        def average(v):
            return (v[:-1]+v[1:])/2.
        with np.errstate(invalid="ignore", divide="ignore"):
            interval_values = {
                "native_dgroup_dgeom": native_ng,
                "path_dgroup_dgeom": native_path_ng,
                "output_ng_average": average(trace["ng_out"]),
                "inv_p_average": average(1./np.sqrt(p2)),
                "inv_output_n_average": average(1./trace["n_out"]),
                "inv_input_n_average": average(1./np.sqrt(input_n2)),
                "offshell_hamiltonian_ng_candidate": average((p2+coefficient*trace["ne_out"])/np.sqrt(p2)),
                "native_ng_minus_output_ng": native_ng-average(trace["ng_out"]),
                "native_ng_minus_inv_p": native_ng-average(1./np.sqrt(p2)),
                "native_ng_minus_inv_output_n": native_ng-average(1./trace["n_out"]),
                "native_ng_minus_inv_input_n": native_ng-average(1./np.sqrt(input_n2)),
                "native_ng_minus_offshell_candidate": native_ng-average((p2+coefficient*trace["ne_out"])/np.sqrt(p2)),
            }
        summary["native_group_interval_segments"] = {
            name: {key: _p1_stats(v[stable & (segment[:-1] == name)])
                   for key, v in interval_values.items()}
            for name in ("below", "first_cell", "interior", "above")}
        summary["native_group_difference_semantics"] = (
            "Delta(native group)/Delta(native geometric); endpoint means are diagnostics, "
            "not exact quadrature or reconstructed group states")
        values.update({key: np.r_[np.where(stable, v, np.nan), np.nan]
                       for key, v in interval_values.items()})
        summary["index_lag_diagnostic"] = {}
        for lag in (-1, 0, 1):
            index = np.arange(max(0, -lag), min(s.size, s.size-lag))
            use = selected[index] & (index > 0) & (index+lag < s.size-1)
            index = index[use]
            summary["index_lag_diagnostic"][str(lag)] = {
                name: {"p2_minus_output_n2": _p1_stats((p2[index]-n2_out[index+lag])[
                    segment[index] == name]),
                       "state_vs_shifted_path_position_m": _p1_stats(np.linalg.norm(
                           xyz[index]-trace["path_xyz"][index+lag], axis=1)[segment[index] == name]*1000.)}
                for name in ("below", "first_cell", "interior")}
        summary["native_terminal"] = {
            "ray_label": trace["ray_label"], "reported_event_kind": trace["termination"],
            "state_height_km": float(height[-1]), "path_height_km": float(trace["path_height"][-1]),
            "state_xyz_km": xyz[-1], "path_xyz_km": trace["path_xyz"][-1],
            "path_state_position_m": float(values["path_state_position_m"][-1]),
            "path_state_geometric_m": float(values["path_state_geometric_m"][-1]),
            "path_state_group_m": float(values["path_state_group_m"][-1]),
            "state_total_group_km": float(trace["group_absolute"][-1]),
            "path_total_group_km": float(trace["path_group_absolute"][-1]),
            "reported_event_group_km": trace["reported_event_group_km"],
            "reported_event_geometric_km": trace["reported_event_geometric_km"]}
    else:
        values["R_exported_n2"] = p2-trace["n2_native"]
        summary["native_terminal"] = {"termination": trace["termination"],
                                       "state_height_km": float(height[-1]),
                                       "state_xyz_km": xyz[-1],
                                       "group_total_km": trace["group_total_km"],
                                       "pointwise_native_group_available": False}
    summary["segments"] = {name: {key: _p1_stats(v[selected & (segment == name)])
                                    for key, v in values.items() if v.shape == s.shape}
                           for name in ("below", "first_cell", "interior", "above")}
    values["segment"] = segment
    values["selected"] = selected
    if job["p1_case"] == "uniform" and "initial_p" not in job:
        exact = _p1_uniform_exact(job, s[selected])
        position_error = np.linalg.norm(xyz[selected]-exact["xyz"], axis=1)*1000.
        p_error = np.linalg.norm(p[selected]-exact["p"], axis=1)
        summary["independent_analytic_native"] = {
            "solution": exact["solution"], "no_fitted_coefficients": True,
            "crossing_s_km": exact["crossing_s_km"], "layer_n2": exact["layer_n2"],
            "position_error_m": _p1_stats(position_error), "p_vector_error": _p1_stats(p_error),
            "outgoing_p": exact["outgoing_p"], "tangential_p": exact["tangential_p"]}
        values["analytic_position_error_m"] = np.full(s.size, np.nan)
        values["analytic_p_vector_error"] = np.full(s.size, np.nan)
        values["analytic_position_error_m"][selected] = position_error
        values["analytic_p_vector_error"][selected] = p_error
    return summary, values


def _p1_pair_metrics(a: dict, b: dict, scope: str) -> dict:
    """Position reconstruction is diagnostic until its own bound is established."""
    end = S_KM if scope == "local" else min(float(a["s"][-1]), float(b["s"][-1]))
    samples = np.linspace(0., end, int(np.ceil(end/.25))+1)
    result = norm_metrics(resample(a["s"], a["xyz"], samples)-resample(b["s"], b["xyz"], samples))
    endpoints = [resample(t["s"], t["xyz"], np.array([end]))[0] if scope == "local"
                 else t["xyz"][-1] for t in (a, b)]
    groups = []
    for trace in (a, b):
        if "group" in trace:
            groups.append(float(resample(trace["s"], trace["group"], np.array([end]))[0])
                          if scope == "local" else float(trace["group"][-1]))
        else:
            if scope == "local" and (trace["termination"] != "fixed_geometric_path_completed"
                                     or abs(float(trace["s"][-1])-S_KM) > 1e-8):
                raise ValueError("AH native group total is not the required 100 km section")
            groups.append(float(trace["group_total_km"]))
    return {"rms_m": result["rms_m"], "max_m": result["max_m"],
            "landing_m": float(np.linalg.norm(endpoints[0]-endpoints[1])*1000.),
            "group_rel": abs(groups[0]-groups[1])/max(abs(groups[0]), 1e-15)}


def _p1_flatten(value: dict, prefix: str = "") -> dict:
    result = {}
    for key, item in value.items():
        name = prefix+key
        if isinstance(item, dict):
            result.update(_p1_flatten(item, name+"."))
        elif isinstance(item, (list, tuple, np.ndarray)):
            result[name] = json.dumps(np.asarray(item).tolist(), allow_nan=False)
        else:
            result[name] = item
    return result


def stage_p1_summary(run: Path, *, dataset: str, subset: str, fd: float, analysis_id: str) -> None:
    """Read raw evidence into a new immutable analysis version, including rejects."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", analysis_id):
        raise ValueError("analysis-id must be 1-64 safe characters")
    folder = run / "restart" if dataset == "restart" else run
    analysis = folder / "analysis" / (f"p1_fd{fd:g}_"+analysis_id)
    analysis.mkdir(parents=True, exist_ok=False)
    sources = analysis / "sources"; sources.mkdir()
    frozen = read_json(folder / "p1_manifest.json")
    jobs_doc = read_json(folder / "jobs.json")
    jobs = [j for j in jobs_doc["jobs"] if subset == "all" or j["p1_case"] == subset or j["role"] == subset]
    if not jobs:
        raise ValueError("empty P1 analysis selection")
    ah_dir = folder / "ah" / f"fd{fd:g}"
    evidence = {"p1_manifest": folder / "p1_manifest.json", "jobs": folder / "jobs.json",
                "input": folder / "input.mat", "grids": run / "grids.json"}
    if (folder / "pharlap_summary.json").is_file():
        evidence["pharlap_summary"] = folder / "pharlap_summary.json"
    preparation_root = run
    source_records = {}
    for key, record in frozen["sources"].items():
        path = Path(record["path"])
        snapshot = preparation_root / record["snapshot"] if record.get("snapshot") else path
        source_records[key] = {**record, "frozen_source_verified": bool(snapshot.is_file()
            and sha(snapshot) == record["sha256"]),
            "current_source_matches_preparation": bool(path.is_file() and sha(path) == record["sha256"])}
        if snapshot.is_file() and snapshot.suffix.lower() not in {".mat", ".exe", ".mexw64"}:
            shutil.copyfile(snapshot, sources / ("prepared_"+key+snapshot.suffix))
    for key, path in {"analysis_runner": Path(__file__),
                      "analysis_verifier": ROOT / "tests/_verify_a4_c1m_attribution.py"}.items():
        if path.is_file():
            target = sources / (key+path.suffix); shutil.copyfile(path, target)
            source_records[key] = {"path": str(path), "sha256": sha(path),
                                   "snapshot": str(target.relative_to(analysis)),
                                   "frozen_source_verified": bool(sha(target) == sha(path))}
    for job in jobs:
        evidence["grid_"+job["grid_sha256"][:16]] = Path(job["grid_mat"])
        for label, path in (("ph", folder / "pharlap_raw" / (job["job_id"]+".mat")),
                            ("ah_meta", ah_dir / (job["job_id"]+".json")),
                            ("ah_raw", ah_dir / (job["job_id"]+".npz"))):
            if path.is_file():
                evidence[label+"_"+job["job_id"]] = path
        if "parent_raw" in job:
            evidence["restart_parent_"+job["parent_raw_sha256"][:16]] = Path(job["parent_raw"])
    records = {key: {"path": str(path), "sha256": sha(path)} for key, path in evidence.items()}
    write_json(analysis / "analysis_manifest.json", {
        "schema": "c1m-p1-analysis-freeze-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": dataset, "subset": subset, "analysis_id": analysis_id,
        "sources": source_records, "evidence": records,
        "snapshot_before_analysis": True, "ray_computation": False, "production_modified": False,
        "python_version": sys.version, "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__})
    checks = p1_checks(); write_json(analysis / "analytic_checks.json", checks)
    if not checks["passed"]:
        raise RuntimeError("P1 analytic self-check failed; failed analysis retained")
    point_fields = ["job_id", "solver", "native_index", "s_km", "segment", "height_state_km",
                    "height_path_km", "state_x_km", "state_y_km", "state_z_km",
                    "path_x_km", "path_y_km", "path_z_km", "p_x", "p_y", "p_z", "p2",
                    "Ne_out_m3", "n_out", "ng_out", "group_state_km", "group_path_km",
                    "input_ne_state_m3", "input_ne_path_m3", "input_n2_state", "input_n2_path",
                    "R_input_state", "R_input_path", "R_output_n", "R_output_density",
                    "R_exported_n2", "output_n2_minus_density_n2",
                    "Ne_out_minus_input_state_m3", "Ne_out_minus_input_path_m3",
                    "path_state_position_m", "path_state_geometric_m", "path_state_group_m",
                    "angular_momentum_x_km", "angular_momentum_y_km", "angular_momentum_z_km",
                    "angular_momentum_norm_km", "angular_momentum_vector_delta_km",
                    "angular_momentum_relative_delta", "native_dgroup_dgeom", "path_dgroup_dgeom",
                    "output_ng_average", "inv_p_average", "inv_output_n_average", "inv_input_n_average",
                    "offshell_hamiltonian_ng_candidate", "native_ng_minus_offshell_candidate",
                    "native_ng_minus_output_ng", "native_ng_minus_inv_p",
                    "native_ng_minus_inv_output_n", "native_ng_minus_inv_input_n",
                    "analytic_position_error_m", "analytic_p_vector_error"]
    rows, traces, densities, grids = [], {}, [], {}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    analytic_fig, analytic_axes = plt.subplots(2, 2, figsize=(12, 8))
    with (analysis / "per_point.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=point_fields); writer.writeheader()
        for job in jobs:
            row = {key: job.get(key) for key in ("job_id", "case_id", "p1_case", "frequency_MHz",
                   "elevation_deg", "level", "h_boundary_km", "Nb_m3", "g_m3_per_km",
                   "width_km", "dh_km", "scope", "role", "job_sha256", "grid_sha256")}
            row["raw_restart_diagnostic_only"] = bool(job["p1_case"] == "restart_raw")
            row["physical_acceptance"] = False
            grid = grids.setdefault(job["grid_sha256"], mat(Path(job["grid_mat"])))
            identity_errors = []
            if _p1_job_hash(job)["job_sha256"] != row["job_sha256"]:
                identity_errors.append("job_hash")
            if sha(Path(job["grid_mat"])) != job["grid_sha256"]:
                identity_errors.append("grid_hash")
            if "parent_raw" in job and sha(Path(job["parent_raw"])) != job["parent_raw_sha256"]:
                identity_errors.append("restart_parent_hash")
            ph_path = folder / "pharlap_raw" / (job["job_id"]+".mat")
            ph, ah, ah_meta = None, None, None
            try:
                if ph_path.is_file():
                    ph = _p1_ph_read(ph_path, job)
                    if ph["mex_sha256"] != frozen["sources"]["pharlap_mex"]["sha256"]:
                        raise ValueError("PH raw MEX source mismatch")
                    row["ph_status"] = "raw_available"
                    row["ph_raw_sha256"] = ph["raw_sha256"]
                else:
                    row["ph_status"] = "missing_raw"
            except Exception as exc:
                row.update({"ph_status": "invalid_raw", "ph_error": str(exc)}); ph = None
                identity_errors.append("ph_raw")
            try:
                ah_meta, ah = _p1_ah_read(ah_dir, job)
                row["ah_status"] = ah_meta.get("status") if ah_meta else "missing_result"
                row["ah_termination"] = ah_meta.get("termination_status") if ah_meta else None
                row["ah_failure"] = ah_meta.get("failure_reason", ah_meta.get("error")) if ah_meta else None
                if ah_meta:
                    expected = {alias: frozen["sources"][key]["sha256"] for alias, key in
                                (("runner", "runner"), ("config", "config"),
                                 ("dispersion", "ah_dispersion"), ("tracer", "ah_solver"),
                                 ("provider", "ah_provider"))}
                    if ah_meta.get("source_sha256") != expected or ah_meta.get("spatial_fd_base_km") != fd:
                        raise ValueError("AH execution source/FD mismatch")
                    row["ah_raw_sha256"] = ah_meta.get("raw_npz_sha256")
            except Exception as exc:
                row.update({"ah_status": "invalid_raw", "ah_error": str(exc)}); ah = None
                identity_errors.append("ah_raw")
            row["identity_passed"] = bool(not identity_errors)
            row["identity_errors"] = identity_errors
            for solver, trace in (("ph", ph), ("ah", ah)):
                if trace is None:
                    continue
                diagnostics, values = _p1_native_diagnostics(trace, job, grid, ph=solver == "ph")
                row[solver+"_diagnostics"] = diagnostics
                angular = np.cross(trace["xyz"], trace["p"])
                for index in np.flatnonzero(values["selected"]):
                    point = {"job_id": job["job_id"], "solver": solver, "native_index": int(index),
                             "s_km": trace["s"][index], "segment": values["segment"][index],
                             "height_state_km": trace["height"][index], "Ne_out_m3": trace["ne_out"][index]}
                    point.update(dict(zip(("state_x_km", "state_y_km", "state_z_km"), trace["xyz"][index])))
                    point.update(dict(zip(("p_x", "p_y", "p_z"), trace["p"][index])))
                    point.update(dict(zip(("angular_momentum_x_km", "angular_momentum_y_km",
                                           "angular_momentum_z_km"), angular[index])))
                    if solver == "ph":
                        point.update({"height_path_km": trace["path_height"][index], "n_out": trace["n_out"][index],
                                      "ng_out": trace["ng_out"][index], "group_state_km": trace["group_absolute"][index],
                                      "group_path_km": trace["path_group_absolute"][index]})
                        point.update(dict(zip(("path_x_km", "path_y_km", "path_z_km"), trace["path_xyz"][index])))
                    for key, value in values.items():
                        if key in point_fields and key != "segment":
                            point[key] = value[index] if np.isfinite(value[index]) else None
                    writer.writerow(point)
                if solver == "ph":
                    x_boundary = P1_C*job["Nb_m3"]/(job["frequency_MHz"]*1e6)**2
                    for name in ("below", "first_cell", "interior"):
                        med = diagnostics["segments"][name]["R_output_n"]["median"]
                        densities.append({"job_id": job["job_id"], "segment": name,
                                          "Nb_m3": job["Nb_m3"], "frequency_MHz": job["frequency_MHz"],
                                          "R_median": med, "X_boundary": x_boundary,
                                          "R_over_X_boundary": med/x_boundary if med is not None and x_boundary > 0. else None,
                                          "fitted_coefficients": False})
                if (solver == "ph" and job["level"] == 3 and job["frequency_MHz"] == 10.
                        and job["dh_km"] == 2. and job["role"] != "cache_control"):
                    local = trace["s"] <= S_KM
                    label = job["case_id"]+f" e{job['elevation_deg']:g}"
                    axes.flat[0].plot(trace["s"][local], values["R_output_n"][local], label=label, linewidth=.8)
                    axes.flat[1].plot(trace["s"][local], values["Ne_out_minus_input_path_m3"][local], linewidth=.8)
                    axes.flat[2].plot(trace["s"][local], values["native_ng_minus_inv_p"][local], linewidth=.8)
                    axes.flat[3].plot(trace["s"][local], values["angular_momentum_relative_delta"][local], linewidth=.8)
                    if job["p1_case"] == "uniform":
                        panel = 0 if job["elevation_deg"] < 30. else 1
                        analytic_axes.flat[panel].plot(trace["s"][local], values["analytic_position_error_m"][local],
                                                       label=label, linewidth=.8)
                        analytic_axes.flat[2].plot(trace["s"][local], values["R_output_n"][local], label=label, linewidth=.8)
                        analytic_axes.flat[3].plot(trace["s"][local], values["native_ng_minus_output_ng"][local], linewidth=.8)
            if ph is not None and ah is not None:
                origin = np.asarray(job.get("initial_xyz_km", geocentric_to_ecef(
                    job.get("origin_lat_deg", 30.), job.get("origin_lon_deg", 120.), job["origin_height_km"], R_KM)))
                row["ph_origin_error_m"] = float(np.linalg.norm(ph["xyz"][0]-origin)*1000.)
                row["ah_origin_error_m"] = float(np.linalg.norm(ah["xyz"][0]-origin)*1000.)
                row["initial_p_difference"] = float(np.linalg.norm(ph["p"][0]-ah["p"][0]))
                if job["scope"] == "full":
                    row["full_natural_terminal_comparison"] = {
                        "ph_reported_ground": bool(ph["termination"] == "ground"),
                        "ah_ground": bool(ah["termination"] == "ground"),
                        "state_endpoint_distance_m": float(np.linalg.norm(ph["xyz"][-1]-ah["xyz"][-1])*1000.),
                        "path_endpoint_vs_ah_distance_m": float(np.linalg.norm(ph["path_xyz"][-1]-ah["xyz"][-1])*1000.),
                        "ph_state_minus_ah_group_km": float(ph["group"][-1]-ah["group_total_km"]),
                        "ph_event_minus_ah_group_km": float(ph["reported_event_group_km"]-ah["group_total_km"]),
                        "same_natural_termination": bool(ph["termination"] == ah["termination"])}
                try:
                    row["cross_tool_diagnostic"] = _p1_pair_metrics(ph, ah, job["scope"])
                except ValueError as exc:
                    row["cross_tool_unavailable_reason"] = str(exc)
            elif job["scope"] == "full":
                row["full_natural_terminal_comparison"] = {
                    "status": "not_evaluated", "reason": "both raw native endpoints required"}
            rows.append(row)
            traces[job["job_id"]] = {"job": job, "ph": ph, "ah": ah, "meta": ah_meta}
    labels = ("p2 - output n2", "Ne output - Hermite input at path (m^-3)",
              "native Delta group/Delta geometric - mean 1/|p|", "relative angular momentum vector change")
    for ax, label in zip(axes.flat, labels):
        ax.set_xlabel("native geometric arc (km)"); ax.set_ylabel(label); ax.grid(alpha=.25)
    axes.flat[0].legend(fontsize=5, ncol=2)
    fig.tight_layout(); fig.savefig(analysis / "radial_state_semantics.png", dpi=160); plt.close(fig)
    for ax, label in zip(analytic_axes.flat, ("analytic position error, 22 deg (m)",
                        "analytic position error, 39.5 deg (m)", "p2 - output n2",
                        "native Delta group/Delta geometric - output ng")):
        ax.set_xlabel("native geometric arc (km)"); ax.set_ylabel(label); ax.grid(alpha=.25)
    analytic_axes.flat[0].legend(fontsize=6)
    analytic_fig.tight_layout(); analytic_fig.savefig(analysis / "uniform_spherical_interface.png", dpi=160)
    plt.close(analytic_fig)
    groups = {}
    for row in rows:
        key = (row["case_id"], row["frequency_MHz"], row["elevation_deg"], row["role"])
        groups.setdefault(key, []).append(row)
    decisions = []
    for key, group in groups.items():
        levels = sorted(group, key=lambda item: item["level"])
        reasons = []
        if len(levels) != 3 or [r["level"] for r in levels] != [1, 2, 3]:
            reasons.append("three_frozen_levels_unavailable")
        if any(not row["identity_passed"] for row in levels):
            reasons.append("raw_source_or_input_identity_failed")
        if any(row.get("ah_status") != "completed" or "cross_tool_diagnostic" not in row for row in levels):
            reasons.append("both_tools_complete_with_required_scope_unavailable")
        if any(row.get("initial_p_difference", 1.) > 1e-10 for row in levels):
            reasons.append("common_physical_initial_state_unavailable")
        if levels[0]["p1_case"] == "restart_raw":
            reasons.append("raw_offshell_restart_is_diagnostic_only")
        if levels[0]["scope"] == "full" and any(not row.get(
                "full_natural_terminal_comparison", {}).get("same_natural_termination", False) for row in levels):
            reasons.append("natural_terminal_topology_mismatch")
        decision = {"case": list(key), "status": "not_evaluated", "passed": False,
                    "physical_acceptance": False, "reasons": reasons}
        if not reasons or reasons == ["raw_offshell_restart_is_diagnostic_only"]:
            try:
                native = [traces[row["job_id"]] for row in levels]
                increments = [{side: _p1_pair_metrics(native[a][side], native[b][side], levels[0]["scope"])
                               for side in ("ph", "ah")} for a, b in ((0, 1), (1, 2))]
                zero = {name: 0. for name in (*LIMIT, "landing_m")}
                decision["zero_resampling_budget_diagnostic"] = numeric_decision_v2(
                    levels[-1]["cross_tool_diagnostic"], increments[0], increments[1], zero, eligible=False)
                decision["increments_12"] = increments[0]
                decision["increments_23"] = increments[1]
            except (ValueError, TypeError) as exc:
                reasons.append("numeric_evidence_unavailable: "+str(exc))
        reasons.append("independent_native_trajectory_reconstruction_error_bound_not_established")
        decision["reconstruction_method_difference_in_budget"] = False
        decision["roundoff_floor_used"] = False
        decisions.append(decision)
    flat = [_p1_flatten(row) for row in rows]
    with (analysis / "per_case.csv").open("w", encoding="utf-8", newline="") as stream:
        fields = list(dict.fromkeys(key for row in flat for key in row))
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(flat)
    cache = []
    cache_jobs = [j for j in jobs if j["role"] == "cache_control"]
    for first, later in ((0, 2), (1, 3)):
        if len(cache_jobs) > later:
            a, b = (traces[cache_jobs[i]["job_id"]]["ph"] for i in (first, later))
            cache.append({"jobs": [cache_jobs[first]["job_id"], cache_jobs[later]["job_id"]],
                          "raw_available": bool(a is not None and b is not None),
                          "native_arrays_identical": bool(a is not None and b is not None and all(
                              np.array_equal(a[key], b[key]) for key in
                              ("s", "xyz", "p", "group", "path_xyz", "ne_out", "n_out", "ng_out")))})
    integrity = bool(all(r["frozen_source_verified"] for r in source_records.values())
                     and all(row["identity_passed"] for row in rows)
                     and jobs_doc["input_sha256"] == sha(folder / "input.mat"))
    summary = {"schema": "c1m-p1-state-attribution-v1", "dataset": dataset, "subset": subset,
               "coefficient_C": P1_C, "fitted_coefficients": False, "earth_radius_km": R_KM,
               "evidence_integrity_passed": integrity,
               "source_verified": integrity, "physics_passed": False,
               "analytic_self_checks_passed": checks["passed"],
               "all_requested_ph_raw_available": bool(all(r.get("ph_status") == "raw_available" for r in rows)),
               "all_requested_ah_completed": bool(all(r.get("ah_status") == "completed" for r in rows)),
               "physical_acceptance_passed": False, "numeric_acceptance_passed": False,
               "production_modified": False, "new_ray_computation": False,
               "analysis_manifest_sha256": sha(analysis / "analysis_manifest.json"),
               "case_count": len(rows), "cases": rows, "three_level_decisions": decisions,
               "density_dependence_no_fit": densities, "cache_order_controls": cache,
               "interpretation": ["Radial, spherical, B=0, collisionless input-medium attribution only.",
                   "Native state p, path n/Ne and event outputs retain separate semantics.",
                   "Lag +/-1 and finite native group differences are diagnostics, not a corrected physical model.",
                   "AH NPZ has no native pointwise group path; only its metadata total is used at an eligible endpoint.",
                   "Incomplete/native-rejected and raw offshell restart cases cannot receive physical acceptance.",
                   "Equivalent on-shell Hamiltonians can give different offshell group velocities; raw-restart differences do not attribute an AH physics error.",
                   "PCHIP comparisons are diagnostic until P2 supplies independent reconstruction-error evidence."]}
    write_json(analysis / "summary.json", summary)
    (analysis / "report.txt").write_text(
        f"P1 径向硬入口状态语义分析：{dataset}，{len(rows)} 条记录。\n"
        f"证据完整性：{integrity}；解析自检：{checks['passed']}。\n"
        f"PH 原始数据完整：{summary['all_requested_ph_raw_available']}；"
        f"AH 全部原生完成：{summary['all_requested_ah_completed']}。\n"
        "逐点保存 path/state 坐标、程、波矢、输入 Hermite 介质、输出 Ne/n 和角动量。\n"
        "按低于网格、首单元和内部统计中位数、范围及 RMS；局部固定 s<=100 km，全程 C1 单列自然终点。\n"
        "群程差分排除入口/单元交叉和末端区间，比较输出 ng、1/|p|、输入 1/n。\n"
        "均匀层独立解使用球面切向 p 守恒和法向 p=sqrt(n²-p_t²)，未拟合任何系数。\n"
        "三档严格判据的零重采样预算仅作诊断；重建误差尚无独立证据，物理/数值验收均未评定。\n"
        "原始离壳重启仅诊断。AH 原生拒绝、缺失和失败均保留；本步未运行射线、未修改生产实现。\n",
        encoding="utf-8")
    write_json(analysis / "sha256_manifest.json", {"schema": "c1m-p1-analysis-artifacts-v1",
        "records": [{"path": str(path), "sha256": sha(path)}
                    for path in sorted(analysis.rglob("*")) if path.is_file()]})
    print(f"P1 summary: {len(rows)} cases; integrity={integrity}; analytic={checks['passed']}; {analysis}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=("p2-prepare", "p2-analyze", "p2-exit", "p2-verify", "p1-prepare", "p1-ah", "p1-restart-prepare", "p1-summary", "p0", "a0", "pharlap", "a1", "a2-prepare", "a2-identified-prepare", "a2-identified-holdout-prepare", "a2-extra-prepare", "a2-ah", "a2-summary", "a2-artifacts"))
    parser.add_argument("--p1-subset", default="all")
    parser.add_argument("--p1-dataset", choices=("base", "restart"), default="base")
    parser.add_argument("--p1-fd-km", type=float, default=.125)
    parser.add_argument("--p1-analysis-id", default="v1")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--a2-dataset", choices=("base", "identified", "identified_holdout", "extra"), default="base")
    parser.add_argument("--a2-extra-id")
    parser.add_argument("--a2-subset", choices=("all", "boundary", "smooth", "internal", "affine"), default="all")
    parser.add_argument("--a2-workers", type=int, default=1)
    parser.add_argument("--a2-k2", type=float)
    parser.add_argument("--a2-interp", choices=("linear", "hermite"), default="linear")
    parser.add_argument("--a2-job-id")
    parser.add_argument("--a2-tag", choices=("pilot", "formal"), default="formal")
    parser.add_argument("--a2-fd-base-km", type=float, default=0.5)
    parser.add_argument("--a2-level", type=int)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        parser.error("run-id must be 1-64 letters, digits, underscore or hyphen")
    run = OUT / args.run_id
    if args.stage == "p2-prepare":
        from _a4_c1m_p2_analysis import prepare
        prepare(run)
    elif args.stage == "p2-analyze":
        from _a4_c1m_p2_analysis import analyze
        analyze(run)
    elif args.stage == "p2-exit":
        from _a4_c1m_p2_exit import run as run_exit
        run_exit(run)
    elif args.stage == "p2-verify":
        from _a4_c1m_p2_verify import verify
        if verify(run):
            return 1
    elif args.stage == "p1-prepare":
        stage_p1_prepare(run)
    elif args.stage == "p1-restart-prepare":
        stage_p1_restart_prepare(run)
    elif args.stage == "p1-ah":
        stage_p1_ah(run if args.p1_dataset == "base" else run / "restart",
                    args.a2_workers, args.p1_subset, args.p1_fd_km)
    elif args.stage == "p1-summary":
        stage_p1_summary(run, dataset=args.p1_dataset, subset=args.p1_subset,
                         fd=args.p1_fd_km, analysis_id=args.p1_analysis_id)
    elif args.stage == "p0":
        stage_p0(run)
    elif args.stage == "a0":
        stage_a0(run)
    elif args.stage == "pharlap":
        stage_pharlap(run)
    elif args.stage == "a2-prepare":
        stage_a2_prepare(run)
    elif args.stage == "a2-identified-prepare":
        stage_a2_identified_prepare(run)
    elif args.stage == "a2-identified-holdout-prepare":
        stage_a2_identified_holdout_prepare(run)
    elif args.stage == "a2-extra-prepare":
        stage_a2_extra_prepare(run, dataset=args.a2_dataset,
                               job_id=args.a2_job_id)
    elif args.stage == "a2-ah":
        if args.a2_workers < 1:
            parser.error("a2-workers must be positive")
        stage_a2_ah(run, dataset=args.a2_dataset, subset=args.a2_subset,
                    workers=args.a2_workers, requested_k2=args.a2_k2,
                    interpolation=args.a2_interp,
                    job_id=args.a2_job_id, tag=args.a2_tag,
                    extra_id=args.a2_extra_id,
                    fd_base_km=args.a2_fd_base_km, level=args.a2_level)
    elif args.stage == "a2-summary":
        stage_a2_summary(run, dataset=args.a2_dataset,
                         interpolation=args.a2_interp, requested_k2=args.a2_k2,
                         tag=args.a2_tag, subset=args.a2_subset,
                         extra_id=args.a2_extra_id,
                         fd_base_km=args.a2_fd_base_km)
    elif args.stage == "a2-artifacts":
        stage_a2_artifacts(run)
    else:
        stage_a1(run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
