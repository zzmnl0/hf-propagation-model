"""Generate and freeze equivalent AH3D/PHaRLAP Chapman C1 inputs."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
from utils import spherical_destination


MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SCRIPT = MATLAB_DIR / "a4_c1_prep1_preflight.m"
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
OUT = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep1"
INPUT_DIR = OUT / "inputs"
PREFLIGHT_DIR = OUT / "preflight"
REPORT = OUT / "a4_c1_prep1_validation.json"
MANIFEST = OUT / "a4_c1_prep1_freeze_manifest.json"
FIGURE = OUT / "a4_c1_prep1_input_contract.png"

TARGET_RANGE_KM = 1185.848260273658
TX = (30.0, 120.0, 0.0)
BEARING_DEG = 0.0
FREQUENCY_MHZ = 10.0
FOF2_MHZ = 7.0
HMF2_KM = 350.0
YMF2_KM = 100.0
LATITUDE_MIN_DEG = 28.0
LATITUDE_MAX_DEG = 65.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _array_hash(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array, dtype="<f8")
    digest = hashlib.sha256()
    digest.update(json.dumps(list(value.shape)).encode("ascii"))
    digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist()
    ) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def _code_hash(paths: tuple[Path, ...]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).replace("\\", "/").encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _surface_distance_km(lat1: float, lon1: float, lat2: float, lon2: float,
                         radius_km: float) -> float:
    lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])
    cosine = (np.sin(lat1) * np.sin(lat2)
              + np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1))
    return float(radius_km * np.arccos(np.clip(cosine, -1.0, 1.0)))


def _axes(height_step_km: float, height_stop_km: float
          ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    latitude = np.arange(
        LATITUDE_MIN_DEG, LATITUDE_MAX_DEG + 0.125, 0.25, dtype=np.float64)
    longitude = 118.0 + 0.25 * np.arange(17, dtype=np.float64)
    count = int(round((height_stop_km - 60.0) / height_step_km)) + 1
    height = 60.0 + height_step_km * np.arange(count, dtype=np.float64)
    return latitude, longitude, height


def _chapman_grid(height_km: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    peak = (FOF2_MHZ * 1.0e6 / 8.98) ** 2
    u = (height_km - HMF2_KM) / YMF2_KM
    profile = peak * np.exp(np.clip(
        0.5 * (1.0 - u - np.exp(np.clip(-u, -700.0, 700.0))),
        -745.0, 700.0))
    return np.broadcast_to(profile, shape).copy()


def _build_input(height_step_km: float, height_stop_km: float,
                 variant: str, path: Path,
                 common_hashes: dict) -> dict:
    radius = float(cfg.EARTH_MODEL["radius_km"])
    rx_lat, rx_lon = spherical_destination(
        TX[0], TX[1], BEARING_DEG, TARGET_RANGE_KM, radius)
    latitude, longitude, height = _axes(height_step_km, height_stop_km)
    shape = (latitude.size, longitude.size, height.size)
    ne_m3 = _chapman_grid(height, shape)
    ne_cm3 = ne_m3 / 1.0e6
    collision = np.zeros(shape, dtype=np.float64)
    grid = np.array([
        latitude[0], latitude[1] - latitude[0], latitude.size,
        longitude[0], longitude[1] - longitude[0], longitude.size,
        height[0], height[1] - height[0], height.size,
    ], dtype=np.float64)
    b_latitude = np.arange(
        LATITUDE_MIN_DEG, LATITUDE_MAX_DEG + 0.25, 0.5, dtype=np.float64)
    b_longitude = 118.0 + 0.5 * np.arange(9, dtype=np.float64)
    b_height = 60.0 + 5.0 * np.arange(89, dtype=np.float64)
    b_shape = (b_latitude.size, b_longitude.size, b_height.size)
    b_zero = np.zeros(b_shape, dtype=np.float64)
    b_grid = np.array([
        b_latitude[0], b_latitude[1] - b_latitude[0], b_latitude.size,
        b_longitude[0], b_longitude[1] - b_longitude[0], b_longitude.size,
        b_height[0], b_height[1] - b_height[0], b_height.size,
    ], dtype=np.float64)
    physics = {
        "model": "Chapman_F2_horizontal_uniform",
        "foF2_MHz": FOF2_MHZ,
        "hmF2_km": HMF2_KM,
        "ymF2_km": YMF2_KM,
        "frequency_MHz": FREQUENCY_MHZ,
        "B_T": 0.0,
        "collision_Hz": 0.0,
        "earth_radius_km": radius,
        "tx_geocentric": list(TX),
        "rx_geocentric": [rx_lat, rx_lon, 0.0],
        "target_range_km": TARGET_RANGE_KM,
        "bearing_deg": BEARING_DEG,
    }
    field_hash = _json_hash({
        "physics": physics,
        "grid_variant": variant,
        "latitude_hash": _array_hash(latitude),
        "longitude_hash": _array_hash(longitude),
        "height_hash": _array_hash(height),
        "Ne_m3_hash": _array_hash(ne_m3),
    })
    provider = GeocentricGridNeProvider3D(
        ne_m3, latitude, longitude, height, radius,
        source_metadata={"field_hash": field_hash, **physics})
    provider_metadata = provider.get_metadata()
    provider_hash = provider_metadata["grid_hash"]
    units = {
        "frequency_MHz": "MHz", "earth_radius_m": "m",
        "latitude_deg": "degree_geocentric",
        "longitude_deg": "degree_east",
        "height_km": "km_above_spherical_earth",
        "iono_en_grid_m3": "electrons/m^3",
        "iono_en_grid_cm3": "electrons/cm^3",
        "collision_freq_Hz": "Hz", "Bx_T": "T", "By_T": "T", "Bz_T": "T",
    }
    data = {
        "schema_version": "a4-c1-prep1-input-v1",
        "case_id": f"a4_c1_chapman_{variant}",
        "case_semantics": "equivalent_input_freeze_not_external_physics_validation",
        "grid_variant": variant,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "frequency_MHz": FREQUENCY_MHZ,
        "wave_mode": "magnetoionic_degenerate_single_count",
        "OX_mode": 0.0,
        "nhops": 1.0,
        "tol": np.array([1.0e-8, 0.01, 10.0], dtype=np.float64),
        "earth_model": "sphere_user_radius",
        "earth_radius_m": radius * 1000.0,
        "coordinate_mapping_version": "pharlap_geocentric_spherical_to_ecef_v1",
        "origin_lat_deg": TX[0], "origin_lon_deg": TX[1],
        "origin_height_km": TX[2], "bearing_deg": BEARING_DEG,
        "receiver_lat_deg": rx_lat, "receiver_lon_deg": rx_lon,
        "receiver_height_km": 0.0, "target_range_km": TARGET_RANGE_KM,
        "foF2_MHz": FOF2_MHZ, "hmF2_km": HMF2_KM, "ymF2_km": YMF2_KM,
        "elevation_scan_deg": np.arange(1.0, 89.0 + 0.25, 0.5, dtype=np.float64),
        "legacy_angle_usage": "not_present_not_used",
        "latitude_deg": latitude, "longitude_deg": longitude, "height_km": height,
        "geomag_latitude_deg": b_latitude, "geomag_longitude_deg": b_longitude,
        "geomag_height_km": b_height,
        "iono_en_grid_m3": ne_m3,
        "iono_en_grid_cm3": ne_cm3,
        "iono_en_grid_5_cm3": ne_cm3.copy(),
        "collision_freq_Hz": collision,
        "iono_grid_parms": grid,
        "Bx_T": b_zero.copy(), "By_T": b_zero.copy(), "Bz_T": b_zero.copy(),
        "geomag_grid_parms": b_grid,
        "array_axis_order": "latitude,longitude,height",
        "units_json": json.dumps(units, sort_keys=True),
        "provider_metadata_json": json.dumps(provider_metadata, sort_keys=True),
        "field_hash": field_hash, "provider_hash": provider_hash,
        "Ne_m3_hash": _array_hash(ne_m3), "Ne_cm3_hash": _array_hash(ne_cm3),
        **common_hashes,
    }
    temporary = path.with_suffix(".tmp.mat")
    savemat(temporary, data, do_compression=True, long_field_names=True,
            appendmat=False)
    temporary.replace(path)

    peak = (FOF2_MHZ * 1.0e6 / 8.98) ** 2
    unit_error = float(np.max(np.abs(ne_cm3 * 1.0e6 - ne_m3)) / peak)
    same_point_error = abs(provider.get_Ne_geocentric(
        TX[0], TX[1], HMF2_KM) - peak) / peak
    rx_distance_error = abs(_surface_distance_km(
        TX[0], TX[1], rx_lat, rx_lon, radius) - TARGET_RANGE_KM)
    checks = {
        "shape": ne_m3.shape == shape and min(shape) >= 2,
        "axes_strictly_increasing": all(np.all(np.diff(axis) > 0.0)
                                         for axis in (latitude, longitude, height)),
        "finite_nonnegative": bool(np.all(np.isfinite(ne_m3)) and np.all(ne_m3 >= 0.0)),
        "peak_relative_error": abs(float(np.max(ne_m3)) - peak) / peak <= 1.0e-12,
        "same_point_provider_error": same_point_error <= 1.0e-12,
        "unit_roundtrip_error": unit_error <= 1.0e-12,
        "rx_surface_distance_error": rx_distance_error <= 1.0e-9,
        "longitude_domain": longitude[0] <= min(TX[1], rx_lon)
                            and longitude[-1] >= max(TX[1], rx_lon),
        "pharlap_grid_limits": all(
            actual <= limit for actual, limit in zip(shape, (701, 701, 401)))
            and all(actual <= limit for actual, limit in zip(b_shape, (101, 101, 201))),
        "zero_B_collision": not np.any(collision) and not np.any(b_zero),
    }
    return {
        "variant": variant, "path": str(path), "sha256": _sha256(path),
        "shape": list(shape), "geomag_shape": list(b_shape),
        "height_step_km": height_step_km,
        "field_hash": field_hash, "provider_hash": provider_hash,
        "Ne_m3_hash": data["Ne_m3_hash"], "Ne_cm3_hash": data["Ne_cm3_hash"],
        "provider_metadata": provider_metadata,
        "metrics": {
            "peak_Ne_m3": peak, "unit_roundtrip_relative_error": unit_error,
            "same_point_provider_relative_error": same_point_error,
            "rx_surface_distance_error_km": rx_distance_error,
            "receiver_geocentric": [rx_lat, rx_lon, 0.0],
        },
        "checks": checks,
    }


def _mat_string(value) -> str:
    return str(np.asarray(value).reshape(-1)[0])


def _roundtrip_checks(item: dict) -> dict:
    data = loadmat(item["path"], squeeze_me=False, struct_as_record=False)
    return {
        "mat_schema": _mat_string(data["schema_version"]) == "a4-c1-prep1-input-v1",
        "mat_axis_order": _mat_string(data["array_axis_order"])
                          == "latitude,longitude,height",
        "mat_Ne_m3_hash": _array_hash(data["iono_en_grid_m3"]) == item["Ne_m3_hash"],
        "mat_Ne_cm3_hash": _array_hash(data["iono_en_grid_cm3"]) == item["Ne_cm3_hash"],
        "mat_file_hash": _sha256(Path(item["path"])) == item["sha256"],
    }


def _plot_input_contract(data: dict, source_hash: str) -> dict:
    latitude = np.asarray(data["latitude_deg"]).reshape(-1)
    height = np.asarray(data["height_km"]).reshape(-1)
    ne = np.asarray(data["iono_en_grid_m3"], dtype=float)
    section = ne[:, ne.shape[1] // 2, :].T / 1.0e11
    tx_index = int(np.argmin(np.abs(latitude - TX[0])))
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5), constrained_layout=True)
    mesh = axes[0].pcolormesh(latitude, height, section, shading="auto",
                              cmap="jet", alpha=0.5)
    axes[0].axvline(TX[0], color="black", linestyle="--", linewidth=1.0,
                   label="TX latitude")
    axes[0].set(xlabel="Geocentric latitude (deg)", ylabel="Physical height (km)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.colorbar(mesh, ax=axes[0], label="Ne (1e11 m^-3)")
    axes[1].plot(section[:, tx_index], height, color="#0072B2", linewidth=1.7)
    axes[1].axhline(HMF2_KM, color="black", linestyle="--", linewidth=1.0,
                    label="hmF2")
    axes[1].set(xlabel="Ne (1e11 m^-3)", ylabel="Physical height (km)")
    axes[1].legend(frameon=False, fontsize=8)
    for label, axis in zip(("A", "B"), axes):
        axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top",
                  fontweight="bold")
    fig.suptitle(
        "A4-C1 Prep-1 input contract | Chapman F2 frozen grid\n"
        "one-way reference input | 10 MHz | B=0 | collision=0 | not ray validation",
        fontsize=9.5)
    fig.savefig(FIGURE, dpi=300, facecolor="white", metadata={
        "Title": "A4-C1 Prep-1 frozen Chapman input contract",
        "Description": "Visualization of the frozen Ne input; not ray validation",
    })
    plt.close(fig)
    image = mpimg.imread(FIGURE)
    checks = {
        "exists": FIGURE.is_file() and FIGURE.stat().st_size > 0,
        "finite_pixels": bool(np.all(np.isfinite(image))),
        "minimum_dimensions": bool(image.shape[0] >= 900 and image.shape[1] >= 1800),
        "nonblank_pixel_variance": bool(float(np.std(image[..., :3])) >= 1.0e-3),
        "source_hash_recorded": len(source_hash) == 64,
    }
    return {"path": str(FIGURE), "sha256": _sha256(FIGURE),
            "source_sha256": source_hash, "pixel_shape": list(image.shape),
            "pixel_std": float(np.std(image[..., :3])), "checks": checks}


def _run_matlab(main_path: Path, refined_path: Path, output_mat: Path,
                output_json: Path, log_path: Path) -> subprocess.CompletedProcess:
    quote = lambda path: str(path).replace("'", "''")
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); "
        f"a4_c1_prep1_preflight('{quote(main_path)}','{quote(refined_path)}',"
        f"'{quote(output_mat)}','{quote(output_json)}')")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run([str(MATLAB), "-batch", expression], env=env,
                            text=True, capture_output=True, timeout=300,
                            errors="replace")
    log_path.write_text((result.stdout or "") + "\n" + (result.stderr or ""),
                        encoding="utf-8")
    return result


def _verify_existing_manifest(allow_rebuild: bool) -> int | None:
    if not MANIFEST.is_file():
        return None
    frozen = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks = {
        "status": frozen.get("status") == "passed_input_contract_frozen",
        "artifacts_unchanged": all(
            Path(item["path"]).is_file() and _sha256(Path(item["path"])) == item["sha256"]
            for item in frozen.get("frozen_artifacts", [])),
    }
    if all(checks.values()):
        print("A4-C1-Prep-1 frozen: OK")
        return 0
    if allow_rebuild:
        print("A4-C1-Prep-1 frozen artifacts changed; explicit rebuild enabled")
        return None
    print("A4-C1-Prep-1 frozen artifacts changed; rebuild refused")
    return 1


def run(allow_rebuild: bool = False) -> int:
    existing = _verify_existing_manifest(allow_rebuild)
    if existing is not None:
        return existing
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREFLIGHT_DIR.mkdir(parents=True, exist_ok=True)
    prerequisites = {
        "matlab": MATLAB.is_file(), "pharlap_mex": MEX.is_file(),
        "matlab_script": MATLAB_SCRIPT.is_file(),
        "earth_model_sphere": cfg.EARTH_MODEL.get("name") == "sphere",
    }
    if not all(prerequisites.values()):
        _atomic_json(REPORT, {"status": "blocked_prerequisite", "checks": prerequisites})
        return 2

    code_paths = (Path(__file__), ROOT / "models" / "ne_provider_geocentric.py",
                  ROOT / "utils.py", ROOT / "config.py")
    common_hashes = {
        "python_generator_sha256": _sha256(Path(__file__)),
        "matlab_script_sha256": _sha256(MATLAB_SCRIPT),
        "pharlap_mex_sha256": _sha256(MEX),
        "code_hash": _code_hash(code_paths),
    }
    main = _build_input(
        2.0, 500.0, "main_2km", INPUT_DIR / "a4_c1_chapman_main_2km.mat",
        common_hashes)
    refined = _build_input(
        1.0, 460.0, "refined_1km", INPUT_DIR / "a4_c1_chapman_refined_1km.mat",
        common_hashes)
    main["checks"].update(_roundtrip_checks(main))
    refined["checks"].update(_roundtrip_checks(refined))
    main_data = loadmat(main["path"])
    refined_data = loadmat(refined["path"])
    figure = _plot_input_contract(main_data, main["sha256"])
    shared_grid_checks = {
        "same_latitude_axis": np.array_equal(
            main_data["latitude_deg"], refined_data["latitude_deg"]),
        "same_longitude_axis": np.array_equal(
            main_data["longitude_deg"], refined_data["longitude_deg"]),
        "refined_contains_shared_main_heights": np.array_equal(
            main_data["height_km"][:, :201], refined_data["height_km"][:, ::2]),
        "same_Ne_on_shared_heights": np.array_equal(
            main_data["iono_en_grid_m3"][:, :, :201],
            refined_data["iono_en_grid_m3"][:, :, ::2]),
    }
    preflight_mat = PREFLIGHT_DIR / "a4_c1_prep1_preflight.mat"
    preflight_json = PREFLIGHT_DIR / "a4_c1_prep1_preflight.json"
    log = OUT / "a4_c1_prep1_matlab.log"
    process = _run_matlab(Path(main["path"]), Path(refined["path"]),
                          preflight_mat, preflight_json, log)
    matlab_result = (json.loads(preflight_json.read_text(encoding="utf-8"))
                     if preflight_json.is_file() else {})
    checks = {
        **prerequisites,
        "main_input": all(main["checks"].values()),
        "refined_input": all(refined["checks"].values()),
        "shared_grid_contract": all(shared_grid_checks.values()),
        "visualization": all(figure["checks"].values()),
        "matlab_returncode": process.returncode == 0,
        "matlab_preflight": matlab_result.get("status")
                            == "passed_input_contract_preflight",
        "matlab_selected_pharlap_mex": str(
            matlab_result.get("raytrace_3d_sp_path", "")).lower().endswith(".mexw64"),
    }
    status = "passed_input_contract_frozen" if all(checks.values()) else "failed"
    report = {
        "schema_version": "a4-c1-prep1-validation-v1",
        "status": status,
        "case_semantics": "input_contract_only_not_ray_or_external_physics_validation",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks, "shared_grid_checks": shared_grid_checks,
        "inputs": [main, refined], "common_hashes": common_hashes,
        "software": {
            "matlab_executable": str(MATLAB), "pharlap_root": str(PHARLAP),
            "matlab_script": str(MATLAB_SCRIPT), "matlab_preflight": matlab_result,
        },
        "artifacts": {"preflight_mat": str(preflight_mat),
                      "preflight_json": str(preflight_json), "run_log": str(log),
                      "figure": figure},
        "limitations": [
            "Prep-1 does not trace rays or compare AH3D with PHaRLAP",
            "Legacy local-coordinate high/low angles are absent from the frozen input",
        ],
    }
    _atomic_json(REPORT, report)
    if status != "passed_input_contract_frozen":
        print(f"A4-C1-Prep-1: FAIL report={REPORT}")
        return 1
    frozen_paths = [Path(main["path"]), Path(refined["path"]), preflight_mat,
                    preflight_json, FIGURE, MATLAB_SCRIPT, MEX, Path(__file__)]
    manifest = {
        "schema_version": "a4-c1-prep1-freeze-v1", "status": status,
        "case_semantics": report["case_semantics"],
        "validation_report": {"path": str(REPORT), "sha256": _sha256(REPORT)},
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in frozen_paths],
    }
    _atomic_json(MANIFEST, manifest)
    print(f"A4-C1-Prep-1: OK report={REPORT}")
    return 0


if __name__ == "__main__":
    arguments = sys.argv[1:]
    if any(argument != "--rebuild-frozen" for argument in arguments):
        raise SystemExit("usage: _run_a4_c1_prep1.py [--rebuild-frozen]")
    raise SystemExit(run(allow_rebuild="--rebuild-frozen" in arguments))
