"""Run the A4-C0 PHaRLAP/AH3D vacuum interface and coordinate precheck."""

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
from models.ray_tracer_3d import RefractiveIndex3D
from models.ray_tracer_ah_3d import RefractiveIndexAH3D, trace_single_ray_ah_3d


MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MATLAB_SCRIPTS = Path(r"D:\Matlab\PHaRLAP")
OUTPUT = ROOT / "output2" / "a4_pharlap"
CASE_ID = "a4_c0_vacuum"
R_KM = 6371.0
ORIGIN_LAT_DEG = 30.0
ORIGIN_LON_DEG = 120.0
ELEVATION_DEG = 35.0
BEARING_DEG = 37.0
FREQUENCY_MHZ = 10.0
COMPARE_PATH_KM = 250.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True),
                         encoding="utf-8")
    temporary.replace(path)


def _matlab_quote(path: Path) -> str:
    return str(path).replace("'", "''")


def _as_1d(value) -> np.ndarray:
    return np.asarray(value, dtype=float).reshape(-1)


def _ecef_from_lat_lon_height(lat_deg, lon_deg, height_km) -> np.ndarray:
    lat = np.radians(_as_1d(lat_deg))
    lon = np.radians(_as_1d(lon_deg))
    radius = R_KM + _as_1d(height_km)
    return np.column_stack((
        radius * np.cos(lat) * np.cos(lon),
        radius * np.cos(lat) * np.sin(lon),
        radius * np.sin(lat),
    ))


def _ah_local_to_ecef(trajectory: np.ndarray) -> np.ndarray:
    x, y, z_cart = trajectory[:, 0], trajectory[:, 1], trajectory[:, 2]
    rho = np.hypot(x, y)
    sigma = rho / R_KM
    bearing = np.arctan2(y, x)
    lat1, lon1 = np.radians([ORIGIN_LAT_DEG, ORIGIN_LON_DEG])
    sin_lat = (np.sin(lat1) * np.cos(sigma)
               + np.cos(lat1) * np.sin(sigma) * np.cos(bearing))
    lat = np.arcsin(np.clip(sin_lat, -1.0, 1.0))
    lon = lon1 + np.arctan2(
        np.sin(bearing) * np.sin(sigma) * np.cos(lat1),
        np.cos(sigma) - np.sin(lat1) * np.sin(lat),
    )
    height = z_cart - rho * rho / (2.0 * R_KM)
    radius = R_KM + height
    return np.column_stack((radius * np.cos(lat) * np.cos(lon),
                            radius * np.cos(lat) * np.sin(lon),
                            radius * np.sin(lat)))


def _cumulative_distance(points: np.ndarray) -> np.ndarray:
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]


def _resample(points: np.ndarray, distance: np.ndarray,
              target: np.ndarray) -> np.ndarray:
    keep = np.r_[True, np.diff(distance) > 1.0e-10]
    return np.column_stack([
        np.interp(target, distance[keep], points[keep, axis]) for axis in range(3)
    ])


def _launch_frame() -> tuple[np.ndarray, np.ndarray]:
    lat, lon = np.radians([ORIGIN_LAT_DEG, ORIGIN_LON_DEG])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                   np.sin(lat)])
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
                      np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    elevation, bearing = np.radians([ELEVATION_DEG, BEARING_DEG])
    direction = (np.cos(elevation) * (np.cos(bearing) * north
                                      + np.sin(bearing) * east)
                 + np.sin(elevation) * up)
    return R_KM * up, direction / np.linalg.norm(direction)


def _ecef_to_enu(points: np.ndarray) -> np.ndarray:
    lat, lon = np.radians([ORIGIN_LAT_DEG, ORIGIN_LON_DEG])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
                      np.cos(lat)])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                   np.sin(lat)])
    origin = R_KM * up
    return (points - origin) @ np.column_stack((east, north, up))


def _generate_input(path: Path) -> dict:
    lat_start, lat_step, n_lat = 15.0, 0.5, 61
    lon_start, lon_step, n_lon = 105.0, 0.5, 61
    height_start, height_step, n_height = 60.0, 5.0, 101
    shape = (n_lat, n_lon, n_height)
    zero = np.zeros(shape, dtype=np.float64)
    grid = np.array([lat_start, lat_step, n_lat, lon_start, lon_step, n_lon,
                     height_start, height_step, n_height], dtype=np.float64)
    field_contract = {
        "case_id": CASE_ID,
        "shape": shape,
        "electron_density_unit": "electrons/cm^3",
        "collision_frequency_unit": "Hz",
        "geomagnetic_field_unit": "T",
        "axis_order": "latitude,longitude,height",
        "Ne": 0.0,
        "B": 0.0,
        "collision": 0.0,
    }
    data = {
        "schema_version": "a4-c0-v1",
        "case_id": CASE_ID,
        "case_semantics": "vacuum_interface_coordinate_precheck_not_physics_gate",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "frequency_MHz": FREQUENCY_MHZ,
        "wave_mode": "no_field",
        "OX_mode": 0.0,
        "nhops": 1.0,
        "tol": np.array([1.0e-8, 0.01, 10.0], dtype=np.float64),
        "earth_model": "sphere_user_radius",
        "earth_radius_m": R_KM * 1000.0,
        "coordinate_mapping_version": "a4-c0-local-earth-flattening-to-spherical-v1",
        "origin_lat_deg": ORIGIN_LAT_DEG,
        "origin_lon_deg": ORIGIN_LON_DEG,
        "origin_height_km": 0.0,
        "elevation_deg": ELEVATION_DEG,
        "bearing_deg": BEARING_DEG,
        "iono_en_grid_cm3": zero,
        "iono_en_grid_5_cm3": zero.copy(),
        "collision_freq_Hz": zero.copy(),
        "iono_grid_parms": grid,
        "Bx_T": zero.copy(),
        "By_T": zero.copy(),
        "Bz_T": zero.copy(),
        "geomag_grid_parms": grid.copy(),
        "field_hash": _stable_hash(field_contract),
        "units_json": json.dumps({
            "frequency_MHz": "MHz", "earth_radius_m": "m",
            "origin_height_km": "km", "angles": "degree",
            "iono_en_grid_cm3": "electrons/cm^3", "collision_freq_Hz": "Hz",
            "Bx_T": "T", "By_T": "T", "Bz_T": "T",
        }, sort_keys=True),
        "array_axis_order": "latitude,longitude,height",
    }
    savemat(path, data, do_compression=True, long_field_names=True)
    return {"field_hash": data["field_hash"], "shape": list(shape)}


def _run_matlab(input_path: Path, preflight_mat: Path, preflight_json: Path,
                raw_path: Path, log_path: Path) -> subprocess.CompletedProcess:
    expression = (
        f"addpath('{_matlab_quote(MATLAB_SCRIPTS)}'); "
        f"a4_pharlap_preflight('{_matlab_quote(input_path)}', "
        f"'{_matlab_quote(preflight_mat)}', '{_matlab_quote(preflight_json)}', "
        f"'{_matlab_quote(raw_path)}')"
    )
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run([str(MATLAB), "-batch", expression], env=env,
                            text=True, capture_output=True, timeout=600)
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    return result


def _run_ah3d(raw_path: Path) -> dict:
    x = np.linspace(-50.0, 1800.0, 186)
    z = np.linspace(0.0, 700.0, 141)
    n3d = RefractiveIndex3D(np.zeros((x.size, z.size)), x, z, FREQUENCY_MHZ,
                            earth_flat=True)
    model = RefractiveIndexAH3D(
        n3d, FREQUENCY_MHZ, "O",
        geomag={**cfg.GEOMAG, "fH_MHz": 0.0}, integration_profile="reference")
    ray = trace_single_ray_ah_3d(
        (0.0, 0.0), ELEVATION_DEG, BEARING_DEG, model, FREQUENCY_MHZ,
        {**cfg.RT, "z_stop_km": 550.0, "max_path_km": 1800.0}, t=0.0)
    trajectory = np.asarray(ray.pop("trajectory_3d"), dtype=float)
    record = {key: (value.item() if isinstance(value, np.generic) else value)
              for key, value in ray.items()}
    record.update({
        "case_id": CASE_ID,
        "case_semantics": "AH3D_v2_vacuum_earth_flattening_coordinate_precheck",
        "trajectory_3d": trajectory.tolist(),
        "coordinate_model": "local_cartesian_earth_flattening",
        "magnetic_field_semantics": "B=0_magnetoionic_degenerate_single_count",
    })
    _atomic_json(raw_path, record)
    return record


def _compare(pharlap_raw: Path, ah: dict) -> tuple[dict, dict, dict]:
    raw = loadmat(pharlap_raw, simplify_cells=True)
    path = raw["ray_path_data"]
    ray_data = raw["ray_data"]
    ph_ecef = _ecef_from_lat_lon_height(path["lat"], path["lon"], path["height"])
    ph_distance = _as_1d(path["geometric_distance"])
    ph_distance -= ph_distance[0]
    ah_local = np.asarray(ah["trajectory_3d"], dtype=float)
    ah_ecef = _ah_local_to_ecef(ah_local)
    ah_distance = _cumulative_distance(ah_ecef)
    common = min(COMPARE_PATH_KM, float(ph_distance[-1]), float(ah_distance[-1]))
    sample = np.linspace(0.0, common, 201)
    ph_sample = _resample(ph_ecef, ph_distance, sample)
    ah_sample = _resample(ah_ecef, ah_distance, sample)
    residual = np.linalg.norm(ph_sample - ah_sample, axis=1)
    origin, direction = _launch_frame()
    analytic = origin[None, :] + sample[:, None] * direction[None, :]
    ph_analytic = np.linalg.norm(ph_sample - analytic, axis=1)
    ah_analytic = np.linalg.norm(ah_sample - analytic, axis=1)
    label = int(_as_1d(ray_data["ray_label"])[-1])
    label_map = {1: "ground", -2: "penetrated", -3: "grid_exit",
                 -4: "max_points", -100: "catastrophic_error"}
    ph_status = label_map.get(label, f"unmapped_{label}")
    ah_status = str(ah["branch_status"])
    broad_status_ok = (ph_status in {"penetrated", "grid_exit"}
                       and ah_status == "upper_boundary")
    metrics = {
        "comparison_path_km": common,
        "trajectory_rms_3d_km": float(np.sqrt(np.mean(residual * residual))),
        "trajectory_max_3d_km": float(np.max(residual)),
        "pharlap_analytic_chord_rms_km": float(np.sqrt(np.mean(ph_analytic ** 2))),
        "pharlap_analytic_chord_max_km": float(np.max(ph_analytic)),
        "ah3d_analytic_chord_rms_km": float(np.sqrt(np.mean(ah_analytic ** 2))),
        "ah3d_analytic_chord_max_km": float(np.max(ah_analytic)),
        "origin_error_km": float(np.linalg.norm(ph_sample[0] - ah_sample[0])),
        "pharlap_geometric_path_km": float(ph_distance[-1]),
        "ah3d_ecef_reconstructed_geometric_path_km": float(ah_distance[-1]),
        "ah3d_local_geometric_path_km": float(ah["geom_path_km"]),
        "pharlap_group_path_km": float(_as_1d(path["group_range"])[-1]),
        "ah3d_group_path_km": float(ah["group_path_km"]),
        "pharlap_vacuum_group_vs_geom_relative": float(abs(
            _as_1d(path["group_range"])[-1] - ph_distance[-1])
            / max(ph_distance[-1], 1.0e-12)),
        "ah3d_vacuum_group_vs_local_geom_relative": float(abs(
            float(ah["group_path_km"]) - float(ah["geom_path_km"]))
            / max(float(ah["geom_path_km"]), 1.0e-12)),
        "pharlap_ray_label": label,
        "pharlap_status": ph_status,
        "ah3d_status": ah_status,
    }
    checks = {
        "finite_coordinates": bool(np.all(np.isfinite(ph_ecef))
                                   and np.all(np.isfinite(ah_ecef))),
        "non_catastrophic_pharlap": ph_status != "catastrophic_error",
        "non_ground_exit_semantics_preserved": broad_status_ok,
        "origin_matches_without_registration": metrics["origin_error_km"] <= 1.0e-6,
        "coordinate_direction_and_scale": (
            metrics["trajectory_rms_3d_km"] <= 2.0
            and metrics["trajectory_max_3d_km"] <= 5.0),
        "pharlap_vacuum_straight_line": metrics["pharlap_analytic_chord_max_km"] <= 0.1,
        "ah3d_vacuum_mapping_consistent": metrics["ah3d_analytic_chord_max_km"] <= 5.0,
        "comparison_interval_complete": common >= COMPARE_PATH_KM - 1.0e-6,
    }
    return metrics, checks, {
        "ph_ecef": ph_ecef, "ph_distance": ph_distance,
        "ah_ecef": ah_ecef, "ah_distance": ah_distance,
        "sample": sample, "ph_sample": ph_sample, "ah_sample": ah_sample,
        "ph_status": ph_status, "ah_status": ah_status,
    }


def _style_axis(ax) -> None:
    ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)


def _plot_c0(data: dict, figures_dir: Path, source_hashes: dict) -> dict:
    figures_dir.mkdir(parents=True, exist_ok=True)
    ph_ecef, ah_ecef = data["ph_ecef"], data["ah_ecef"]
    ph_s, ah_s = data["ph_distance"], data["ah_distance"]
    ph_enu, ah_enu = _ecef_to_enu(ph_ecef), _ecef_to_enu(ah_ecef)
    ph_h = np.linalg.norm(ph_ecef, axis=1) - R_KM
    ah_h = np.linalg.norm(ah_ecef, axis=1) - R_KM
    origin, direction = _launch_frame()

    def analytic(s):
        ecef = origin[None, :] + s[:, None] * direction[None, :]
        return _ecef_to_enu(ecef), np.linalg.norm(ecef, axis=1) - R_KM

    ph_analytic_enu, ph_analytic_h = analytic(ph_s)
    ah_analytic_enu, ah_analytic_h = analytic(ah_s)
    max_s = 1.03 * max(float(ph_s[-1]), float(ah_s[-1]))
    max_h = 1.05 * max(float(np.max(ph_h)), float(np.max(ah_h)))
    all_e = np.r_[ph_enu[:, 0], ah_enu[:, 0]]
    all_n = np.r_[ph_enu[:, 1], ah_enu[:, 1]]
    e_pad = max(5.0, 0.05 * float(np.ptp(all_e)))
    n_pad = max(5.0, 0.05 * float(np.ptp(all_n)))
    plan_limits = ((float(np.min(all_e) - e_pad), float(np.max(all_e) + e_pad)),
                   (float(np.min(all_n) - n_pad), float(np.max(all_n) + n_pad)))
    specs = (
        ("pharlap", ph_s, ph_h, ph_enu, ph_analytic_h, ph_analytic_enu,
         "#0072B2", "PHaRLAP raw", data["ph_status"],
         "PHaRLAP 4.7.4 raytrace_3d_sp | spherical Earth"),
        ("ah3d", ah_s, ah_h, ah_enu, ah_analytic_h, ah_analytic_enu,
         "#D55E00", "AH3D raw", data["ah_status"],
         "AH3D V2 DOP853 reference | local earth-flattening"),
    )
    paths = {}
    for key, s, height, enu, analytic_h, analytic_enu, color, label, status, solver in specs:
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 5.0), constrained_layout=True)
        markevery = max(1, len(s) // 20)
        axes[0].plot(s, analytic_h, "--", color="#000000", linewidth=1.2,
                     label="Analytic vacuum chord")
        axes[0].plot(s, height, color=color, marker="o", markersize=2.8,
                     markevery=markevery, linewidth=1.5, label=label)
        axes[0].set(xlabel="Geometric path (km)", ylabel="Physical height (km)",
                    xlim=(0.0, max_s), ylim=(0.0, max_h))
        axes[0].legend(frameon=False, fontsize=8)
        axes[1].plot(analytic_enu[:, 0], analytic_enu[:, 1], "--", color="#000000",
                     linewidth=1.2, label="Analytic vacuum chord")
        axes[1].plot(enu[:, 0], enu[:, 1], color=color, marker="o", markersize=2.8,
                     markevery=markevery, linewidth=1.5, label=label)
        axes[1].set(xlabel="East (km)", ylabel="North (km)",
                    xlim=plan_limits[0], ylim=plan_limits[1])
        axes[1].set_aspect("equal", adjustable="box")
        axes[1].legend(frameon=False, fontsize=8)
        for panel, ax in zip(("A", "B"), axes):
            _style_axis(ax)
            ax.text(0.01, 0.98, panel, transform=ax.transAxes, va="top",
                    fontweight="bold")
        fig.suptitle(
            f"A4-C0 | one-way vacuum | {solver}\n"
            "10 MHz | OX=0 | B=0 T | collision=0 Hz\n"
            f"status={status} | scattering/power not evaluated",
            fontsize=9.2,
        )
        path = figures_dir / f"trajectory_c0_{key}.png"
        fig.savefig(path, dpi=300, facecolor="white", metadata={
            "Title": f"A4-C0 {key} vacuum trajectory",
            "Description": "One-way vacuum interface/coordinate precheck; not physical validation",
        })
        plt.close(fig)
        paths[key] = path

    sample = data["sample"]
    ph_common, ah_common = data["ph_sample"], data["ah_sample"]
    ph_common_enu, ah_common_enu = _ecef_to_enu(ph_common), _ecef_to_enu(ah_common)
    ph_common_h = np.linalg.norm(ph_common, axis=1) - R_KM
    ah_common_h = np.linalg.norm(ah_common, axis=1) - R_KM
    residual_enu = ah_common_enu - ph_common_enu
    residual_3d = np.linalg.norm(ah_common - ph_common, axis=1)
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5), constrained_layout=True)
    axes[0].plot(sample, ph_common_h, color="#0072B2", linewidth=1.6,
                 label="PHaRLAP")
    axes[0].plot(sample, ah_common_h, "--", color="#D55E00", linewidth=1.6,
                 label="AH3D")
    axes[0].set(xlabel="Common geometric path (km)", ylabel="Physical height (km)")
    axes[1].plot(ph_common_enu[:, 0], ph_common_enu[:, 1], color="#0072B2",
                 linewidth=1.6, label="PHaRLAP")
    axes[1].plot(ah_common_enu[:, 0], ah_common_enu[:, 1], "--", color="#D55E00",
                 linewidth=1.6, label="AH3D")
    axes[1].set(xlabel="East (km)", ylabel="North (km)")
    axes[1].set_aspect("equal", adjustable="box")
    axes[2].plot(sample, residual_3d, color="#000000", linewidth=1.8,
                 label="3D norm")
    for component, color, label in zip(range(3), ("#009E73", "#CC79A7", "#E69F00"),
                                       ("East", "North", "Up")):
        axes[2].plot(sample, np.abs(residual_enu[:, component]), color=color,
                     linewidth=1.0, label=f"|{label}|")
    axes[2].set(xlabel="Common geometric path (km)", ylabel="AH3D-PHaRLAP error (km)")
    for panel, ax in zip(("A", "B", "C"), axes):
        _style_axis(ax)
        ax.legend(frameon=False, fontsize=8)
        ax.text(0.01, 0.98, panel, transform=ax.transAxes, va="top", fontweight="bold")
    fig.suptitle(
        "A4-C0 one-way vacuum comparison | no trajectory registration\n"
        "PHaRLAP sphere vs AH3D earth-flattening | 10 MHz | B=0 | collision=0 | "
        "first 250 km | not an ionospheric physics gate",
        fontsize=10,
    )
    comparison = figures_dir / "trajectory_c0_comparison.png"
    fig.savefig(comparison, dpi=300, facecolor="white", metadata={
        "Title": "A4-C0 PHaRLAP-AH3D vacuum comparison",
        "Description": "No translation, rotation, scaling, or fitted registration applied",
    })
    plt.close(fig)
    paths["comparison"] = comparison

    manifest = {}
    for key, path in paths.items():
        image = mpimg.imread(path)
        checks = {
            "exists": path.is_file(),
            "finite_pixels": bool(np.all(np.isfinite(image))),
            "minimum_dimensions": bool(image.shape[0] >= 900 and image.shape[1] >= 1800),
            "nonblank_pixel_variance": bool(float(np.std(image[..., :3])) >= 1.0e-3),
            "source_hashes_recorded": all(
                isinstance(value, str) and len(value) == 64
                for value in (source_hashes.values() if key == "comparison"
                              else (source_hashes[key],))),
        }
        manifest[key] = {
            "path": str(path), "sha256": _sha256(path),
            "pixel_shape": list(image.shape), "pixel_std": float(np.std(image[..., :3])),
            "checks": checks,
            "source_sha256": source_hashes if key == "comparison" else {
                key: source_hashes[key]
            },
        }
    return manifest


def run() -> int:
    input_dir = OUTPUT / "inputs"
    raw_dir = OUTPUT / "raw"
    preflight_dir = OUTPUT / "preflight"
    for directory in (input_dir, raw_dir, preflight_dir):
        directory.mkdir(parents=True, exist_ok=True)
    input_path = input_dir / f"{CASE_ID}.mat"
    raw_pharlap = raw_dir / f"{CASE_ID}_pharlap.mat"
    raw_ah3d = raw_dir / f"{CASE_ID}_ah3d.json"
    preflight_mat = preflight_dir / "pharlap_preflight.mat"
    preflight_json = preflight_dir / "pharlap_preflight.json"
    log_path = OUTPUT / "run.log"
    report_path = OUTPUT / "a4_c0_validation.json"

    prerequisites = {
        "matlab_exists": MATLAB.is_file(),
        "pharlap_mex_exists": (PHARLAP / "mex" / "raytrace_3d_sp.mexw64").is_file(),
        "matlab_script_exists": (MATLAB_SCRIPTS / "a4_pharlap_preflight.m").is_file(),
    }
    if not all(prerequisites.values()):
        report = {"schema_version": "a4-c0-validation-v1", "case_id": CASE_ID,
                  "status": "blocked_pharlap_runtime", "prerequisites": prerequisites}
        _atomic_json(report_path, report)
        return 2

    input_metadata = _generate_input(input_path)
    process = _run_matlab(input_path, preflight_mat, preflight_json,
                          raw_pharlap, log_path)
    if process.returncode != 0 or not preflight_json.exists() or not raw_pharlap.exists():
        report = {
            "schema_version": "a4-c0-validation-v1", "case_id": CASE_ID,
            "status": "blocked_pharlap_runtime", "prerequisites": prerequisites,
            "matlab_returncode": process.returncode, "run_log": str(log_path),
        }
        _atomic_json(report_path, report)
        return 2

    preflight = json.loads(preflight_json.read_text(encoding="utf-8"))
    ah = _run_ah3d(raw_ah3d)
    metrics, checks, plot_data = _compare(raw_pharlap, ah)
    figures = _plot_c0(plot_data, OUTPUT / "figures", {
        "pharlap": _sha256(raw_pharlap), "ah3d": _sha256(raw_ah3d),
    })
    checks.update({
        "preflight_passed": preflight.get("status") == "passed",
        "input_shape_and_units_validated": bool(preflight.get("input_validated")),
        "mex_precedes_help_file": str(preflight.get("raytrace_3d_sp_path", "")).lower().endswith(".mexw64"),
        "figure_evidence_complete": all(
            all(item["checks"].values()) for item in figures.values()),
    })
    report = {
        "schema_version": "a4-c0-validation-v1",
        "case_id": CASE_ID,
        "case_semantics": "interface_and_coordinate_precheck_not_physics_validation",
        "status": "passed" if all(checks.values()) else "failed",
        "external_validation_claim": "none_C0_only",
        "prerequisites": prerequisites,
        "input": {
            **input_metadata, "path": str(input_path), "sha256": _sha256(input_path),
            "units": "Ne=electrons/cm^3,B=T,collision=Hz,distance=km",
        },
        "software": {
            "matlab": str(MATLAB), "pharlap_root": str(PHARLAP),
            "project_matlab_script": str(MATLAB_SCRIPTS / "a4_pharlap_preflight.m"),
            "pharlap_mex_sha256": _sha256(PHARLAP / "mex" / "raytrace_3d_sp.mexw64"),
            "matlab_script_sha256": _sha256(MATLAB_SCRIPTS / "a4_pharlap_preflight.m"),
            "preflight": preflight,
        },
        "configuration": {
            "frequency_MHz": FREQUENCY_MHZ, "elevation_deg": ELEVATION_DEG,
            "bearing_deg": BEARING_DEG, "earth_radius_m": R_KM * 1000.0,
            "OX_mode": 0, "collision_frequency_Hz": 0.0,
            "geomagnetic_field_T": 0.0, "NRT_NUM_THREADS": 1,
            "pharlap_tol": [1.0e-8, 0.01, 10.0],
            "ah3d_profile": "reference",
        },
        "metrics": metrics,
        "checks": checks,
        "figures": figures,
        "artifacts": {
            "preflight_mat": str(preflight_mat), "preflight_json": str(preflight_json),
            "pharlap_raw": str(raw_pharlap), "ah3d_raw": str(raw_ah3d),
            "run_log": str(log_path),
        },
        "limitations": [
            "C0 validates runtime, units, axes, coordinates, and status mapping only",
            "C0 does not validate ionospheric reflection, O/X splitting, P2P roots, scattering, or power",
            "AH3D earth-flattening and PHaRLAP spherical equations remain distinct",
            "AH3D vacuum group/local-geometric path mismatch is diagnostic here and must be tested by the C1 formal group-delay gate",
        ],
    }
    _atomic_json(report_path, report)
    print(f"A4-C0 status={report['status']} output={report_path}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(run())
