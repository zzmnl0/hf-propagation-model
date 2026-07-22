"""C0R-7: direct PHaRLAP/AH3D ECEF vacuum comparison."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.ne_provider_geocentric import VacuumGeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    COORDINATE_MODEL,
    RefractiveIndexAH3DECEF,
    integrate_fixed_path_ah3d_ecef,
)
from utils import enu_basis, geocentric_to_ecef


R_KM = 6371.0
TX_LAT_DEG = 30.0
TX_LON_DEG = 120.0
ELEVATION_DEG = 35.0
BEARING_DEG = 37.0
FREQUENCY_MHZ = 10.0
COMPARE_PATH_KM = 250.0
OUT = ROOT / "output2" / "a4_pharlap"
INPUT = OUT / "inputs" / "a4_c0_vacuum.mat"
PH_RAW = OUT / "raw" / "a4_c0_vacuum_pharlap.mat"
C0_REPORT = OUT / "a4_c0_validation.json"
AH_RAW = OUT / "raw" / "a4_c0r7_vacuum_ah3d_ecef.json"
REPORT = OUT / "a4_c0r7_validation.json"
FIGURES = OUT / "figures_c0r7"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def _resample(points: np.ndarray, distance: np.ndarray,
              sample: np.ndarray) -> np.ndarray:
    keep = np.r_[True, np.diff(distance) > 1.0e-12]
    return np.column_stack([
        np.interp(sample, distance[keep], points[keep, axis])
        for axis in range(3)
    ])


def _validate_frozen_pharlap() -> tuple[dict, dict]:
    if not all(path.is_file() for path in (INPUT, PH_RAW, C0_REPORT)):
        raise RuntimeError("frozen C0 PHaRLAP input/raw/report is incomplete")
    c0 = json.loads(C0_REPORT.read_text(encoding="utf-8"))
    input_hash = _sha256(INPUT)
    ph_hash = _sha256(PH_RAW)
    if c0.get("status") != "passed":
        raise RuntimeError("frozen C0 report did not pass")
    if input_hash != c0["input"]["sha256"]:
        raise RuntimeError("frozen PHaRLAP input hash drift; rerun MATLAB explicitly")
    bound_ph_hash = c0["figures"]["pharlap"]["source_sha256"]["pharlap"]
    if ph_hash != bound_ph_hash:
        raise RuntimeError("frozen PHaRLAP raw MAT hash drift; rerun MATLAB explicitly")
    expected = {
        "frequency_MHz": FREQUENCY_MHZ,
        "OX_mode": 0.0,
        "earth_radius_m": R_KM * 1000.0,
        "origin_lat_deg": TX_LAT_DEG,
        "origin_lon_deg": TX_LON_DEG,
        "elevation_deg": ELEVATION_DEG,
        "bearing_deg": BEARING_DEG,
    }
    inp = loadmat(INPUT, simplify_cells=True)
    if any(float(inp[key]) != value for key, value in expected.items()):
        raise RuntimeError("frozen PHaRLAP scalar configuration does not match C0R-7")
    for key in ("iono_en_grid_cm3", "iono_en_grid_5_cm3",
                "collision_freq_Hz", "Bx_T", "By_T", "Bz_T"):
        if np.any(np.asarray(inp[key], dtype=float) != 0.0):
            raise RuntimeError(f"C0R-7 requires zero field: {key}")
    return c0, {"input": input_hash, "pharlap_raw": ph_hash}


def _run_ah3d() -> dict:
    provider = VacuumGeocentricNeProvider3D(R_KM)
    model = RefractiveIndexAH3DECEF(
        provider, earth_radius_km=R_KM, freq_MHz=FREQUENCY_MHZ,
        wave_mode="O", gyro_MHz=0.0, b_ecef_hat=(0.0, 0.0, 1.0))
    result = integrate_fixed_path_ah3d_ecef(
        model, tx_lat_deg=TX_LAT_DEG, tx_lon_deg=TX_LON_DEG,
        tx_height_km=0.0, elevation_deg=ELEVATION_DEG,
        bearing_deg=BEARING_DEG, path_length_km=COMPARE_PATH_KM,
        profile="reference")
    raw = {
        "case_semantics": "one_way_vacuum_strict_ecef_coordinate_gate",
        "coordinate_model": result["coordinate_model"],
        "earth_radius_km": R_KM,
        "frequency_MHz": FREQUENCY_MHZ,
        "wave_mode_semantics": "B0_magnetoionic_degenerate_single_count",
        "termination_status": result["termination_status"],
        "geom_path_km": result["geom_path_km"],
        "group_path_km": result["group_path_km"],
        "hamiltonian_residual_max": result["hamiltonian_residual_max"],
        "provider_hash": result["provider_hash"],
        "coordinate_mapping_hash": result["coordinate_mapping_hash"],
        "trajectory_geom_path_km": np.asarray(
            result["trajectory_geom_path_km"], dtype=float).tolist(),
        "trajectory_ecef_km": np.asarray(
            result["trajectory_ecef_km"], dtype=float).tolist(),
    }
    _json_write(AH_RAW, raw)
    return raw


def _load_and_compare(ah: dict) -> tuple[dict, dict]:
    raw = loadmat(PH_RAW, simplify_cells=True)
    path = raw["ray_path_data"]
    ph_ecef = geocentric_to_ecef(
        path["lat"], path["lon"], path["height"], R_KM)
    ph_s = np.asarray(path["geometric_distance"], dtype=float).reshape(-1)
    ph_s -= ph_s[0]
    ah_ecef = np.asarray(ah["trajectory_ecef_km"], dtype=float)
    ah_s = np.asarray(ah["trajectory_geom_path_km"], dtype=float)
    common = min(COMPARE_PATH_KM, float(ph_s[-1]), float(ah_s[-1]))
    sample = np.linspace(0.0, common, 201)
    ph_sample = _resample(ph_ecef, ph_s, sample)
    ah_sample = _resample(ah_ecef, ah_s, sample)
    error = np.linalg.norm(ah_sample - ph_sample, axis=1)
    ph_group = float(np.asarray(path["group_range"]).reshape(-1)[-1])
    ph_geom = float(ph_s[-1])
    ah_group = float(ah["group_path_km"])
    ah_geom = float(ah["geom_path_km"])
    metrics = {
        "comparison_path_km": common,
        "trajectory_rms_3d_km": float(np.sqrt(np.mean(error * error))),
        "trajectory_max_3d_km": float(np.max(error)),
        "origin_error_km": float(np.linalg.norm(ph_sample[0] - ah_sample[0])),
        "pharlap_group_geom_relative_difference": abs(ph_group - ph_geom) / ph_geom,
        "ah3d_group_geom_relative_difference": abs(ah_group - ah_geom) / ah_geom,
        "ah3d_hamiltonian_residual_max": float(ah["hamiltonian_residual_max"]),
    }
    data = {
        "ph_ecef": ph_ecef, "ph_s": ph_s,
        "ah_ecef": ah_ecef, "ah_s": ah_s,
        "sample": sample, "ph_sample": ph_sample, "ah_sample": ah_sample,
    }
    return metrics, data


def _enu(points: np.ndarray) -> np.ndarray:
    origin = geocentric_to_ecef(TX_LAT_DEG, TX_LON_DEG, 0.0, R_KM)
    return (points - origin) @ enu_basis(TX_LAT_DEG, TX_LON_DEG, R_KM)


def _style(ax) -> None:
    ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)


def _plot_solver(name: str, ecef: np.ndarray, distance: np.ndarray,
                 color: str, source_hash: str) -> Path:
    enu = _enu(ecef)
    height = np.linalg.norm(ecef, axis=1) - R_KM
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    axes[0].plot(distance, height, color=color, linewidth=1.6)
    axes[0].set(xlabel="Geometric path (km)", ylabel="Physical height (km)")
    axes[1].plot(enu[:, 0], enu[:, 1], color=color, linewidth=1.6)
    axes[1].set(xlabel="East from TX (km)", ylabel="North from TX (km)")
    axes[1].set_aspect("equal", adjustable="box")
    for ax in axes:
        _style(ax)
    fig.suptitle(
        f"A4-C0R-7 {name} | one-way vacuum | strict spherical/ECEF\n"
        "10 MHz | elevation 35 deg | bearing 37 deg | Ne=0 | B=0 | collision=0",
        fontsize=9.5)
    path = FIGURES / f"trajectory_c0r7_{name.lower()}.png"
    fig.savefig(path, dpi=300, facecolor="white", metadata={
        "Title": f"A4-C0R-7 {name} vacuum trajectory",
        "Description": f"source_sha256={source_hash}; no registration applied",
    })
    plt.close(fig)
    return path


def _plot_all(data: dict, hashes: dict) -> dict:
    FIGURES.mkdir(parents=True, exist_ok=True)
    paths = {
        "pharlap": _plot_solver("PHaRLAP", data["ph_ecef"], data["ph_s"],
                                "#0072B2", hashes["pharlap_raw"]),
        "ah3d": _plot_solver("AH3D", data["ah_ecef"], data["ah_s"],
                             "#D55E00", hashes["ah3d_raw"]),
    }
    sample = data["sample"]
    ph = data["ph_sample"]
    ah = data["ah_sample"]
    ph_enu, ah_enu = _enu(ph), _enu(ah)
    ph_h = np.linalg.norm(ph, axis=1) - R_KM
    ah_h = np.linalg.norm(ah, axis=1) - R_KM
    error = np.linalg.norm(ah - ph, axis=1)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)
    axes[0].plot(sample, ph_h, color="#0072B2", label="PHaRLAP")
    axes[0].plot(sample, ah_h, "--", color="#D55E00", label="AH3D")
    axes[0].set(xlabel="Common geometric path (km)", ylabel="Physical height (km)")
    axes[1].plot(ph_enu[:, 0], ph_enu[:, 1], color="#0072B2", label="PHaRLAP")
    axes[1].plot(ah_enu[:, 0], ah_enu[:, 1], "--", color="#D55E00", label="AH3D")
    axes[1].set(xlabel="East from TX (km)", ylabel="North from TX (km)")
    axes[1].set_aspect("equal", adjustable="box")
    axes[2].plot(sample, error, color="#000000")
    axes[2].set(xlabel="Common geometric path (km)", ylabel="Direct ECEF error (km)")
    for ax in axes:
        _style(ax)
        ax.legend(frameon=False, fontsize=8) if ax is not axes[2] else None
    fig.suptitle(
        "A4-C0R-7 PHaRLAP vs AH3D | direct common-ECEF comparison | no registration\n"
        "one-way vacuum | spherical Earth R=6371 km | first 250 km",
        fontsize=9.5)
    paths["comparison"] = FIGURES / "trajectory_c0r7_comparison.png"
    fig.savefig(paths["comparison"], dpi=300, facecolor="white", metadata={
        "Title": "A4-C0R-7 direct PHaRLAP-AH3D vacuum comparison",
        "Description": "No translation, rotation, scaling, fitting, or registration",
    })
    plt.close(fig)
    manifest = {}
    for name, path in paths.items():
        pixels = mpimg.imread(path)
        manifest[name] = {
            "path": str(path), "sha256": _sha256(path),
            "pixel_shape": list(pixels.shape),
            "pixel_std": float(np.std(pixels[..., :3])),
            "checks": {
                "exists": path.is_file(),
                "finite_pixels": bool(np.all(np.isfinite(pixels))),
                "nonblank": float(np.std(pixels[..., :3])) >= 1.0e-3,
                "minimum_size": pixels.shape[0] >= 900 and pixels.shape[1] >= 1800,
            },
        }
    return manifest


def run() -> int:
    try:
        c0, hashes = _validate_frozen_pharlap()
        ah = _run_ah3d()
        hashes["ah3d_raw"] = _sha256(AH_RAW)
        metrics, plot_data = _load_and_compare(ah)
        figures = _plot_all(plot_data, hashes)
        checks = {
            "frozen_pharlap_input_and_raw_hashes_match": True,
            "pharlap_runtime_preflight_passed": c0["software"]["preflight"]["status"] == "passed",
            "strict_ecef_ah3d_used": ah["coordinate_model"] == COORDINATE_MODEL,
            "ah3d_fixed_path_completed": ah["termination_status"] == "fixed_geometric_path_completed",
            "comparison_interval_complete": metrics["comparison_path_km"] >= 250.0 - 1.0e-9,
            "origin_matches_without_registration": metrics["origin_error_km"] <= 1.0e-9,
            "trajectory_rms_gate": metrics["trajectory_rms_3d_km"] <= 1.0e-3,
            "trajectory_max_gate": metrics["trajectory_max_3d_km"] <= 1.0e-2,
            "pharlap_group_geom_gate": metrics["pharlap_group_geom_relative_difference"] <= 1.0e-4,
            "ah3d_group_geom_gate": metrics["ah3d_group_geom_relative_difference"] <= 1.0e-4,
            "ah3d_hamiltonian_gate": metrics["ah3d_hamiltonian_residual_max"] <= 1.0e-8,
            "figure_evidence_complete": all(
                all(item["checks"].values()) for item in figures.values()),
        }
        report = {
            "schema_version": "a4-c0r7-pharlap-vacuum-v1",
            "stage": "A4-C0R-7",
            "status": "passed" if all(checks.values()) else "failed",
            "case_semantics": "vacuum_coordinate_and_path_metric_gate_only",
            "external_validation_claim": "PHaRLAP_vacuum_coordinate_contract_only",
            "configuration": {
                "Ne_m^-3": 0.0, "B_T": 0.0, "collision_Hz": 0.0,
                "frequency_MHz": FREQUENCY_MHZ,
                "elevation_deg": ELEVATION_DEG, "bearing_deg": BEARING_DEG,
                "earth_radius_km": R_KM, "comparison_path_km": COMPARE_PATH_KM,
                "ah3d_integrator": "DOP853_reference",
                "coordinate_model": COORDINATE_MODEL,
                "registration": "none",
            },
            "metrics": metrics, "checks": checks, "figures": figures,
            "source_sha256": {
                **hashes,
                "test": _sha256(Path(__file__)),
                "ecef_ah3d": _sha256(ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
                "provider": _sha256(ROOT / "models" / "ne_provider_geocentric.py"),
                "geometry": _sha256(ROOT / "utils.py"),
            },
            "claim_boundary": (
                "C0R-7 validates the one-way vacuum spherical/ECEF coordinate and "
                "path-metric contract only; it does not validate ionospheric "
                "reflection, O/X splitting, P2P mode completeness, PE, scattering, "
                "or absolute power."
            ),
        }
    except Exception as exc:
        report = {
            "schema_version": "a4-c0r7-pharlap-vacuum-v1",
            "stage": "A4-C0R-7", "status": "blocked",
            "case_semantics": "vacuum_coordinate_and_path_metric_gate_only",
            "error": f"{type(exc).__name__}: {exc}",
        }
        checks = {"execution_completed": False}
    _json_write(REPORT, report)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-7: {failures} fail(s); status={report['status']}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(run())
