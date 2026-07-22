"""Run and freeze Prep-5 1-km branch convergence and complete signatures."""

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

from models.point_to_point_ah_3d_ecef import _trajectory_rms_km
from tests._run_a4_c1_prep3 import _load_context, _solve_roots
from tests._run_a4_c1_prep4 import _turning_count
from utils import enu_basis, geocentric_to_ecef


BASE = ROOT / "output2" / "a4_pharlap" / "c1_prep"
PREP1 = BASE / "prep1"
PREP4 = BASE / "prep4"
OUT = BASE / "prep5"
CHECKPOINT = OUT / "checkpoints"
PREP1_REPORT = PREP1 / "a4_c1_prep1_validation.json"
PREP1_MANIFEST = PREP1 / "a4_c1_prep1_freeze_manifest.json"
PREP4_REPORT = PREP4 / "a4_c1_prep4_validation.json"
PREP4_MANIFEST = PREP4 / "a4_c1_prep4_freeze_manifest.json"
PREP4_NPZ = PREP4 / "a4_c1_prep4_branch_trajectories.npz"

MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SOURCE = ROOT / "tests" / "matlab" / "a4_c1_prep5_pharlap.m"
MATLAB_SCRIPT = MATLAB_DIR / "a4_c1_prep5_pharlap.m"
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
PH_RAW = CHECKPOINT / "pharlap_refined_1km_raw.mat"
PH_JSON = CHECKPOINT / "pharlap_refined_1km_summary.json"
PH_LOG = CHECKPOINT / "pharlap_refined_1km.log"
PH_CHECKPOINT = CHECKPOINT / "pharlap_checkpoint.json"

CONVERGENCE_CSV = OUT / "a4_c1_prep5_convergence.csv"
SIGNATURE_CSV = OUT / "a4_c1_prep5_frozen_signatures.csv"
TRAJECTORY_NPZ = OUT / "a4_c1_prep5_trajectories.npz"
FIGURE = OUT / "a4_c1_prep5_grid_convergence.png"
REPORT = OUT / "a4_c1_prep5_validation.json"
MANIFEST = OUT / "a4_c1_prep5_freeze_manifest.json"

AH_SOLVER = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
AH_PROVIDER = ROOT / "models" / "ne_provider_geocentric.py"
AH_PREP3 = ROOT / "tests" / "_run_a4_c1_prep3.py"
CONFIG = ROOT / "config.py"

LIMITS = {
    "elevation_abs_deg": 0.02,
    "apex_abs_km": 0.5,
    "tau_relative_fraction": 0.002,
    "endpoint_residual_km": 1.0,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _array_hash(value) -> str:
    array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
    digest = hashlib.sha256(str(array.shape).encode("ascii"))
    digest.update(array.tobytes())
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _verify_manifest(path: Path, status: str) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != status:
        raise RuntimeError(f"unexpected prerequisite status: {path}")
    for item in manifest.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"frozen prerequisite changed: {artifact}")
    return manifest


def _source_hashes(refined_hash: str, prep4_hash: str) -> dict:
    return {
        "python": _sha256(Path(__file__)),
        "refined_input": refined_hash,
        "prep4_report": prep4_hash,
        "prep4_manifest": _sha256(PREP4_MANIFEST),
        "matlab_source": _sha256(MATLAB_SOURCE),
        "matlab_deployed": _sha256(MATLAB_SCRIPT),
        "pharlap_mex": _sha256(MEX),
        "ah_solver": _sha256(AH_SOLVER),
        "ah_provider": _sha256(AH_PROVIDER),
        "ah_prep3_driver": _sha256(AH_PREP3),
        "config": _sha256(CONFIG),
    }


def _verify_existing(allow_rebuild: bool, sources: dict | None = None) -> int | None:
    if not MANIFEST.is_file():
        return None
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    valid = (manifest.get("status") == "passed_1km_grid_stability_and_signatures"
             and sources is not None and manifest.get("source_sha256") == sources)
    for item in manifest.get("frozen_artifacts", []):
        path = Path(item["path"])
        valid = valid and path.is_file() and _sha256(path) == item["sha256"]
    if valid:
        print("A4-C1-Prep-5 frozen: OK")
        return 0
    if allow_rebuild:
        return None
    print("A4-C1-Prep-5 frozen/source hash mismatch; rebuild refused")
    return 1


def _pharlap_brackets(prep4: dict) -> list[list[float]]:
    items = sorted((item for item in prep4["signatures"]
                    if item["solver"] == "PHaRLAP"),
                   key=lambda item: item["elevation_deg"])
    return [[float(item["bracket_audit"]["bracket_low_deg"]),
             float(item["bracket_audit"]["bracket_high_deg"])] for item in items]


def _ah_brackets(prep4: dict) -> list[tuple[str, tuple[float, float]]]:
    items = sorted((item for item in prep4["signatures"]
                    if item["solver"] == "AH3D"),
                   key=lambda item: item["elevation_deg"])
    return [(item["branch_label"], (
        float(item["bracket_audit"]["bracket_low_deg"]),
        float(item["bracket_audit"]["bracket_high_deg"]))) for item in items]


def _run_matlab(input_path: Path, brackets: list[list[float]]) -> subprocess.CompletedProcess:
    quote = lambda value: str(value).replace("'", "''")
    bracket_text = "[" + ";".join(f"{low:.12g},{high:.12g}"
                                    for low, high in brackets) + "]"
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); "
        f"a4_c1_prep5_pharlap('{quote(input_path)}','{quote(PH_RAW)}',"
        f"'{quote(PH_JSON)}',{bracket_text})")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run([str(MATLAB), "-batch", expression], env=env,
                            text=True, capture_output=True, timeout=7200,
                            errors="replace")
    PH_LOG.write_text((result.stdout or "") + "\n" + (result.stderr or ""),
                      encoding="utf-8")
    return result


def _load_pharlap_paths(radius_km: float) -> list[dict]:
    paths = []
    with h5py.File(PH_RAW, "r") as source:
        references = source["root_raw/ray_path_data"][()].reshape(-1)
        for reference in references:
            group = source[reference]
            latitude = np.asarray(group["lat"], dtype=float).reshape(-1)
            longitude = np.asarray(group["lon"], dtype=float).reshape(-1)
            height = np.asarray(group["height"], dtype=float).reshape(-1)
            path = np.asarray(group["geometric_distance"], dtype=float).reshape(-1)
            path -= path[0]
            paths.append({
                "trajectory_ecef_km": geocentric_to_ecef(
                    latitude, longitude, height, radius_km),
                "trajectory_geom_path_km": path,
                "trajectory_geocentric": np.column_stack((latitude, longitude, height)),
                "turning_count": _turning_count(height),
            })
    return paths


def _ah_checkpoint_paths(label: str) -> tuple[Path, Path]:
    return (CHECKPOINT / f"ah3d_{label}.json", CHECKPOINT / f"ah3d_{label}.npz")


def _load_ah_checkpoint(label: str, sources: dict, allow_rebuild: bool):
    metadata_path, data_path = _ah_checkpoint_paths(label)
    if not metadata_path.is_file() and not data_path.is_file():
        return None
    if not metadata_path.is_file() or not data_path.is_file():
        if allow_rebuild:
            return None
        raise RuntimeError(f"incomplete AH3D checkpoint: {label}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (metadata.get("source_sha256") != sources
            or metadata.get("trajectory_npz_sha256") != _sha256(data_path)):
        if allow_rebuild:
            return None
        raise RuntimeError(f"AH3D checkpoint hash mismatch: {label}")
    with np.load(data_path, allow_pickle=False) as data:
        trajectory = {key: np.asarray(data[key], dtype=float) for key in data.files}
    return metadata["root"], trajectory, metadata["trace_hashes"]


def _write_ah_checkpoint(label: str, root: dict, trajectory: dict,
                         trace_hashes: dict, sources: dict) -> None:
    metadata_path, data_path = _ah_checkpoint_paths(label)
    np.savez_compressed(data_path, **trajectory)
    _atomic_json(metadata_path, {
        "schema_version": "a4-c1-prep5-ah3d-checkpoint-v1",
        "branch_label": label,
        "source_sha256": sources,
        "trajectory_npz_sha256": _sha256(data_path),
        "trace_hashes": trace_hashes,
        "root": root,
    })


def _solve_ah_worker(payload):
    label, input_path, bracket = payload
    context = _load_context(input_path)
    roots, rays = _solve_roots(context, [tuple(bracket)])
    if len(roots) != 1 or len(rays) != 1:
        raise RuntimeError(f"AH3D refined root failed: {label}")
    root, ray = roots[0], rays[0]
    root["branch_label"] = label
    trajectory = {
        "trajectory_ecef_km": np.asarray(ray["trajectory_ecef_km"], dtype=float),
        "trajectory_geom_path_km": np.asarray(ray["trajectory_geom_path_km"], dtype=float),
        "trajectory_geocentric": np.asarray(
            ray["trajectory_geocentric_lat_lon_height"], dtype=float),
    }
    trace_hashes = {
        "source_provider_hash": context["source_provider_hash"],
        "trace_provider_hash": context["trace_provider_hash"],
        "coordinate_mapping_hash": ray["coordinate_mapping_hash"],
    }
    return label, root, trajectory, trace_hashes


def _load_baseline(prep4: dict) -> dict:
    signatures = {(item["solver"], item["branch_label"]): item
                  for item in prep4["signatures"]}
    result = {}
    with np.load(PREP4_NPZ, allow_pickle=False) as data:
        for solver in ("PHaRLAP", "AH3D"):
            for label in ("low_ray", "high_ray"):
                prefix = f"{solver.lower()}_{label}"
                result[(solver, label)] = {
                    "root": signatures[(solver, label)],
                    "trajectory": {
                        "trajectory_ecef_km": np.asarray(
                            data[f"{prefix}_ecef_km"], dtype=float),
                        "trajectory_geom_path_km": np.asarray(
                            data[f"{prefix}_geom_path_km"], dtype=float),
                    },
                }
    return result


def _complete_signature(solver: str, label: str, root: dict, trajectory: dict,
                        input_record: dict, trace_hashes: dict, data,
                        source_hashes: dict) -> dict:
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    rx = np.array([float(np.asarray(data[name]).item()) for name in
                   ("receiver_lat_deg", "receiver_lon_deg", "receiver_height_km")])
    terminal = np.asarray(trajectory["trajectory_ecef_km"], dtype=float)[-1]
    receiver = geocentric_to_ecef(*rx, radius)
    residual_enu = (terminal - receiver) @ enu_basis(rx[0], rx[1], radius)
    terminal_geocentric = np.asarray(
        trajectory.get("trajectory_geocentric", [[root["landing_lat_deg"],
                                                    root["landing_lon_deg"], 0.0]]),
        dtype=float)[-1]
    return {
        "solver": solver,
        "branch_id": label,
        "wave_branch_id": "magnetoionic_degenerate_single_count",
        "grid_variant": "refined_1km",
        "elevation_deg": float(root["elevation_deg"]),
        "bearing_deg": float(np.asarray(data["bearing_deg"]).item()),
        "termination": str(root["status"]),
        "turning_count": int(root.get("turning_count",
            _turning_count(trajectory["trajectory_geocentric"][:, 2]))),
        "apex_km": float(root["apex_km"]),
        "landing_geocentric_lat_lon_height": terminal_geocentric.tolist(),
        "landing_ecef_km": terminal.tolist(),
        "landing_rx_enu_residual_km": residual_enu.tolist(),
        "along_track_residual_km": float(root["along_track_residual_km"]),
        "cross_track_km": float(root["cross_track_km"]),
        "trajectory_ecef_sha256": _array_hash(trajectory["trajectory_ecef_km"]),
        "trajectory_point_count": int(len(trajectory["trajectory_ecef_km"])),
        "geometric_path_km": float(root["geometric_path_km"]),
        "group_path_km": float(root["group_path_km"]),
        "tau_ms": float(root["tau_ms"]),
        "provider_hashes": trace_hashes,
        "input_path": input_record["path"],
        "input_sha256": input_record["sha256"],
        "code_sha256": source_hashes,
    }


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _range_height(trajectory: dict, radius: float, tx_ecef: np.ndarray):
    ecef = np.asarray(trajectory["trajectory_ecef_km"], dtype=float)
    unit = ecef / np.linalg.norm(ecef, axis=1)[:, None]
    tx_unit = tx_ecef / np.linalg.norm(tx_ecef)
    return (radius * np.arccos(np.clip(unit @ tx_unit, -1.0, 1.0)),
            np.linalg.norm(ecef, axis=1) - radius)


def _plot(rows: list[dict], refined: dict, baseline: dict, data,
          source_hashes: dict) -> dict:
    labels = [f"{row['solver']}\n{row['branch_id']}" for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    metrics = (
        ("elevation_abs_difference_deg", LIMITS["elevation_abs_deg"],
         "|1 km - 2 km| elevation (deg)"),
        ("apex_abs_difference_km", LIMITS["apex_abs_km"],
         "|1 km - 2 km| apex (km)"),
        ("tau_relative_difference_percent", LIMITS["tau_relative_fraction"] * 100.0,
         "Relative group-delay change (%)"),
    )
    colors = ["#0072B2", "#56B4E9", "#D55E00", "#E69F00"]
    for axis, (key, threshold, title) in zip(axes.flat[:3], metrics):
        values = [row[key] for row in rows]
        axis.bar(x, values, color=colors)
        axis.axhline(threshold, color="black", linestyle="--", label="acceptance limit")
        axis.set_xticks(x, labels, rotation=20, ha="right")
        axis.set_ylabel(title)
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize=8)

    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    tx = np.array([float(np.asarray(data[name]).item()) for name in
                   ("origin_lat_deg", "origin_lon_deg", "origin_height_km")])
    tx_ecef = geocentric_to_ecef(*tx, radius)
    target = float(np.asarray(data["target_range_km"]).item())
    height_grid = np.asarray(data["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(data["iono_en_grid_m3"], dtype=float)
    profile = ne[ne.shape[0] // 2, ne.shape[1] // 2, :] / 1.0e11
    ranges = np.linspace(0.0, target * 1.12, 350)
    mesh = axes[1, 1].pcolormesh(
        ranges, height_grid, np.broadcast_to(profile[:, None],
        (height_grid.size, ranges.size)), cmap="jet", shading="auto", alpha=0.5)
    style = {("PHaRLAP", "low_ray"): "#0072B2",
             ("PHaRLAP", "high_ray"): "#56B4E9",
             ("AH3D", "low_ray"): "#D55E00",
             ("AH3D", "high_ray"): "#E69F00"}
    for key, current in refined.items():
        solver, label = key
        for trajectory, linestyle, grid in (
                (baseline[key]["trajectory"], "--", "2 km"),
                (current["trajectory"], "-", "1 km")):
            ray_range, ray_height = _range_height(trajectory, radius, tx_ecef)
            axes[1, 1].plot(ray_range, ray_height, color=style[key],
                            linestyle=linestyle, linewidth=1.4,
                            label=f"{solver} {label} {grid}")
    axes[1, 1].axvline(target, color="black", linestyle=":", label="RX")
    axes[1, 1].set(xlabel="Great-circle range from TX (km)",
                   ylabel="Physical height (km)",
                   title="2 km dashed vs 1 km solid trajectories")
    axes[1, 1].set_ylim(0.0, 460.0)
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=6, ncol=2)
    fig.colorbar(mesh, ax=axes[1, 1], label="Ne (1e11 m^-3)")
    for label, axis in zip("ABCD", axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top",
                  fontsize=13, fontweight="bold")
    fig.suptitle(
        "A4-C1 Prep-5 | frozen Chapman branch grid-stability audit\n"
        "one-way | 10 MHz | B=0 O/X degenerate | 1 km vs 2 km | not A4-C1 validation",
        fontsize=12)
    fig.savefig(FIGURE, dpi=300, facecolor="white", metadata={
        "Title": "A4-C1 Prep-5 grid convergence",
        "Description": "Within-solver 1-km versus 2-km stability; not cross-solver acceptance",
    })
    plt.close(fig)
    image = mpimg.imread(FIGURE)
    checks = {
        "exists": FIGURE.is_file() and FIGURE.stat().st_size > 0,
        "finite_pixels": bool(np.all(np.isfinite(image))),
        "minimum_dimensions": bool(image.shape[0] >= 1800 and image.shape[1] >= 2500),
        "nonblank_pixel_variance": bool(float(np.std(image[..., :3])) >= 1.0e-3),
        "source_hashes_recorded": all(len(value) == 64 for value in source_hashes.values()),
    }
    return {"path": str(FIGURE), "sha256": _sha256(FIGURE),
            "source_sha256": source_hashes, "pixel_shape": list(image.shape),
            "pixel_std": float(np.std(image[..., :3])), "checks": checks}


def run(*, allow_rebuild: bool = False, workers: int = 2) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    CHECKPOINT.mkdir(parents=True, exist_ok=True)
    try:
        _verify_manifest(PREP1_MANIFEST, "passed_input_contract_frozen")
        _verify_manifest(PREP4_MANIFEST, "passed_branch_identification_and_joint_dedup")
        prep1 = json.loads(PREP1_REPORT.read_text(encoding="utf-8"))
        prep4 = json.loads(PREP4_REPORT.read_text(encoding="utf-8"))
        refined_input = next(item for item in prep1["inputs"]
                             if item["variant"] == "refined_1km")
        input_path = Path(refined_input["path"])
        if _sha256(input_path) != refined_input["sha256"]:
            raise RuntimeError("Prep-1 refined input hash mismatch")
        sources = _source_hashes(refined_input["sha256"], _sha256(PREP4_REPORT))
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_prerequisite", "error": str(error)})
        return 2
    existing = _verify_existing(allow_rebuild, sources)
    if existing is not None:
        return existing

    ph_sources = {key: sources[key] for key in (
        "refined_input", "prep4_report", "prep4_manifest", "matlab_source",
        "matlab_deployed", "pharlap_mex")}
    ph_reused = False
    if PH_CHECKPOINT.is_file() and PH_RAW.is_file() and PH_JSON.is_file():
        checkpoint = json.loads(PH_CHECKPOINT.read_text(encoding="utf-8"))
        valid = (checkpoint.get("source_sha256") == ph_sources
                 and checkpoint.get("raw_sha256") == _sha256(PH_RAW)
                 and checkpoint.get("summary_sha256") == _sha256(PH_JSON))
        if valid:
            ph_reused = True
        elif not allow_rebuild:
            _atomic_json(REPORT, {"status": "blocked_checkpoint_hash",
                                  "error": "PHaRLAP checkpoint hash mismatch"})
            return 2
    if not ph_reused:
        process = _run_matlab(input_path, _pharlap_brackets(prep4))
        summary = json.loads(PH_JSON.read_text(encoding="utf-8")) if PH_JSON.is_file() else {}
        if process.returncode != 0 or summary.get("status") != "passed_pharlap_1km_branch_resolve":
            _atomic_json(REPORT, {"status": "failed_pharlap_refined_roots",
                                  "returncode": process.returncode,
                                  "summary": summary})
            return 1
        _atomic_json(PH_CHECKPOINT, {
            "schema_version": "a4-c1-prep5-pharlap-checkpoint-v1",
            "source_sha256": ph_sources,
            "raw_sha256": _sha256(PH_RAW),
            "summary_sha256": _sha256(PH_JSON),
        })

    ah_sources = {key: sources[key] for key in (
        "python", "refined_input", "prep4_report", "prep4_manifest",
        "ah_solver", "ah_provider", "ah_prep3_driver", "config")}
    ah_results = {}
    missing = []
    try:
        for label, bracket in _ah_brackets(prep4):
            loaded = _load_ah_checkpoint(label, ah_sources, allow_rebuild)
            if loaded is None:
                missing.append((label, str(input_path), bracket))
            else:
                root, trajectory, trace_hashes = loaded
                ah_results[label] = {"root": root, "trajectory": trajectory,
                                     "trace_hashes": trace_hashes, "reused": True}
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_checkpoint_hash", "error": str(error)})
        return 2
    if missing:
        try:
            with ProcessPoolExecutor(max_workers=min(max(1, workers), len(missing))) as executor:
                futures = {executor.submit(_solve_ah_worker, item): item[0]
                           for item in missing}
                for future in as_completed(futures):
                    label, root, trajectory, trace_hashes = future.result()
                    _write_ah_checkpoint(
                        label, root, trajectory, trace_hashes, ah_sources)
                    ah_results[label] = {"root": root, "trajectory": trajectory,
                                         "trace_hashes": trace_hashes, "reused": False}
        except Exception as error:
            _atomic_json(REPORT, {"status": "failed_ah3d_refined_roots",
                                  "error": f"{type(error).__name__}: {error}"})
            return 1

    ph_summary = json.loads(PH_JSON.read_text(encoding="utf-8"))
    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    ph_paths = _load_pharlap_paths(radius)
    ph_roots = list(ph_summary["roots"])
    refined = {}
    for index, label in enumerate(("low_ray", "high_ray")):
        refined[("PHaRLAP", label)] = {
            "root": {**ph_roots[index], "branch_label": label},
            "trajectory": ph_paths[index],
            "trace_hashes": {"provider_hash": refined_input["provider_hash"],
                             "field_hash": refined_input["field_hash"]},
        }
        refined[("AH3D", label)] = ah_results[label]
    baseline = _load_baseline(prep4)

    convergence = []
    signatures = []
    trajectory_payload = {}
    for key in (("PHaRLAP", "low_ray"), ("PHaRLAP", "high_ray"),
                ("AH3D", "low_ray"), ("AH3D", "high_ray")):
        solver, label = key
        old = baseline[key]
        current = refined[key]
        root = current["root"]
        tau_relative = abs(float(root["tau_ms"]) - float(old["root"]["tau_ms"])) / abs(
            float(old["root"]["tau_ms"]))
        row = {
            "solver": solver,
            "branch_id": label,
            "elevation_2km_deg": float(old["root"]["elevation_deg"]),
            "elevation_1km_deg": float(root["elevation_deg"]),
            "elevation_abs_difference_deg": abs(
                float(root["elevation_deg"]) - float(old["root"]["elevation_deg"])),
            "apex_2km_km": float(old["root"]["apex_km"]),
            "apex_1km_km": float(root["apex_km"]),
            "apex_abs_difference_km": abs(
                float(root["apex_km"]) - float(old["root"]["apex_km"])),
            "tau_2km_ms": float(old["root"]["tau_ms"]),
            "tau_1km_ms": float(root["tau_ms"]),
            "tau_relative_difference": tau_relative,
            "tau_relative_difference_percent": tau_relative * 100.0,
            "trajectory_rms_km": _trajectory_rms_km(
                old["trajectory"], current["trajectory"]),
            "geometric_path_abs_difference_km": abs(
                float(root["geometric_path_km"])
                - float(old["root"]["geometric_path_km"])),
            "group_path_abs_difference_km": abs(
                float(root["group_path_km"])
                - float(old["root"]["group_path_km"])),
            "along_track_residual_km": float(root["along_track_residual_km"]),
            "cross_track_km": float(root["cross_track_km"]),
            "termination": str(root["status"]),
        }
        row["passed"] = bool(
            row["termination"] == "ground"
            and row["elevation_abs_difference_deg"] <= LIMITS["elevation_abs_deg"]
            and row["apex_abs_difference_km"] <= LIMITS["apex_abs_km"]
            and row["tau_relative_difference"] <= LIMITS["tau_relative_fraction"]
            and abs(row["along_track_residual_km"]) <= LIMITS["endpoint_residual_km"]
            and abs(row["cross_track_km"]) <= LIMITS["endpoint_residual_km"])
        convergence.append(row)
        signatures.append(_complete_signature(
            solver, label, root, current["trajectory"], refined_input,
            current["trace_hashes"], data, sources))
        prefix = f"{solver.lower()}_{label}"
        trajectory_payload[f"{prefix}_1km_ecef_km"] = current["trajectory"]["trajectory_ecef_km"]
        trajectory_payload[f"{prefix}_1km_geom_path_km"] = current["trajectory"]["trajectory_geom_path_km"]
        trajectory_payload[f"{prefix}_2km_ecef_km"] = old["trajectory"]["trajectory_ecef_km"]
        trajectory_payload[f"{prefix}_2km_geom_path_km"] = old["trajectory"]["trajectory_geom_path_km"]

    _write_csv(CONVERGENCE_CSV, convergence, list(convergence[0]))
    signature_rows = [{key: value for key, value in item.items()
                       if not isinstance(value, (dict, list))} for item in signatures]
    _write_csv(SIGNATURE_CSV, signature_rows, list(signature_rows[0]))
    np.savez_compressed(TRAJECTORY_NPZ, **trajectory_payload)
    plot_hashes = {
        "refined_input": refined_input["sha256"],
        "prep4_report": _sha256(PREP4_REPORT),
        "convergence_csv": _sha256(CONVERGENCE_CSV),
        "signature_csv": _sha256(SIGNATURE_CSV),
        "trajectory_npz": _sha256(TRAJECTORY_NPZ),
    }
    figure = _plot(convergence, refined, baseline, data, plot_hashes)

    checks = {
        "prep1_frozen": True,
        "prep4_frozen": True,
        "refined_input_hash": _sha256(input_path) == refined_input["sha256"],
        "matlab_source_matches_deployment": sources["matlab_source"] == sources["matlab_deployed"],
        "two_roots_per_solver": len(ph_roots) == len(ah_results) == 2,
        "all_roots_ground_and_endpoint_valid": all(row["termination"] == "ground"
            and abs(row["along_track_residual_km"]) <= 1.0
            and abs(row["cross_track_km"]) <= 1.0 for row in convergence),
        "all_elevation_changes_within_0p02deg": all(
            row["elevation_abs_difference_deg"] <= LIMITS["elevation_abs_deg"]
            for row in convergence),
        "all_apex_changes_within_0p5km": all(
            row["apex_abs_difference_km"] <= LIMITS["apex_abs_km"]
            for row in convergence),
        "all_tau_changes_within_0p2percent": all(
            row["tau_relative_difference"] <= LIMITS["tau_relative_fraction"]
            for row in convergence),
        "all_four_complete_signatures": len(signatures) == 4 and all(
            len(item["trajectory_ecef_sha256"]) == 64
            and len(item["landing_ecef_km"]) == 3
            and len(item["landing_rx_enu_residual_km"]) == 3
            and item["turning_count"] == 1 for item in signatures),
        "checkpoint_hash_contract_recorded": all(
            Path(path).is_file() for label in ("low_ray", "high_ray")
            for path in _ah_checkpoint_paths(label)) and PH_CHECKPOINT.is_file(),
        "prep4_labels_preserved": [item["branch_id"] for item in signatures]
            == ["low_ray", "high_ray", "low_ray", "high_ray"],
        "a4_c1_cross_solver_acceptance_not_applied": True,
        "visualization": all(figure["checks"].values()),
    }
    status = "passed_1km_grid_stability_and_signatures" if all(checks.values()) else "failed"
    report = {
        "schema_version": "a4-c1-prep5-validation-v1",
        "status": status,
        "case_semantics": (
            "within_solver_1km_vs_2km_grid_stability_and_complete_signatures_"
            "not_a4_c1_cross_solver_error_acceptance"),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "limits": LIMITS,
        "checkpoint_provenance": {
            "pharlap_reused": ph_reused,
            "ah3d_reused": {label: value["reused"] for label, value in ah_results.items()},
            "mismatch_semantics": "reject_unless_explicit_rebuild_frozen",
        },
        "refined_input": refined_input,
        "convergence": convergence,
        "frozen_signatures": signatures,
        "source_sha256": sources,
        "artifacts": {
            "pharlap_raw": str(PH_RAW),
            "pharlap_summary": str(PH_JSON),
            "convergence_csv": str(CONVERGENCE_CSV),
            "signature_csv": str(SIGNATURE_CSV),
            "trajectory_npz": str(TRAJECTORY_NPZ),
            "figure": figure,
        },
        "limitations": [
            "Prep-5 tests vertical-grid stability within each solver, not PHaRLAP-AH3D agreement",
            "Low/high labels remain topology-based rather than a direct optical-path Hessian test",
            "B=0 Chapman stability does not validate O/X splitting or complex 3-D mode completeness",
        ],
    }
    _atomic_json(REPORT, report)
    artifacts = [PH_RAW, PH_JSON, PH_LOG, PH_CHECKPOINT,
                 *[path for label in ("low_ray", "high_ray")
                   for path in _ah_checkpoint_paths(label)],
                 CONVERGENCE_CSV, SIGNATURE_CSV, TRAJECTORY_NPZ, FIGURE, REPORT]
    manifest = {
        "schema_version": "a4-c1-prep5-freeze-v1",
        "status": status,
        "case_semantics": report["case_semantics"],
        "source_sha256": sources,
        "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                             for path in artifacts],
    }
    _atomic_json(MANIFEST, manifest)
    print(f"A4-C1-Prep-5: {'OK' if status.startswith('passed') else 'FAIL'} "
          f"branches={len(convergence)} report={REPORT}")
    return 0 if status.startswith("passed") else 1


if __name__ == "__main__":
    arguments = sys.argv[1:]
    allow_rebuild = "--rebuild-frozen" in arguments
    workers = 2
    for argument in arguments:
        if argument.startswith("--workers="):
            workers = int(argument.split("=", 1)[1])
        elif argument != "--rebuild-frozen":
            raise SystemExit(
                "usage: _run_a4_c1_prep5.py [--rebuild-frozen] [--workers=N]")
    raise SystemExit(run(allow_rebuild=allow_rebuild, workers=workers))
