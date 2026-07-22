"""Run and freeze the independent AH3D full-elevation Chapman Prep-3 scan."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)
from utils import rx_enu_residual


PREP1 = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep1"
PREP1_REPORT = PREP1 / "a4_c1_prep1_validation.json"
PREP1_MANIFEST = PREP1 / "a4_c1_prep1_freeze_manifest.json"
OUT = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep3"
CSV = OUT / "a4_c1_prep3_scan.csv"
ROOT_CSV = OUT / "a4_c1_prep3_roots.csv"
ROOT_NPZ = OUT / "a4_c1_prep3_root_trajectories.npz"
REPORT = OUT / "a4_c1_prep3_validation.json"
MANIFEST = OUT / "a4_c1_prep3_freeze_manifest.json"
FIGURE = OUT / "a4_c1_prep3_ah3d_full_domain_scan.png"
SOLVER_SOURCE = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
PROVIDER_SOURCE = ROOT / "models" / "ne_provider_geocentric.py"

PROFILE = "research"
H_STOP_KM = 500.0
MAX_PATH_KM = 10000.0
HAMILTONIAN_LIMIT = 1.0e-8
ROOT_RESIDUAL_LIMIT_KM = 1.0
COARSE_ANGLES = np.arange(1.0, 89.0 + 0.25, 0.5)

_WORKER_CONTEXT = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _mat_text(value) -> str:
    return str(np.asarray(value).reshape(-1)[0])


def _verify_prep1() -> tuple[Path, str, dict]:
    manifest = json.loads(PREP1_MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "passed_input_contract_frozen":
        raise RuntimeError("Prep-1 is not frozen")
    for item in manifest["frozen_artifacts"]:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"Prep-1 frozen artifact changed: {path}")
    report = json.loads(PREP1_REPORT.read_text(encoding="utf-8"))
    main = next(item for item in report["inputs"] if item["variant"] == "main_2km")
    path = Path(main["path"])
    if _sha256(path) != main["sha256"]:
        raise RuntimeError("Prep-1 main input hash mismatch")
    return path, main["sha256"], main


def _load_context(input_path: str | Path) -> dict:
    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    latitude = np.asarray(data["latitude_deg"], dtype=float).reshape(-1)
    longitude = np.asarray(data["longitude_deg"], dtype=float).reshape(-1)
    source_height = np.asarray(data["height_km"], dtype=float).reshape(-1)
    source_ne = np.asarray(data["iono_en_grid_m3"], dtype=float)
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    source_provider = GeocentricGridNeProvider3D(
        source_ne, latitude, longitude, source_height, radius,
        source_metadata={"field_hash": _mat_text(data["field_hash"])},
    )
    if source_provider.grid_hash != _mat_text(data["provider_hash"]):
        raise RuntimeError("Prep-1 Provider hash cannot be reproduced")

    step = float(source_height[1] - source_height[0])
    # Chapman ghost nodes keep centered ECEF stencils inside their numerical
    # domain at the ground and 500-km events. Frozen 60--500 km nodes stay exact.
    lower_height = np.arange(-step, source_height[0], step, dtype=float)
    upper_height = source_height[-1] + step * np.arange(1.0, 3.0)
    foF2_MHz = float(np.asarray(data["foF2_MHz"]).item())
    hmF2_km = float(np.asarray(data["hmF2_km"]).item())
    ymF2_km = float(np.asarray(data["ymF2_km"]).item())
    continuation_height = np.r_[lower_height, upper_height]
    u = (continuation_height - hmF2_km) / ymF2_km
    continuation_profile = (foF2_MHz * 1.0e6 / 8.98) ** 2 * np.exp(np.clip(
        0.5 * (1.0 - u - np.exp(np.clip(-u, -700.0, 700.0))),
        -745.0, 700.0,
    ))
    lower_profile = continuation_profile[:lower_height.size]
    upper_profile = continuation_profile[lower_height.size:]
    height = np.r_[lower_height, source_height, upper_height]
    ne = np.concatenate((
        np.broadcast_to(
            lower_profile,
            (latitude.size, longitude.size, lower_height.size),
        ).copy(),
        source_ne,
        np.broadcast_to(
            upper_profile,
            (latitude.size, longitude.size, upper_height.size),
        ).copy(),
    ), axis=2)
    provider = GeocentricGridNeProvider3D(
        ne, latitude, longitude, height, radius,
        source_metadata={
            "source_input_sha256": _sha256(Path(input_path)),
            "source_grid_hash": source_provider.grid_hash,
            "boundary_continuation_semantics": (
                "same-Chapman 2-km continuation outside the frozen height domain; "
                "all frozen source nodes unchanged"
            ),
        },
    )
    tx = np.array([
        float(np.asarray(data["origin_lat_deg"]).item()),
        float(np.asarray(data["origin_lon_deg"]).item()),
        float(np.asarray(data["origin_height_km"]).item()),
    ])
    rx = np.array([
        float(np.asarray(data["receiver_lat_deg"]).item()),
        float(np.asarray(data["receiver_lon_deg"]).item()),
        float(np.asarray(data["receiver_height_km"]).item()),
    ])
    return {
        "provider": provider,
        "source_provider_hash": source_provider.grid_hash,
        "trace_provider_hash": provider.grid_hash,
        "source_height": source_height,
        "source_ne": source_ne,
        "source_start": int(lower_height.size),
        "height": height,
        "ne": ne,
        "latitude": latitude,
        "longitude": longitude,
        "radius": radius,
        "tx": tx,
        "rx": rx,
        "bearing_deg": float(np.asarray(data["bearing_deg"]).item()),
        "frequency_MHz": float(np.asarray(data["frequency_MHz"]).item()),
        "target_range_km": float(np.asarray(data["target_range_km"]).item()),
    }


def _model(context: dict) -> RefractiveIndexAH3DECEF:
    return RefractiveIndexAH3DECEF(
        context["provider"], earth_radius_km=context["radius"],
        freq_MHz=context["frequency_MHz"], wave_mode="O", gyro_MHz=0.0,
        b_ecef_hat=(0.0, 0.0, 1.0), collision_s=0.0,
        spatial_fd_base_km=0.5,
        spatial_fd_scales=(
            0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0),
        integration_profile=PROFILE,
    )


def _surface_distance_km(lat1: float, lon1: float, lat2: float, lon2: float,
                         radius_km: float) -> float:
    lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])
    dot = (np.sin(lat1) * np.sin(lat2)
           + np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1))
    return float(radius_km * np.arccos(np.clip(dot, -1.0, 1.0)))


def _trace(model: RefractiveIndexAH3DECEF, context: dict,
           elevation_deg: float) -> tuple[dict, dict]:
    ray = trace_ray_ah3d_ecef(
        model, tx_lat_deg=context["tx"][0], tx_lon_deg=context["tx"][1],
        tx_height_km=context["tx"][2], elevation_deg=float(elevation_deg),
        bearing_deg=context["bearing_deg"], h_stop_km=H_STOP_KM,
        max_path_km=MAX_PATH_KM, profile=PROFILE,
    )
    constraint_ok = (
        ray["failure_reason"] is None
        and np.isfinite(ray["hamiltonian_residual_max"])
        and ray["hamiltonian_residual_max"] <= HAMILTONIAN_LIMIT
        and ray["spatial_fd_rejection_count"] == 0
        and ray["spatial_fd_stability_ratio_max"] <= 1.0
    )
    landing = ray["landing_geocentric_lat_lon_height"]
    along = cross = ground_range = None
    if ray["termination_status"] == "ground" and landing is not None:
        terminal = np.asarray(ray["trajectory_ecef_km"], dtype=float)[-1]
        east, north = rx_enu_residual(
            terminal, context["rx"][0], context["rx"][1], context["rx"][2],
            context["radius"],
        )
        bearing = np.radians(context["bearing_deg"])
        along = float(east * np.sin(bearing) + north * np.cos(bearing))
        cross = float(east * np.cos(bearing) - north * np.sin(bearing))
        ground_range = _surface_distance_km(
            context["tx"][0], context["tx"][1], landing[0], landing[1],
            context["radius"],
        )
    record = {
        "elevation_deg": float(elevation_deg),
        "status": ray["termination_status"],
        "constraint_ok": bool(constraint_ok),
        "ground_range_km": ground_range,
        "along_track_residual_km": along,
        "cross_track_km": cross,
        "apex_km": float(np.max(ray["trajectory_physical_height_km"])),
        "turning_count": int(ray["turning_point_count"]),
        "geometric_path_km": float(ray["geom_path_km"]),
        "group_path_km": float(ray["group_path_km"]),
        "tau_ms": float(ray["tau_ms"]),
        "landing_lat_deg": None if landing is None else float(landing[0]),
        "landing_lon_deg": None if landing is None else float(landing[1]),
        "path_point_count": int(len(ray["trajectory_ecef_km"])),
        "hamiltonian_residual_max": float(ray["hamiltonian_residual_max"]),
        "projection_correction_max": float(ray["dispersion_projection_correction_max"]),
        "accepted_steps": int(ray["accepted_steps"]),
        "nfev": int(ray["nfev"]),
        "spatial_fd_relative_spread_max": float(
            ray["spatial_fd_relative_spread_max"]),
        "spatial_fd_absolute_spread_max": float(
            ray["spatial_fd_absolute_spread_max"]),
        "spatial_fd_stability_ratio_max": float(
            ray["spatial_fd_stability_ratio_max"]),
        "spatial_fd_selected_platform_counts": json.dumps(
            ray["spatial_fd_selected_platform_counts"], sort_keys=True),
        "failure_reason": ray["failure_reason"],
    }
    return record, ray


def _init_worker(input_path: str) -> None:
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = _load_context(input_path)


def _worker_trace(elevation_deg: float) -> dict:
    try:
        return _trace(_model(_WORKER_CONTEXT), _WORKER_CONTEXT, elevation_deg)[0]
    except Exception as error:
        return {
            "elevation_deg": float(elevation_deg), "status": "catastrophic_error",
            "constraint_ok": False, "ground_range_km": None,
            "along_track_residual_km": None, "cross_track_km": None,
            "apex_km": None, "turning_count": None, "geometric_path_km": None,
            "group_path_km": None, "tau_ms": None, "landing_lat_deg": None,
            "landing_lon_deg": None, "path_point_count": 0,
            "hamiltonian_residual_max": None, "projection_correction_max": None,
            "accepted_steps": 0, "nfev": 0, "spatial_fd_relative_spread_max": None,
            "spatial_fd_absolute_spread_max": None,
            "spatial_fd_stability_ratio_max": None,
            "spatial_fd_selected_platform_counts": "{}",
            "failure_reason": f"{type(error).__name__}: {error}",
        }


def _worker_trace_full(elevation_deg: float) -> tuple[dict, dict]:
    return _trace(_model(_WORKER_CONTEXT), _WORKER_CONTEXT, elevation_deg)


def _scan_angles(input_path: Path, angles: np.ndarray, workers: int) -> list[dict]:
    values = [float(value) for value in np.asarray(angles)]
    if workers == 1:
        context = _load_context(input_path)
        return [_trace(_model(context), context, angle)[0] for angle in values]
    with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker,
            initargs=(str(input_path),)) as executor:
        return list(executor.map(_worker_trace, values, chunksize=1))


def _adaptive_angles(coarse: list[dict]) -> np.ndarray:
    intervals = set()
    for left, right in zip(coarse, coarse[1:]):
        left_ground = left["status"] == "ground" and left["constraint_ok"]
        right_ground = right["status"] == "ground" and right["constraint_ok"]
        near_target = any(
            item["along_track_residual_km"] is not None
            and abs(item["along_track_residual_km"]) <= 250.0
            for item in (left, right)
        )
        sign_change = (
            left_ground and right_ground
            and left["along_track_residual_km"] * right["along_track_residual_km"] <= 0.0
        )
        if left["status"] != right["status"] or near_target or sign_change:
            intervals.add((left["elevation_deg"], right["elevation_deg"]))
    ground = [item for item in coarse
              if item["status"] == "ground" and item["constraint_ok"]]
    for a, b, c in zip(ground, ground[1:], ground[2:]):
        if (b["elevation_deg"] - a["elevation_deg"] <= 0.51
                and c["elevation_deg"] - b["elevation_deg"] <= 0.51):
            left_slope = b["ground_range_km"] - a["ground_range_km"]
            right_slope = c["ground_range_km"] - b["ground_range_km"]
            if left_slope * right_slope <= 0.0:
                intervals.update(((a["elevation_deg"], b["elevation_deg"]),
                                  (b["elevation_deg"], c["elevation_deg"])))
    angles = set()
    coarse_set = {round(item["elevation_deg"], 10) for item in coarse}
    for low, high in intervals:
        for angle in np.arange(low, high + 0.05, 0.1):
            key = round(float(angle), 10)
            if key not in coarse_set:
                angles.add(key)
    return np.asarray(sorted(angles), dtype=float)


def _root_brackets(records: list[dict]) -> list[tuple[float, float]]:
    ordered = sorted(records, key=lambda item: item["elevation_deg"])
    brackets = []
    for left, right in zip(ordered, ordered[1:]):
        valid = all(item["status"] == "ground" and item["constraint_ok"]
                    and item["along_track_residual_km"] is not None
                    for item in (left, right))
        if (valid and right["elevation_deg"] - left["elevation_deg"] <= 0.11
                and left["along_track_residual_km"]
                * right["along_track_residual_km"] <= 0.0):
            brackets.append((left["elevation_deg"], right["elevation_deg"]))
    return brackets


def _solve_roots(context: dict, brackets: list[tuple[float, float]]) -> tuple[list[dict], list[dict]]:
    cache: dict[float, tuple[dict, dict]] = {}

    def evaluate(angle: float) -> tuple[dict, dict]:
        key = round(float(angle), 12)
        if key not in cache:
            cache[key] = _trace(_model(context), context, angle)
        record, ray = cache[key]
        if (record["status"] != "ground" or not record["constraint_ok"]
                or record["along_track_residual_km"] is None):
            raise ValueError(f"root trial left valid ground branch at {angle:.12g} deg")
        return record, ray

    roots = []
    rays = []
    for low, high in brackets:
        try:
            angle = brentq(
                lambda value: evaluate(value)[0]["along_track_residual_km"],
                low, high, xtol=1.0e-10, rtol=1.0e-12, maxiter=60,
            )
            record, ray = evaluate(angle)
        except (ValueError, RuntimeError):
            continue
        if any(abs(angle - item["elevation_deg"]) <= 1.0e-6 for item in roots):
            continue
        label = f"elevation_order_candidate_{len(roots) + 1}"
        roots.append({
            **record, "candidate_label": label,
            "bracket_low_deg": float(low), "bracket_high_deg": float(high),
            "root_solver": "bounded_Brent_RX_ENU_along_track_ground_only",
        })
        rays.append(ray)
    return roots, rays


def _write_csv(path: Path, records: list[dict]) -> None:
    fields = [
        "scan_stage", "elevation_deg", "candidate_label", "status",
        "constraint_ok", "ground_range_km", "along_track_residual_km",
        "cross_track_km", "apex_km", "turning_count", "geometric_path_km",
        "group_path_km", "tau_ms", "landing_lat_deg", "landing_lon_deg",
        "path_point_count", "hamiltonian_residual_max", "projection_correction_max",
        "accepted_steps", "nfev", "spatial_fd_relative_spread_max",
        "spatial_fd_absolute_spread_max", "spatial_fd_stability_ratio_max",
        "spatial_fd_selected_platform_counts", "failure_reason", "bracket_low_deg",
        "bracket_high_deg", "root_solver",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def _trajectory_range_height(ray: dict, context: dict) -> tuple[np.ndarray, np.ndarray]:
    geocentric = np.asarray(ray["trajectory_geocentric_lat_lon_height"], dtype=float)
    ranges = np.array([
        _surface_distance_km(
            context["tx"][0], context["tx"][1], row[0], row[1], context["radius"])
        for row in geocentric
    ])
    return ranges, geocentric[:, 2]


def _plot_scan(coarse: list[dict], adaptive: list[dict], roots: list[dict],
               rays: list[dict], context: dict, source_hashes: dict[str, str]) -> dict:
    colors = {
        "ground": "#009E73", "upper_boundary": "#E69F00",
        "grid_exit": "#D55E00", "max_path": "#CC79A7",
        "physics_domain_rejected": "#000000", "constraint_failed": "#0072B2",
        "integration_failed": "#777777",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    for status in sorted({item["status"] for item in coarse}):
        selected = [item for item in coarse if item["status"] == status]
        color = colors.get(status, "#777777")
        ground_y = [item["ground_range_km"] for item in selected]
        axes[0, 0].scatter([item["elevation_deg"] for item in selected], ground_y,
                           s=15, color=color, label=status)
        axes[0, 1].scatter([item["elevation_deg"] for item in selected],
                           [item["apex_km"] for item in selected], s=15,
                           color=color, label=status)
    ground_adaptive = [item for item in adaptive if item["ground_range_km"] is not None]
    if ground_adaptive:
        axes[0, 0].scatter(
            [item["elevation_deg"] for item in ground_adaptive],
            [item["ground_range_km"] for item in ground_adaptive], s=8, marker="x",
            color="#0072B2", alpha=0.65, label="independent adaptive 0.1 deg")
    axes[0, 0].axhline(context["target_range_km"], color="black", linestyle="--",
                       linewidth=1.0, label="RX target range")
    for index, root in enumerate(roots):
        for axis, key in ((axes[0, 0], "ground_range_km"), (axes[0, 1], "apex_km")):
            axis.scatter(root["elevation_deg"], root[key], s=90, marker="*",
                         color="#D55E00", edgecolor="black", zorder=5,
                         label="AH3D bounded ground roots" if index == 0 else None)

    max_range = max(
        [context["target_range_km"] * 1.15]
        + [float(np.max(_trajectory_range_height(ray, context)[0])) for ray in rays]
    )
    height = context["height"]
    profile = context["ne"][context["ne"].shape[0] // 2,
                            context["ne"].shape[1] // 2, :] / 1.0e11
    background = np.broadcast_to(profile[:, None], (height.size, 240))
    range_axis = np.linspace(0.0, max_range, background.shape[1])
    mesh = axes[1, 0].pcolormesh(range_axis, height, background, shading="auto",
                                 cmap="jet", alpha=0.5)
    for index, ray in enumerate(rays):
        ranges, heights = _trajectory_range_height(ray, context)
        axes[1, 0].plot(ranges, heights, linewidth=1.8,
                        label=f"candidate {index + 1}")
    axes[1, 0].axvline(context["target_range_km"], color="black", linestyle="--",
                       linewidth=1.0, label="RX")
    fig.colorbar(mesh, ax=axes[1, 0], label="Ne (1e11 m^-3)")

    statuses = sorted({item["status"] for item in coarse})
    counts = [sum(item["status"] == status for item in coarse) for status in statuses]
    axes[1, 1].bar(np.arange(len(statuses)), counts,
                   color=[colors.get(status, "#777777") for status in statuses])
    axes[1, 1].set_xticks(np.arange(len(statuses)), statuses, rotation=20, ha="right")
    for x, count in enumerate(counts):
        axes[1, 1].text(x, count, str(count), ha="center", va="bottom", fontsize=8)

    axes[0, 0].set(xlabel="Launch elevation (deg)",
                   ylabel="Terminal ground range (km)", xlim=(1.0, 89.0))
    axes[0, 1].set(xlabel="Launch elevation (deg)", ylabel="Apex height (km)",
                   xlim=(1.0, 89.0))
    axes[1, 0].set(xlabel="Great-circle range from TX (km)",
                   ylabel="Physical height (km)", xlim=(0.0, max_range))
    axes[1, 1].set(xlabel="AH3D termination status", ylabel="Coarse ray count")
    for label, axis in zip(("A", "B", "C", "D"), axes.flat):
        axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top",
                  fontweight="bold")
    for axis in axes.flat[:3]:
        axis.legend(frameon=False, fontsize=7, loc="best")
    fig.suptitle(
        "A4-C1 Prep-3 | independent ECEF AH3D full-elevation scan\n"
        "one-way | Chapman F2 frozen grid | 10 MHz | O/X degenerate (B=0) | "
        "collision=0 | DOP853 research | not PHaRLAP comparison",
        fontsize=9.5,
    )
    fig.savefig(FIGURE, dpi=300, facecolor="white", metadata={
        "Title": "A4-C1 Prep-3 independent AH3D full-domain scan",
        "Description": "AH3D-only candidate-root validation; not A4-C1 comparison",
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


def _verify_existing(allow_rebuild: bool) -> int | None:
    if not MANIFEST.is_file():
        return None
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    ok = (manifest.get("status") == "passed_ah3d_full_domain_scan"
          and all(Path(item["path"]).is_file()
                  and _sha256(Path(item["path"])) == item["sha256"]
                  for item in manifest.get("frozen_artifacts", [])))
    if ok:
        print("A4-C1-Prep-3 frozen: OK")
        return 0
    if allow_rebuild:
        print("A4-C1-Prep-3 frozen artifacts changed; explicit rebuild enabled")
        return None
    print("A4-C1-Prep-3 frozen artifacts changed; rebuild refused")
    return 1


def _load_resume(input_hash: str, trace_provider_hash: str) -> tuple[
        list[dict], list[dict], list[dict], dict] | None:
    if not REPORT.is_file() or not CSV.is_file():
        return None
    previous = json.loads(REPORT.read_text(encoding="utf-8"))
    if (previous.get("status") not in {"passed_ah3d_full_domain_scan", "failed"}
            or previous.get("prep1_input", {}).get("sha256") != input_hash
            or previous.get("hashes", {}).get("scan_csv") != _sha256(CSV)):
        return None
    previous_hashes = previous.get("hashes", {})
    if previous.get("resume_contract_version") == 1 and (
            previous_hashes.get("python") != _sha256(Path(__file__))
            or previous_hashes.get("solver_source") != _sha256(SOLVER_SOURCE)
            or previous_hashes.get("provider_source") != _sha256(PROVIDER_SOURCE)
            or previous_hashes.get("trace_grid_with_boundary_continuation")
            != trace_provider_hash):
        return None
    string_fields = {
        "scan_stage", "candidate_label", "status", "failure_reason",
        "spatial_fd_selected_platform_counts", "root_solver",
    }
    integer_fields = {"turning_count", "path_point_count", "accepted_steps", "nfev"}
    records = []
    with CSV.open(newline="", encoding="utf-8") as stream:
        for raw in csv.DictReader(stream):
            record = {}
            for key, value in raw.items():
                if key in string_fields:
                    record[key] = value or None
                elif key == "constraint_ok":
                    record[key] = value.lower() == "true"
                elif value == "":
                    record[key] = None
                elif key in integer_fields:
                    record[key] = int(value)
                else:
                    record[key] = float(value)
            relative_spread = record.get("spatial_fd_relative_spread_max")
            record.setdefault("spatial_fd_absolute_spread_max", None)
            record.setdefault(
                "spatial_fd_stability_ratio_max",
                None if relative_spread is None else min(1.0, relative_spread / 0.02),
            )
            records.append(record)
    coarse = [item for item in records if item["scan_stage"] == "coarse_0.5_deg"]
    adaptive = [item for item in records if item["scan_stage"] == "adaptive_0.1_deg"]
    roots = previous.get("roots", [])
    if (len(coarse) != COARSE_ANGLES.size or len(roots) < 2
            or not np.allclose(
                [item["elevation_deg"] for item in coarse], COARSE_ANGLES,
                atol=1.0e-12)):
        return None
    provenance = {
        "previous_report_sha256": _sha256(REPORT),
        "previous_scan_csv_sha256": _sha256(CSV),
        "previous_trace_provider_sha256": previous_hashes.get(
            "trace_grid_with_boundary_continuation",
            previous_hashes.get("trace_grid_with_lower_continuation")),
        "current_trace_provider_sha256": trace_provider_hash,
        "prior_resume_provenance": previous.get("resume_provenance"),
        "pre_upper_continuation_ground_reuse": bool(
            (previous.get("resume_provenance") or {}).get(
                "reused_ground_records", 0)),
        "reused_ground_records": sum(
            item["status"] == "ground" and item["constraint_ok"]
            for item in records),
        "reused_ground_max_apex_km": max(
            item["apex_km"] for item in records
            if item["status"] == "ground" and item["apex_km"] is not None),
    }
    return coarse, adaptive, roots, provenance


def _retrace_roots(input_path: Path, previous_roots: list[dict],
                   workers: int) -> tuple[list[dict], list[dict]]:
    angles = [float(item["elevation_deg"]) for item in previous_roots]
    with ProcessPoolExecutor(
            max_workers=min(workers, len(angles)), initializer=_init_worker,
            initargs=(str(input_path),)) as executor:
        traced = list(executor.map(_worker_trace_full, angles, chunksize=1))
    roots, rays = [], []
    for previous, (record, ray) in zip(previous_roots, traced):
        roots.append({
            **record,
            "candidate_label": previous["candidate_label"],
            "bracket_low_deg": previous["bracket_low_deg"],
            "bracket_high_deg": previous["bracket_high_deg"],
            "root_solver": previous["root_solver"] + "+post_repair_fixed_angle_recheck",
        })
        rays.append(ray)
    return roots, rays


def run(*, allow_rebuild: bool = False, workers: int | None = None) -> int:
    existing = _verify_existing(allow_rebuild)
    if existing is not None:
        return existing
    OUT.mkdir(parents=True, exist_ok=True)
    prerequisites = {
        "prep1_report": PREP1_REPORT.is_file(),
        "prep1_manifest": PREP1_MANIFEST.is_file(),
        "production_profile_is_research": (
            cfg.UPDATE2.get("AH3D_INTEGRATION_PROFILE") == PROFILE),
        "production_solver_is_ecef_dop853": (
            cfg.UPDATE2.get("AH3D_SOLVER") == "ecef_dop853"),
    }
    if not all(prerequisites.values()):
        _atomic_json(REPORT, {"status": "blocked_prerequisite", "checks": prerequisites})
        return 2
    try:
        input_path, input_hash, prep1_input = _verify_prep1()
        context = _load_context(input_path)
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_prep1_hash", "error": str(error)})
        return 2

    worker_count = max(1, min(int(workers or min(os.cpu_count() or 1, 8)), 8))
    resume = (_load_resume(input_hash, context["trace_provider_hash"])
              if allow_rebuild else None)
    resume_provenance = None
    if resume is None:
        coarse = _scan_angles(input_path, COARSE_ANGLES, worker_count)
        adaptive_angles = _adaptive_angles(coarse)
        adaptive = (_scan_angles(input_path, adaptive_angles, worker_count)
                    if adaptive_angles.size else [])
        brackets = _root_brackets(coarse + adaptive)
        roots, rays = _solve_roots(context, brackets)
    else:
        coarse, adaptive, previous_roots, resume_provenance = resume
        records = coarse + adaptive
        retry_indices = [
            index for index, item in enumerate(records)
            if item["status"] not in {"ground", "upper_boundary"}
            or not item["constraint_ok"]
        ]
        retraced = _scan_angles(
            input_path,
            np.asarray([records[index]["elevation_deg"] for index in retry_indices]),
            worker_count,
        )
        for index, replacement in zip(retry_indices, retraced):
            replacement["scan_stage"] = records[index]["scan_stage"]
            records[index] = replacement
        coarse = records[:len(coarse)]
        adaptive = records[len(coarse):]
        brackets = _root_brackets(coarse + adaptive)
        roots, rays = _retrace_roots(input_path, previous_roots, worker_count)
        resume_provenance["recomputed_records"] = len(retry_indices)
    for record in coarse:
        record["scan_stage"] = "coarse_0.5_deg"
    for record in adaptive:
        record["scan_stage"] = "adaptive_0.1_deg"
    _write_csv(CSV, coarse + adaptive)
    _write_csv(ROOT_CSV, [{**root, "scan_stage": "bounded_root"} for root in roots])
    trajectory_payload = {}
    for index, ray in enumerate(rays, start=1):
        trajectory_payload[f"candidate_{index}_ecef_km"] = ray["trajectory_ecef_km"]
        trajectory_payload[f"candidate_{index}_geocentric"] = (
            ray["trajectory_geocentric_lat_lon_height"])
        trajectory_payload[f"candidate_{index}_geom_path_km"] = (
            ray["trajectory_geom_path_km"])
    np.savez_compressed(ROOT_NPZ, **trajectory_payload)
    source_hashes = {
        "scan_csv": _sha256(CSV), "root_csv": _sha256(ROOT_CSV),
        "root_trajectories": _sha256(ROOT_NPZ), "input": input_hash,
    }
    figure = _plot_scan(coarse, adaptive, roots, rays, context, source_hashes)

    expected = COARSE_ANGLES
    coarse_angles = np.array([item["elevation_deg"] for item in coarse])
    valid_roots = [
        item["status"] == "ground" and item["constraint_ok"]
        and abs(item["along_track_residual_km"]) <= ROOT_RESIDUAL_LIMIT_KM
        and abs(item["cross_track_km"]) <= ROOT_RESIDUAL_LIMIT_KM
        for item in roots
    ]
    status_counts = {
        status: sum(item["status"] == status for item in coarse + adaptive)
        for status in sorted({item["status"] for item in coarse + adaptive})
    }
    scan_records = coarse + adaptive
    allowed_scan_statuses = {"ground", "upper_boundary"}
    resume_equivalence = resume_provenance is None or (
        (not resume_provenance["pre_upper_continuation_ground_reuse"]
         and resume_provenance["previous_trace_provider_sha256"]
         == resume_provenance["current_trace_provider_sha256"])
        or (resume_provenance["pre_upper_continuation_ground_reuse"]
            and resume_provenance["reused_ground_max_apex_km"]
            <= H_STOP_KM - 10.0)
    )
    checks = {
        **prerequisites,
        "prep1_frozen_hash": True,
        "source_grid_hash_matches_prep1": (
            context["source_provider_hash"] == prep1_input["provider_hash"]),
        "frozen_source_nodes_unchanged": (
            np.array_equal(
                context["source_height"],
                context["height"][context["source_start"]:
                                  context["source_start"] + context["source_height"].size])
            and np.array_equal(
                context["source_ne"],
                context["ne"][:, :, context["source_start"]:
                              context["source_start"] + context["source_height"].size])),
        "chapman_continuation_only_outside_source_grid": (
            context["height"][0] == -float(context["source_height"][1]
                                            - context["source_height"][0])
            and context["height"][-1] == context["source_height"][-1] + 2.0 * float(
                context["source_height"][1] - context["source_height"][0])
            and np.all(np.isfinite(context["ne"]))
            and np.all(context["ne"] >= 0.0)),
        "complete_coarse_domain": (
            len(coarse) == expected.size and np.allclose(coarse_angles, expected, atol=1.0e-12)),
        "no_unrecorded_exception": all(item["status"] for item in coarse + adaptive),
        "only_ground_or_upper_boundary_status": all(
            item["status"] in allowed_scan_statuses for item in scan_records),
        "no_physics_domain_or_catastrophic_rejection": all(
            item["status"] not in {"physics_domain_rejected", "catastrophic_error"}
            for item in scan_records),
        "all_scan_samples_constraint_valid": all(
            item["constraint_ok"] for item in scan_records),
        "resume_lineage_physical_equivalence": resume_equivalence,
        "ground_and_upper_boundary_both_present": (
            {item["status"] for item in scan_records} == allowed_scan_statuses),
        "adaptive_scan_present": len(adaptive) > 0,
        "ground_only_brackets": all(
            low < high and high - low <= 0.11 for low, high in brackets),
        "at_least_two_roots": len(roots) >= 2,
        "all_roots_ground_constraint_and_within_1km": bool(valid_roots) and all(valid_roots),
        "all_ground_samples_constraint_valid": all(
            item["constraint_ok"] for item in coarse + adaptive
            if item["status"] == "ground"),
        "root_trajectories_saved": ROOT_NPZ.is_file() and ROOT_NPZ.stat().st_size > 0,
        "visualization": all(figure["checks"].values()),
        "legacy_and_pharlap_roots_not_used": True,
    }
    status = "passed_ah3d_full_domain_scan" if all(checks.values()) else "failed"
    report = {
        "schema_version": "a4-c1-prep3-validation-v1",
        "resume_contract_version": 1, "status": status,
        "case_semantics": "ah3d_independent_full_domain_scan_not_pharlap_comparison_or_c1_validation",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "configuration": {
            "elevation_domain_deg": [1.0, 89.0], "coarse_step_deg": 0.5,
            "adaptive_step_deg": 0.1, "bearing_deg": context["bearing_deg"],
            "profile": PROFILE, "integrator": "DOP853", "h_stop_km": H_STOP_KM,
            "max_path_km": MAX_PATH_KM, "wave_mode": "O/X_degenerate_single_count",
            "spatial_fd_base_km": 0.5,
            "spatial_fd_scales": [
                0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0],
            "spatial_fd_stability_rtol": 0.02,
            "spatial_fd_stability_atol_n2_per_km": 2.0e-5,
            "B_T": 0.0, "collision_Hz": 0.0, "workers": worker_count,
            "boundary_continuation_semantics": (
                "same-Chapman 2-km continuation below ground and through 504 km; "
                "frozen 60-500 km Ne nodes unchanged"),
        },
        "prep1_input": {"path": str(input_path), "sha256": input_hash},
        "hashes": {
            "python": _sha256(Path(__file__)),
            "solver_source": _sha256(SOLVER_SOURCE),
            "provider_source": _sha256(PROVIDER_SOURCE),
            "source_grid": context["source_provider_hash"],
            "trace_grid_with_boundary_continuation": context["trace_provider_hash"],
            **source_hashes,
        },
        "status_counts": status_counts,
        "resume_provenance": resume_provenance,
        "coarse_ray_count": len(coarse), "adaptive_ray_count": len(adaptive),
        "root_brackets_deg": brackets, "roots": roots,
        "metrics": {
            "maximum_hamiltonian_residual": max(
                item["hamiltonian_residual_max"] for item in coarse + adaptive
                if item["hamiltonian_residual_max"] is not None),
            "maximum_spatial_fd_relative_spread": max(
                item["spatial_fd_relative_spread_max"] for item in coarse + adaptive
                if item["spatial_fd_relative_spread_max"] is not None),
            "maximum_spatial_fd_stability_ratio": max(
                item["spatial_fd_stability_ratio_max"] for item in coarse + adaptive
                if item["spatial_fd_stability_ratio_max"] is not None),
            "maximum_ground_cross_track_km": max(
                [abs(item["cross_track_km"]) for item in coarse + adaptive
                 if item["cross_track_km"] is not None] or [0.0]),
        },
        "artifacts": {"scan_csv": str(CSV), "root_csv": str(ROOT_CSV),
                      "root_trajectories_npz": str(ROOT_NPZ), "figure": figure},
        "limitations": [
            "Prep-3 is AH3D-only and does not compare trajectories with PHaRLAP",
            "Candidates are ordered by elevation; low/high variational semantics remain Prep-4 work",
            "The B=0 horizontally uniform case does not validate O/X splitting or complex 3-D fields",
            "Same-formula boundary continuation is numerical stencil support only; the frozen PHaRLAP 60-500 km Ne grid is unchanged",
        ],
    }
    _atomic_json(REPORT, report)
    if status != "passed_ah3d_full_domain_scan":
        print(f"A4-C1-Prep-3: FAIL roots={len(roots)} report={REPORT}")
        return 1
    artifacts = [CSV, ROOT_CSV, ROOT_NPZ, FIGURE, Path(__file__), SOLVER_SOURCE,
                 PROVIDER_SOURCE, input_path]
    _atomic_json(MANIFEST, {
        "schema_version": "a4-c1-prep3-freeze-v1", "status": status,
        "case_semantics": report["case_semantics"],
        "validation_report": {"path": str(REPORT), "sha256": _sha256(REPORT)},
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in artifacts],
    })
    print(f"A4-C1-Prep-3: OK roots={len(roots)} report={REPORT}")
    return 0


if __name__ == "__main__":
    arguments = sys.argv[1:]
    allow_rebuild = "--rebuild-frozen" in arguments
    workers = None
    for argument in arguments:
        if argument.startswith("--workers="):
            workers = int(argument.split("=", 1)[1])
        elif argument != "--rebuild-frozen":
            raise SystemExit(
                "usage: _run_a4_c1_prep3.py [--rebuild-frozen] [--workers=N]")
    raise SystemExit(run(allow_rebuild=allow_rebuild, workers=workers))
