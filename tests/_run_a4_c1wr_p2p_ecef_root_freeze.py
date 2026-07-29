"""Freeze current strict spherical-ECEF cubic-Chapman P2P roots."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.ne_provider_geocentric import GeocentricGridNeProvider3D  # noqa: E402
from models.point_to_point_ah_3d_ecef import solve_p2p_ray_ah3d_ecef  # noqa: E402
from models.ray_tracer_ah_3d_ecef import (  # noqa: E402
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)
from tests._run_a4_c1_prep3 import _load_context  # noqa: E402
from utils import enu_basis, geocentric_to_ecef, rx_enu_residual  # noqa: E402


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr_p2pf"
INPUT = (ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep1" /
         "inputs" / "a4_c1_chapman_refined_1km.mat")
C1WR_MANIFEST = (ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr" /
                  "a4_c1wr_manifest.json")
PREP5_REPORT = (ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep5" /
                "a4_c1_prep5_validation.json")
P0C_REFERENCE = (ROOT / "output2" / "a3r_benchmark" /
                  "a3r_p0c_reference_preparation.json")
REPORT = OUT / "a4_c1wr_p2pf_result.json"
CSV = OUT / "a4_c1wr_p2pf_roots.csv"
NPZ = OUT / "a4_c1wr_p2pf_trajectories.npz"
LOG = OUT / "a4_c1wr_p2pf_run.log"
MANIFEST = OUT / "a4_c1wr_p2pf_manifest.json"
FIG_TRAJECTORIES = OUT / "a4_c1wr_p2pf_root_trajectories.png"
FIG_CLOSURE = OUT / "a4_c1wr_p2pf_solver_closure.png"
FIG_MIGRATION = OUT / "a4_c1wr_p2pf_root_migration.png"
FIG_GATES = OUT / "a4_c1wr_p2pf_acceptance_gates.png"

SEMANTICS = "current_strict_spherical_ECEF_cubic_Chapman_P2P_reference"
BRANCHES = {
    "low_ray": (21.9, 22.0),
    "high_ray": (39.5, 39.6),
}
PREP5_SEEDS = {
    "low_ray": 21.94562872092976,
    "high_ray": 39.58701138438118,
}
LEGACY_ROOTS = {
    "low_ray": 35.33980751304579,
    "high_ray": 39.918776604143765,
}
PROFILES = ("research", "reference")
H_STOP_KM = 500.0
MAX_PATH_KM = 10000.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic)
        else item.tolist(),
    ) + "\n", encoding="utf-8")


def _verify_frozen_manifest(path: Path, expected_status: str) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != expected_status:
        raise RuntimeError(f"unexpected prerequisite status: {path}")
    artifacts = manifest.get("artifact_sha256", {})
    if not artifacts:
        artifacts = {
            item["path"]: item["sha256"]
            for item in manifest.get("frozen_artifacts", [])
        }
    for name, expected in artifacts.items():
        artifact = Path(name)
        if not artifact.is_absolute():
            artifact = ROOT / artifact
        if not artifact.is_file() or _sha256(artifact) != expected:
            raise RuntimeError(f"frozen prerequisite changed: {artifact}")
    return manifest


def _load_cubic_context() -> dict:
    context = _load_context(INPUT)
    source = context["provider"]
    reference_profile = source.Ne_m3[0, 0, :]
    if not np.array_equal(
            source.Ne_m3,
            np.broadcast_to(reference_profile, source.Ne_m3.shape)):
        raise RuntimeError(
            "P2PF exact horizontal factorization requires a horizontally "
            "uniform frozen Chapman field")
    latitude_indices = np.unique(np.linspace(
        0, source.latitude_deg.size - 1, 4, dtype=int))
    longitude_indices = np.unique(np.linspace(
        0, source.longitude_deg.size - 1, 4, dtype=int))
    if latitude_indices.size != 4 or longitude_indices.size != 4:
        raise RuntimeError("cubic factorization requires four horizontal nodes")
    factorized_ne = source.Ne_m3[latitude_indices][:, longitude_indices, :]
    context["provider"] = GeocentricGridNeProvider3D(
        factorized_ne,
        source.latitude_deg[latitude_indices],
        source.longitude_deg[longitude_indices],
        source.height_km,
        source.earth_radius_km,
        fill_value_m3=source.fill_value_m3,
        longitude_periodic_complete=source.longitude_periodic_complete,
        source_metadata={
            **source.source_metadata,
            "p2pf_semantics": SEMANTICS,
            "frozen_refined_1km_input_sha256": _sha256(INPUT),
            "horizontal_factorization": (
                "exact_constant_field_149x17_to_4x4;"
                "all_source_horizontal_nodes_verified_bitwise_equal"),
        },
        interpolation="cubic",
    )
    context["source_trace_grid_shape"] = list(source.Ne_m3.shape)
    context["factorized_trace_grid_shape"] = list(factorized_ne.shape)
    context["linear_trace_provider_hash"] = context["trace_provider_hash"]
    context["trace_provider_hash"] = context["provider"].grid_hash
    return context


def _model(context: dict, profile: str) -> RefractiveIndexAH3DECEF:
    return RefractiveIndexAH3DECEF(
        context["provider"],
        earth_radius_km=context["radius"],
        freq_MHz=context["frequency_MHz"],
        wave_mode="O",
        gyro_MHz=0.0,
        b_ecef_hat=(0.0, 0.0, 1.0),
        collision_s=0.0,
        spatial_fd_base_km=0.25,
        spatial_fd_scales=(0.5, 1.0, 2.0),
        integration_profile=profile,
    )


def _trace(context: dict, elevation_deg: float, profile: str) -> dict:
    return trace_ray_ah3d_ecef(
        _model(context, profile),
        tx_lat_deg=float(context["tx"][0]),
        tx_lon_deg=float(context["tx"][1]),
        tx_height_km=float(context["tx"][2]),
        elevation_deg=float(elevation_deg),
        bearing_deg=float(context["bearing_deg"]),
        h_stop_km=H_STOP_KM,
        max_path_km=MAX_PATH_KM,
        profile=profile,
    )


def _ray_record(context: dict, ray: dict, *, branch: str, profile: str,
                solver: str, elevation_deg: float, elapsed_s: float) -> dict:
    terminal = np.asarray(ray["trajectory_ecef_km"], dtype=float)[-1]
    receiver = geocentric_to_ecef(*context["rx"], context["radius"])
    enu = rx_enu_residual(
        terminal, *context["rx"], context["radius"])
    return {
        "branch": branch,
        "profile": profile,
        "solver": solver,
        "elevation_deg": float(elevation_deg),
        "bearing_deg": float(ray["bearing_tx_deg"]),
        "termination_status": ray["termination_status"],
        "turning_point_count": int(ray["turning_point_count"]),
        "rx_ecef_error_km": float(np.linalg.norm(terminal - receiver)),
        "rx_enu_east_km": float(enu[0]),
        "rx_enu_north_km": float(enu[1]),
        "rx_enu_error_km": float(np.linalg.norm(enu)),
        "apex_km": float(ray["apex_height_km"]),
        "geom_path_km": float(ray["geom_path_km"]),
        "group_path_km": float(ray["group_path_km"]),
        "tau_ms": float(ray["tau_ms"]),
        "hamiltonian_residual_max": float(ray["hamiltonian_residual_max"]),
        "spatial_fd_stability_ratio_max": float(
            ray["spatial_fd_stability_ratio_max"]),
        "spatial_fd_rejection_count": int(ray["spatial_fd_rejection_count"]),
        "accepted_steps": int(ray["accepted_steps"]),
        "nfev": int(ray["nfev"]),
        "elapsed_s": float(elapsed_s),
        "failure_reason": ray["failure_reason"],
    }


def _valid_root_ray(ray: dict) -> bool:
    return bool(
        ray["termination_status"] == "ground"
        and ray["turning_point_count"] == 1
        and ray["failure_reason"] is None
        and ray["spatial_fd_rejection_count"] == 0
        and np.isfinite(ray["hamiltonian_residual_max"])
        and ray["hamiltonian_residual_max"] <= 1.0e-8
        and np.isfinite(ray["spatial_fd_stability_ratio_max"])
        and ray["spatial_fd_stability_ratio_max"] <= 1.0
    )


def _solve_brent(context: dict, branch: str, bracket: tuple[float, float],
                 profile: str) -> tuple[dict, dict, dict]:
    cache: dict[float, tuple[dict, float]] = {}

    def evaluate(angle: float) -> tuple[dict, float]:
        key = round(float(angle), 12)
        if key not in cache:
            started = time.perf_counter()
            ray = _trace(context, key, profile)
            cache[key] = ray, time.perf_counter() - started
        ray, elapsed = cache[key]
        if not _valid_root_ray(ray):
            raise ValueError(
                f"{profile}/{branch} left valid ground branch at {angle:.12g}")
        terminal = np.asarray(ray["trajectory_ecef_km"], dtype=float)[-1]
        east, north = rx_enu_residual(
            terminal, *context["rx"], context["radius"])
        bearing = np.radians(context["bearing_deg"])
        return ray, float(east * np.sin(bearing) + north * np.cos(bearing))

    low, high = bracket
    low_ray, low_value = evaluate(low)
    high_ray, high_value = evaluate(high)
    if low_value * high_value > 0.0:
        raise RuntimeError(
            f"input_or_code_drift_blocked:{profile}/{branch}:"
            f"endpoint signs {low_value:.9g},{high_value:.9g}")
    started = time.perf_counter()
    angle = brentq(
        lambda value: evaluate(value)[1],
        low, high, xtol=1.0e-6, rtol=1.0e-12, maxiter=60)
    solve_elapsed = time.perf_counter() - started
    ray, _ = evaluate(angle)
    total_trace_elapsed = sum(value[1] for value in cache.values())
    record = _ray_record(
        context, ray, branch=branch, profile=profile, solver="bounded_Brent",
        elevation_deg=angle, elapsed_s=total_trace_elapsed)
    diagnostics = {
        "bracket_deg": [float(low), float(high)],
        "endpoint_residual_km": [low_value, high_value],
        "root_xtol_deg": 1.0e-6,
        "function_evaluations": len(cache),
        "forward_integration_elapsed_s": total_trace_elapsed,
        "brent_overhead_elapsed_s": max(0.0, solve_elapsed),
        "endpoint_status": [
            low_ray["termination_status"], high_ray["termination_status"]],
    }
    return record, ray, diagnostics


def _solve_trf(context: dict, branch: str, seed: float) -> tuple[dict | None,
                                                                  dict | None,
                                                                  float]:
    started = time.perf_counter()
    ray = solve_p2p_ray_ah3d_ecef(
        _model(context, "reference"),
        tx_geocentric=context["tx"],
        rx_geocentric=context["rx"],
        elevation_seed_deg=float(seed),
        bearing_seed_deg=float(context["bearing_deg"]),
        h_stop_km=H_STOP_KM,
        max_path_km=MAX_PATH_KM,
        elevation_bounds_deg=BRANCHES[branch],
        bearing_half_width_deg=0.1,
        jacobian_scales=(0.5, 1.0, 2.0),
        tol_km=1.0,
        max_nfev=6,
        profile="reference",
        seed_semantics=SEMANTICS,
    )
    elapsed = time.perf_counter() - started
    if ray is None:
        return None, None, elapsed
    record = _ray_record(
        context, ray, branch=branch, profile="reference", solver="bounded_TRF_2D",
        elevation_deg=ray["elevation_tx_deg"], elapsed_s=elapsed)
    record.update({
        "p2p_nfev": int(ray["p2p_nfev"]),
        "p2p_njev": int(ray["p2p_njev"]),
        "jacobian_diagnostics": ray["p2p_jacobian_reference"],
    })
    return record, ray, elapsed


def _range_height(context: dict, ray: dict) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(ray["trajectory_ecef_km"], dtype=float)
    unit = points / np.linalg.norm(points, axis=1)[:, None]
    tx = geocentric_to_ecef(*context["tx"], context["radius"])
    tx_unit = tx / np.linalg.norm(tx)
    ranges = context["radius"] * np.arccos(
        np.clip(unit @ tx_unit, -1.0, 1.0))
    return ranges, np.linalg.norm(points, axis=1) - context["radius"]


def _plot_trajectories(context: dict, reference_rays: dict,
                       legacy_rays: dict) -> None:
    figure, axis = plt.subplots(figsize=(10.5, 6.2), constrained_layout=True)
    heights = np.asarray(context["source_height"], dtype=float)
    source_ne = np.asarray(context["source_ne"], dtype=float)
    profile = source_ne[source_ne.shape[0] // 2,
                        source_ne.shape[1] // 2, :] / 1.0e11
    ranges = np.linspace(0.0, context["target_range_km"] * 1.1, 360)
    mesh = axis.pcolormesh(
        ranges, heights,
        np.broadcast_to(profile[:, None], (heights.size, ranges.size)),
        shading="auto", cmap="jet", alpha=0.5)
    colors = {"low_ray": "#0072B2", "high_ray": "#D55E00"}
    for branch in BRANCHES:
        for ray, style, label in (
                (reference_rays[branch], "-", f"current cubic {branch}"),
                (legacy_rays[branch], "--", f"legacy root forward {branch}")):
            distance, height = _range_height(context, ray)
            axis.plot(distance, height, color=colors[branch], linestyle=style,
                      linewidth=1.8 if style == "-" else 1.2, label=label)
    axis.axvline(context["target_range_km"], color="black", linestyle=":",
                 label="RX")
    axis.set(
        xlim=(0.0, context["target_range_km"] * 1.1), ylim=(0.0, 460.0),
        xlabel="Great-circle range from TX (km)",
        ylabel="Physical height (km)",
        title="P2PF current cubic ECEF roots vs legacy-coordinate roots")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.colorbar(mesh, ax=axis, label="Ne (1e11 m^-3)")
    figure.savefig(FIG_TRAJECTORIES, dpi=240, facecolor="white")
    plt.close(figure)


def _plot_closure(rows: list[dict]) -> None:
    by_key = {(row["branch"], row["profile"], row["solver"]): row for row in rows}
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    labels = list(BRANCHES)
    x = np.arange(len(labels))
    width = 0.25
    series = [
        ("research", "bounded_Brent", "#56B4E9"),
        ("reference", "bounded_Brent", "#0072B2"),
        ("reference", "bounded_TRF_2D", "#D55E00"),
    ]
    for index, (profile, solver, color) in enumerate(series):
        selected = [by_key.get((branch, profile, solver)) for branch in labels]
        angles = [np.nan if row is None else row["elevation_deg"] for row in selected]
        residuals = [np.nan if row is None else row["rx_enu_error_km"]
                     for row in selected]
        offset = (index - 1) * width
        axes[0].bar(x + offset, angles, width, color=color,
                    label=f"{profile} {solver}")
        axes[1].bar(x + offset, residuals, width, color=color,
                    label=f"{profile} {solver}")
    axes[0].set_ylabel("Elevation (deg)")
    axes[1].set_ylabel("RX ENU residual (km)")
    axes[1].axhline(1.0, color="black", linestyle="--", label="1 km gate")
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_title("Research/reference and Brent/TRF angles")
    axes[1].set_title("Receiver closure")
    figure.savefig(FIG_CLOSURE, dpi=240, facecolor="white")
    plt.close(figure)


def _plot_migration(reference_rows: dict, prep5: dict, legacy_rows: dict) -> None:
    metrics = [
        ("elevation_deg", "Elevation difference (deg)"),
        ("apex_km", "Apex difference (km)"),
        ("geom_path_km", "Geometric-path difference (km)"),
        ("tau_ms", "Group-delay difference (ms)"),
    ]
    prep5_rows = {
        item["branch_id"]: item for item in prep5["frozen_signatures"]
        if item["solver"] == "AH3D"
    }
    prep5_keys = {
        "elevation_deg": "elevation_deg",
        "apex_km": "apex_km",
        "geom_path_km": "geometric_path_km",
        "tau_ms": "tau_ms",
    }
    figure, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), constrained_layout=True)
    x = np.arange(len(BRANCHES))
    for axis, (key, title) in zip(axes.flat, metrics):
        versus_linear = [
            reference_rows[branch][key] - prep5_rows[branch][prep5_keys[key]]
            for branch in BRANCHES]
        versus_legacy = [
            reference_rows[branch][key] - legacy_rows[branch][key]
            for branch in BRANCHES]
        axis.bar(x - 0.18, versus_linear, 0.36, color="#0072B2",
                 label="current cubic - Prep5 linear")
        axis.bar(x + 0.18, versus_legacy, 0.36, color="#D55E00",
                 label="current cubic - legacy root forward")
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set_xticks(x, BRANCHES)
        axis.set_ylabel(title)
        axis.grid(axis="y", alpha=0.25)
        axis.legend(fontsize=7)
    figure.suptitle("P2PF coordinate/field representation migration")
    figure.savefig(FIG_MIGRATION, dpi=240, facecolor="white")
    plt.close(figure)


def _plot_gates(gates: dict) -> None:
    names = list(gates)
    values = np.array([[1.0 if gates[name] else 0.0] for name in names])
    figure, axis = plt.subplots(
        figsize=(10.8, max(4.8, 0.45 * len(names))),
        constrained_layout=True)
    axis.imshow(values, aspect="auto", vmin=0.0, vmax=1.0,
                cmap=matplotlib.colors.ListedColormap(["#D55E00", "#009E73"]))
    axis.set_xticks([0], ["Acceptance"])
    axis.set_yticks(np.arange(len(names)), names)
    for index, name in enumerate(names):
        axis.text(0, index, "PASS" if gates[name] else "FAIL",
                  color="white", ha="center", va="center", fontweight="bold")
    axis.set_title("A4-C1WR-P2PF acceptance gate matrix")
    figure.savefig(FIG_GATES, dpi=240, facecolor="white")
    plt.close(figure)


def _image_record(path: Path) -> dict:
    image = mpimg.imread(path)
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": _sha256(path),
        "shape": list(image.shape),
        "pixel_std": float(np.std(image[..., :3])),
    }


def _source_hashes() -> dict:
    paths = [
        Path(__file__),
        ROOT / "tests" / "_verify_a4_c1wr_p2p_ecef_root_freeze.py",
        ROOT / "models" / "ne_provider_geocentric.py",
        ROOT / "models" / "point_to_point_ah_3d_ecef.py",
        ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        ROOT / "tests" / "_run_a4_c1_prep3.py",
        INPUT,
        C1WR_MANIFEST,
        PREP5_REPORT,
        P0C_REFERENCE,
    ]
    return {str(path.relative_to(ROOT)): _sha256(path) for path in paths}


def _write_csv(rows: list[dict]) -> None:
    fields = [
        "branch", "profile", "solver", "elevation_deg", "bearing_deg",
        "termination_status", "turning_point_count", "rx_ecef_error_km",
        "rx_enu_east_km", "rx_enu_north_km", "rx_enu_error_km", "apex_km",
        "geom_path_km", "group_path_km", "tau_ms",
        "hamiltonian_residual_max", "spatial_fd_stability_ratio_max",
        "spatial_fd_rejection_count", "accepted_steps", "nfev", "p2p_nfev",
        "p2p_njev", "elapsed_s", "failure_reason",
    ]
    with CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _read_csv() -> list[dict]:
    text_fields = {"branch", "profile", "solver", "termination_status",
                   "failure_reason"}
    integer_fields = {
        "turning_point_count", "spatial_fd_rejection_count", "accepted_steps",
        "nfev", "p2p_nfev", "p2p_njev",
    }
    with CSV.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in list(row.items()):
            if key in text_fields:
                row[key] = value or None
            elif not value:
                row[key] = None
            elif key in integer_fields:
                row[key] = int(value)
            else:
                row[key] = float(value)
    return rows


def _finalize_existing() -> int:
    if not CSV.is_file() or not NPZ.is_file():
        raise RuntimeError("saved P2PF roots and trajectories are required")
    _verify_frozen_manifest(
        C1WR_MANIFEST, "identified_and_repaired_grid_derivative")
    context = _load_cubic_context()
    prep5 = json.loads(PREP5_REPORT.read_text(encoding="utf-8"))
    rows = _read_csv()
    by_key = {(row["branch"], row["profile"], row["solver"]): row for row in rows}
    expected = {
        (branch, profile, "bounded_Brent")
        for branch in BRANCHES for profile in PROFILES
    } | {(branch, "reference", "bounded_TRF_2D") for branch in BRANCHES}
    if set(by_key) != expected:
        raise RuntimeError("saved P2PF root table is incomplete")

    legacy_rows = {}
    if REPORT.is_file():
        prior = json.loads(REPORT.read_text(encoding="utf-8"))
        if (prior.get("semantics") == SEMANTICS
                and prior.get("configuration", {}).get("provider_hash")
                == context["trace_provider_hash"]):
            legacy_rows = prior.get("legacy_forward_rows", {})
    if set(legacy_rows) != set(BRANCHES):
        legacy_rows = {}
        for branch, angle in LEGACY_ROOTS.items():
            started = time.perf_counter()
            ray = _trace(context, angle, "reference")
            legacy_rows[branch] = _ray_record(
                context, ray, branch=branch, profile="reference",
                solver="legacy_coordinate_reference_only", elevation_deg=angle,
                elapsed_s=time.perf_counter() - started)

    reference_rows = {
        branch: by_key[(branch, "reference", "bounded_Brent")]
        for branch in BRANCHES}
    research_rows = {
        branch: by_key[(branch, "research", "bounded_Brent")]
        for branch in BRANCHES}
    trf_rows = {
        branch: by_key[(branch, "reference", "bounded_TRF_2D")]
        for branch in BRANCHES}
    profile_differences = {
        branch: {
            "elevation_abs_deg": abs(
                reference_rows[branch]["elevation_deg"]
                - research_rows[branch]["elevation_deg"]),
            "apex_abs_km": abs(
                reference_rows[branch]["apex_km"]
                - research_rows[branch]["apex_km"]),
            "tau_relative_fraction": abs(
                reference_rows[branch]["tau_ms"]
                - research_rows[branch]["tau_ms"])
                / reference_rows[branch]["tau_ms"],
        }
        for branch in BRANCHES
    }
    solver_differences = {
        branch: {
            "elevation_abs_deg": abs(
                reference_rows[branch]["elevation_deg"]
                - trf_rows[branch]["elevation_deg"]),
            "bearing_abs_deg": abs(
                trf_rows[branch]["bearing_deg"] - context["bearing_deg"]),
            "rx_enu_error_km": trf_rows[branch]["rx_enu_error_km"],
        }
        for branch in BRANCHES
    }
    gates = {
        "two_unique_reference_brent_roots": (
            len({row["elevation_deg"] for row in reference_rows.values()}) == 2),
        "ground_single_turning": all(
            row["termination_status"] == "ground"
            and row["turning_point_count"] == 1 for row in rows),
        "high_apex_and_delay_exceed_low": (
            reference_rows["high_ray"]["apex_km"]
            > reference_rows["low_ray"]["apex_km"]
            and reference_rows["high_ray"]["tau_ms"]
            > reference_rows["low_ray"]["tau_ms"]),
        "rx_ecef_enu_error_le_1_km": all(
            row["rx_ecef_error_km"] <= 1.0
            and row["rx_enu_error_km"] <= 1.0 for row in rows),
        "hamiltonian_residual_le_1e_8": all(
            row["hamiltonian_residual_max"] <= 1.0e-8 for row in rows),
        "spatial_fd_platform_accepted": all(
            row["spatial_fd_stability_ratio_max"] <= 1.0
            and row["spatial_fd_rejection_count"] == 0 for row in rows),
        "brent_trf_angle_le_0_02_deg": all(
            value["elevation_abs_deg"] <= 0.02
            for value in solver_differences.values()),
        "no_branch_exchange": all(
            BRANCHES[branch][0] <= trf_rows[branch]["elevation_deg"]
            <= BRANCHES[branch][1] for branch in BRANCHES),
        "research_reference_elevation_le_0_02_deg": all(
            value["elevation_abs_deg"] <= 0.02
            for value in profile_differences.values()),
        "research_reference_apex_le_0_5_km": all(
            value["apex_abs_km"] <= 0.5
            for value in profile_differences.values()),
        "research_reference_delay_le_0_2_percent": all(
            value["tau_relative_fraction"] <= 0.002
            for value in profile_differences.values()),
        "roots_inside_frozen_brackets": all(
            BRANCHES[row["branch"]][0] <= row["elevation_deg"]
            <= BRANCHES[row["branch"]][1] for row in rows),
    }
    status = (
        "passed_current_ecef_chapman_p2p_root_freeze"
        if all(gates.values())
        else "failed_current_ecef_chapman_p2p_acceptance")
    prep5_rows = {
        item["branch_id"]: item for item in prep5["frozen_signatures"]
        if item["solver"] == "AH3D"}
    migration = {
        branch: {
            "from_prep5_linear": {
                key: reference_rows[branch][target]
                - prep5_rows[branch][source]
                for key, target, source in (
                    ("elevation_deg", "elevation_deg", "elevation_deg"),
                    ("apex_km", "apex_km", "apex_km"),
                    ("geom_path_km", "geom_path_km", "geometric_path_km"),
                    ("group_path_km", "group_path_km", "group_path_km"),
                    ("tau_ms", "tau_ms", "tau_ms"))
            },
            "from_legacy_p0c": {
                "elevation_deg": (
                    reference_rows[branch]["elevation_deg"]
                    - LEGACY_ROOTS[branch]),
                "legacy_forward_rx_error_km": legacy_rows[branch][
                    "rx_enu_error_km"],
            },
        }
        for branch in BRANCHES
    }
    _plot_migration(reference_rows, prep5, legacy_rows)
    _plot_gates(gates)
    source_hashes = _source_hashes()
    brent_diagnostics = {
        f"{row['profile']}_{row['branch']}": {
            "bracket_deg": list(BRANCHES[row["branch"]]),
            "root_xtol_deg": 1.0e-6,
            "forward_integration_elapsed_s": row["elapsed_s"],
            "recovery_note": (
                "root and trajectory checkpoint survived a plotting-only "
                "failure; endpoint values and evaluation count were not persisted"),
        }
        for row in rows if row["solver"] == "bounded_Brent"
    }
    report = {
        "schema_version": "a4-c1wr-p2pf-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "semantics": SEMANTICS,
        "recovery": (
            "finalized_from_saved_CSV_NPZ_after_plot_field_mapping_failure;"
            "Brent_and_TRF_not_recomputed;legacy_forward_only_recomputed"),
        "limitations": [
            "not WGS84 or IGRF production validation",
            "not full two-dimensional angular-domain mode recall",
            "not O/X magnetoionic splitting or complex three-dimensional P2P",
            "legacy roots are forward-traced references only",
        ],
        "configuration": {
            "earth_radius_km": context["radius"],
            "frequency_MHz": context["frequency_MHz"],
            "B": 0,
            "collision": 0,
            "one_way": True,
            "bearing_deg": context["bearing_deg"],
            "brackets_deg": BRANCHES,
            "brent_xtol_deg": 1.0e-6,
            "trf_bearing_half_width_deg": 0.1,
            "trf_max_nfev": 6,
            "spatial_fd_steps_km": [0.125, 0.25, 0.5],
            "provider_interpolation": "cubic",
            "provider_hash": context["trace_provider_hash"],
            "prior_linear_provider_hash": context["linear_trace_provider_hash"],
            "source_trace_grid_shape": context["source_trace_grid_shape"],
            "factorized_trace_grid_shape": context[
                "factorized_trace_grid_shape"],
            "horizontal_factorization": (
                "exact_after_bitwise_horizontal_uniformity_check"),
        },
        "prep5_seeds_deg": PREP5_SEEDS,
        "legacy_coordinate_reference_only_deg": LEGACY_ROOTS,
        "rows": rows,
        "legacy_forward_rows": legacy_rows,
        "brent_diagnostics": brent_diagnostics,
        "profile_differences": profile_differences,
        "solver_differences": solver_differences,
        "migration": migration,
        "gates": gates,
        "timings": {
            "brent_forward_integration_s": sum(
                row["elapsed_s"] for row in rows
                if row["solver"] == "bounded_Brent"),
            "trf_total_s": sum(
                row["elapsed_s"] for row in rows
                if row["solver"] == "bounded_TRF_2D"),
            "legacy_forward_recovery_s": sum(
                row["elapsed_s"] for row in legacy_rows.values()),
            "saved_numeric_run_s_excluding_legacy_and_plots": sum(
                row["elapsed_s"] for row in rows),
        },
        "source_sha256": source_hashes,
        "figures": [
            _image_record(path) for path in (
                FIG_TRAJECTORIES, FIG_CLOSURE, FIG_MIGRATION, FIG_GATES)],
    }
    _write_json(REPORT, report)
    LOG.write_text(
        "\n".join([
            "P2PF numerical roots recovered from saved CSV/NPZ.",
            "Brent and TRF were not recomputed after plotting-only failure.",
            *(f"{row['branch']} {row['profile']} {row['solver']} "
              f"elevation={row['elevation_deg']:.12g} "
              f"rx_error_km={row['rx_enu_error_km']:.12g}"
              for row in rows),
            f"status={status}",
        ]) + "\n", encoding="utf-8")
    artifacts = [REPORT, CSV, NPZ, LOG, FIG_TRAJECTORIES, FIG_CLOSURE,
                 FIG_MIGRATION, FIG_GATES]
    _write_json(MANIFEST, {
        "schema_version": "a4-c1wr-p2pf-manifest-v1",
        "status": status,
        "semantics": SEMANTICS,
        "source_sha256": source_hashes,
        "artifact_sha256": {
            str(path.relative_to(ROOT)): _sha256(path) for path in artifacts},
    })
    print(json.dumps({
        "status": status,
        "roots_deg": {
            branch: reference_rows[branch]["elevation_deg"]
            for branch in BRANCHES},
        "gates": gates,
        "timings": report["timings"],
    }, ensure_ascii=True))
    return 0 if status.startswith("passed") else 1


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    started_all = time.perf_counter()
    log_lines = []
    try:
        _verify_frozen_manifest(
            C1WR_MANIFEST, "identified_and_repaired_grid_derivative")
        if _sha256(INPUT) != (
                "0a978478b52f5363c93cbc77b0ee1868a4f50e488c82c829eb7794f0cf1235a7"):
            raise RuntimeError("refined 1 km input hash mismatch")
        prep5 = json.loads(PREP5_REPORT.read_text(encoding="utf-8"))
        p0c = json.loads(P0C_REFERENCE.read_text(encoding="utf-8"))
        if prep5.get("status") != "passed_1km_grid_stability_and_signatures":
            raise RuntimeError("Prep5 is not frozen")
        if p0c.get("status") != "p0c_1_passed":
            raise RuntimeError("P0C reference is not frozen")
        context = _load_cubic_context()
    except Exception as error:
        status = "input_or_code_drift_blocked"
        _write_json(REPORT, {
            "schema_version": "a4-c1wr-p2pf-v1",
            "status": status,
            "error": f"{type(error).__name__}: {error}",
        })
        LOG.write_text(str(error) + "\n", encoding="utf-8")
        return 2

    rows: list[dict] = []
    rays: dict[tuple[str, str, str], dict] = {}
    brent_diagnostics = {}
    try:
        for profile in PROFILES:
            for branch, bracket in BRANCHES.items():
                record, ray, diagnostics = _solve_brent(
                    context, branch, bracket, profile)
                rows.append(record)
                rays[(branch, profile, "bounded_Brent")] = ray
                brent_diagnostics[f"{profile}_{branch}"] = diagnostics
                log_lines.append(
                    f"{profile} {branch} Brent={record['elevation_deg']:.9f} "
                    f"residual={record['rx_enu_error_km']:.6g} km")
    except RuntimeError as error:
        status = "input_or_code_drift_blocked"
        _write_json(REPORT, {
            "schema_version": "a4-c1wr-p2pf-v1",
            "status": status,
            "semantics": SEMANTICS,
            "error": str(error),
            "brent_diagnostics": brent_diagnostics,
        })
        LOG.write_text("\n".join(log_lines + [str(error)]) + "\n",
                       encoding="utf-8")
        return 2

    trf_runtime = {}
    for branch in BRANCHES:
        seed = next(
            row["elevation_deg"] for row in rows
            if row["branch"] == branch and row["profile"] == "reference")
        record, ray, elapsed = _solve_trf(context, branch, seed)
        trf_runtime[branch] = elapsed
        if record is not None:
            rows.append(record)
            rays[(branch, "reference", "bounded_TRF_2D")] = ray
            log_lines.append(
                f"reference {branch} TRF={record['elevation_deg']:.9f} "
                f"residual={record['rx_enu_error_km']:.6g} km")
        else:
            log_lines.append(f"reference {branch} TRF failed after {elapsed:.3f} s")

    legacy_rows = {}
    legacy_rays = {}
    for branch, angle in LEGACY_ROOTS.items():
        started = time.perf_counter()
        ray = _trace(context, angle, "reference")
        legacy_rays[branch] = ray
        legacy_rows[branch] = _ray_record(
            context, ray, branch=branch, profile="reference",
            solver="legacy_coordinate_reference_only", elevation_deg=angle,
            elapsed_s=time.perf_counter() - started)

    by_key = {(row["branch"], row["profile"], row["solver"]): row for row in rows}
    reference_rows = {
        branch: by_key[(branch, "reference", "bounded_Brent")]
        for branch in BRANCHES
    }
    research_rows = {
        branch: by_key[(branch, "research", "bounded_Brent")]
        for branch in BRANCHES
    }
    trf_rows = {
        branch: by_key.get((branch, "reference", "bounded_TRF_2D"))
        for branch in BRANCHES
    }

    profile_differences = {}
    solver_differences = {}
    for branch in BRANCHES:
        reference = reference_rows[branch]
        research = research_rows[branch]
        profile_differences[branch] = {
            "elevation_abs_deg": abs(
                reference["elevation_deg"] - research["elevation_deg"]),
            "apex_abs_km": abs(reference["apex_km"] - research["apex_km"]),
            "tau_relative_fraction": abs(
                reference["tau_ms"] - research["tau_ms"])
                / reference["tau_ms"],
        }
        trf = trf_rows[branch]
        solver_differences[branch] = None if trf is None else {
            "elevation_abs_deg": abs(
                reference["elevation_deg"] - trf["elevation_deg"]),
            "bearing_abs_deg": abs(trf["bearing_deg"] - context["bearing_deg"]),
            "rx_enu_error_km": trf["rx_enu_error_km"],
        }

    gates = {
        "two_unique_reference_brent_roots": (
            len({row["elevation_deg"] for row in reference_rows.values()}) == 2),
        "ground_single_turning": all(
            row["termination_status"] == "ground"
            and row["turning_point_count"] == 1 for row in rows),
        "high_apex_and_delay_exceed_low": (
            reference_rows["high_ray"]["apex_km"]
            > reference_rows["low_ray"]["apex_km"]
            and reference_rows["high_ray"]["tau_ms"]
            > reference_rows["low_ray"]["tau_ms"]),
        "rx_ecef_enu_error_le_1_km": all(
            row["rx_ecef_error_km"] <= 1.0
            and row["rx_enu_error_km"] <= 1.0 for row in rows),
        "hamiltonian_residual_le_1e_8": all(
            row["hamiltonian_residual_max"] <= 1.0e-8 for row in rows),
        "spatial_fd_platform_accepted": all(
            row["spatial_fd_stability_ratio_max"] <= 1.0
            and row["spatial_fd_rejection_count"] == 0 for row in rows),
        "brent_trf_angle_le_0_02_deg": all(
            value is not None and value["elevation_abs_deg"] <= 0.02
            for value in solver_differences.values()),
        "no_branch_exchange": all(
            trf_rows[branch] is not None
            and BRANCHES[branch][0] <= trf_rows[branch]["elevation_deg"]
            <= BRANCHES[branch][1] for branch in BRANCHES),
        "research_reference_elevation_le_0_02_deg": all(
            value["elevation_abs_deg"] <= 0.02
            for value in profile_differences.values()),
        "research_reference_apex_le_0_5_km": all(
            value["apex_abs_km"] <= 0.5
            for value in profile_differences.values()),
        "research_reference_delay_le_0_2_percent": all(
            value["tau_relative_fraction"] <= 0.002
            for value in profile_differences.values()),
        "roots_inside_frozen_brackets": all(
            BRANCHES[row["branch"]][0] <= row["elevation_deg"]
            <= BRANCHES[row["branch"]][1]
            for row in rows if row["solver"] != "legacy_coordinate_reference_only"),
    }
    status = (
        "passed_current_ecef_chapman_p2p_root_freeze"
        if all(gates.values())
        else "failed_current_ecef_chapman_p2p_acceptance"
    )

    trajectory_arrays = {}
    for (branch, profile, solver), ray in rays.items():
        prefix = f"{branch}_{profile}_{solver}".lower()
        trajectory_arrays[f"{prefix}_ecef_km"] = np.asarray(
            ray["trajectory_ecef_km"], dtype=float)
        trajectory_arrays[f"{prefix}_geom_path_km"] = np.asarray(
            ray["trajectory_geom_path_km"], dtype=float)
    for branch, ray in legacy_rays.items():
        trajectory_arrays[f"{branch}_legacy_ecef_km"] = np.asarray(
            ray["trajectory_ecef_km"], dtype=float)
        trajectory_arrays[f"{branch}_legacy_geom_path_km"] = np.asarray(
            ray["trajectory_geom_path_km"], dtype=float)
    np.savez_compressed(NPZ, **trajectory_arrays)
    _write_csv(rows)
    _plot_trajectories(
        context,
        {branch: rays[(branch, "reference", "bounded_Brent")]
         for branch in BRANCHES},
        legacy_rays)
    _plot_closure(rows)
    _plot_migration(reference_rows, prep5, legacy_rows)
    _plot_gates(gates)

    prep5_rows = {
        item["branch_id"]: item for item in prep5["frozen_signatures"]
        if item["solver"] == "AH3D"
    }
    migration = {
        branch: {
            "from_prep5_linear": {
                key: reference_rows[branch][key] - prep5_rows[branch][key]
                for key in ("elevation_deg", "apex_km", "geom_path_km",
                            "group_path_km", "tau_ms")
            },
            "from_legacy_p0c": {
                "elevation_deg": (
                    reference_rows[branch]["elevation_deg"] - LEGACY_ROOTS[branch]),
                "legacy_forward_rx_error_km": legacy_rows[branch][
                    "rx_enu_error_km"],
            },
        }
        for branch in BRANCHES
    }
    source_hashes = _source_hashes()
    report = {
        "schema_version": "a4-c1wr-p2pf-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "semantics": SEMANTICS,
        "limitations": [
            "not WGS84 or IGRF production validation",
            "not full two-dimensional angular-domain mode recall",
            "not O/X magnetoionic splitting or complex three-dimensional P2P",
            "legacy roots are forward-traced references only",
        ],
        "configuration": {
            "earth_radius_km": context["radius"],
            "frequency_MHz": context["frequency_MHz"],
            "B": 0,
            "collision": 0,
            "one_way": True,
            "bearing_deg": context["bearing_deg"],
            "brackets_deg": BRANCHES,
            "brent_xtol_deg": 1.0e-6,
            "trf_bearing_half_width_deg": 0.1,
            "trf_max_nfev": 6,
            "spatial_fd_steps_km": [0.125, 0.25, 0.5],
            "provider_interpolation": "cubic",
            "provider_hash": context["trace_provider_hash"],
            "prior_linear_provider_hash": context["linear_trace_provider_hash"],
            "source_trace_grid_shape": context["source_trace_grid_shape"],
            "factorized_trace_grid_shape": (
                context["factorized_trace_grid_shape"]),
            "horizontal_factorization": (
                "exact_after_bitwise_horizontal_uniformity_check"),
        },
        "prep5_seeds_deg": PREP5_SEEDS,
        "legacy_coordinate_reference_only_deg": LEGACY_ROOTS,
        "rows": rows,
        "legacy_forward_rows": legacy_rows,
        "brent_diagnostics": brent_diagnostics,
        "profile_differences": profile_differences,
        "solver_differences": solver_differences,
        "migration": migration,
        "gates": gates,
        "timings": {
            "brent_forward_integration_s": sum(
                item["forward_integration_elapsed_s"]
                for item in brent_diagnostics.values()),
            "trf_total_s": sum(trf_runtime.values()),
            "total_s": time.perf_counter() - started_all,
        },
        "source_sha256": source_hashes,
        "figures": [
            _image_record(path) for path in (
                FIG_TRAJECTORIES, FIG_CLOSURE, FIG_MIGRATION, FIG_GATES)],
    }
    _write_json(REPORT, report)
    LOG.write_text("\n".join(log_lines + [
        f"status={status}",
        f"total_s={report['timings']['total_s']:.6f}",
    ]) + "\n", encoding="utf-8")
    artifacts = [REPORT, CSV, NPZ, LOG, FIG_TRAJECTORIES, FIG_CLOSURE,
                 FIG_MIGRATION, FIG_GATES]
    _write_json(MANIFEST, {
        "schema_version": "a4-c1wr-p2pf-manifest-v1",
        "status": status,
        "semantics": SEMANTICS,
        "source_sha256": source_hashes,
        "artifact_sha256": {
            str(path.relative_to(ROOT)): _sha256(path) for path in artifacts},
    })
    print(json.dumps({
        "status": status,
        "roots_deg": {
            branch: reference_rows[branch]["elevation_deg"]
            for branch in BRANCHES},
        "gates": gates,
        "timings": report["timings"],
    }, ensure_ascii=True))
    return 0 if status.startswith("passed") else 1


if __name__ == "__main__":
    if sys.argv[1:] == ["--finalize-existing"]:
        raise SystemExit(_finalize_existing())
    if sys.argv[1:]:
        raise SystemExit(
            "usage: _run_a4_c1wr_p2p_ecef_root_freeze.py "
            "[--finalize-existing]")
    raise SystemExit(run())
