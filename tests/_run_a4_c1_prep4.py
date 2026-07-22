"""Identify, match, deduplicate, and freeze Prep-4 Chapman ray branches."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
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
from utils import ecef_to_geocentric, geocentric_to_ecef


BASE = ROOT / "output2" / "a4_pharlap" / "c1_prep"
PREP2 = BASE / "prep2"
PREP3 = BASE / "prep3"
OUT = BASE / "prep4"
REPORT = OUT / "a4_c1_prep4_validation.json"
MANIFEST = OUT / "a4_c1_prep4_freeze_manifest.json"
SIGNATURE_CSV = OUT / "a4_c1_prep4_branch_signatures.csv"
MATCH_CSV = OUT / "a4_c1_prep4_branch_matching.csv"
TRAJECTORY_NPZ = OUT / "a4_c1_prep4_branch_trajectories.npz"
FIGURE = OUT / "a4_c1_prep4_branch_identification.png"

PREP2_REPORT = PREP2 / "a4_c1_prep2_validation.json"
PREP3_REPORT = PREP3 / "a4_c1_prep3_validation.json"
PREP2_MANIFEST = PREP2 / "a4_c1_prep2_freeze_manifest.json"
PREP3_MANIFEST = PREP3 / "a4_c1_prep3_freeze_manifest.json"
PREP2_SCAN = PREP2 / "a4_c1_prep2_scan.csv"
PREP3_SCAN = PREP3 / "a4_c1_prep3_scan.csv"
PREP2_RAW = PREP2 / "a4_c1_prep2_pharlap_raw.mat"
PREP3_NPZ = PREP3 / "a4_c1_prep3_root_trajectories.npz"
P2P_SOURCE = ROOT / "models" / "point_to_point_ah_3d_ecef.py"

DEDUP = {
    "launch_deg": 0.02,
    "trajectory_rms_km": 1.0,
    "apex_km": 1.0,
    "tau_ms": 0.05,
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
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _verify_manifest(path: Path, status: str) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != status:
        raise RuntimeError(f"unexpected frozen status: {path}")
    for item in manifest.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"frozen artifact changed: {artifact}")
    validation = manifest.get("validation_report", {})
    report = Path(validation["path"])
    if not report.is_file() or _sha256(report) != validation["sha256"]:
        raise RuntimeError(f"frozen validation report changed: {report}")
    return manifest


def _read_csv(path: Path) -> list[dict]:
    records = []
    with path.open(newline="", encoding="utf-8") as stream:
        for raw in csv.DictReader(stream):
            record = dict(raw)
            for key in (
                "elevation_deg", "ground_range_km", "along_track_residual_km",
                "cross_track_km", "apex_km", "geometric_path_km",
                "group_path_km", "tau_ms", "bracket_low_deg",
                "bracket_high_deg",
            ):
                if key in record:
                    record[key] = None if record[key] == "" else float(record[key])
            records.append(record)
    return records


def _dataset(group, name: str) -> np.ndarray:
    return np.asarray(group[name], dtype=float).reshape(-1)


def _turning_count(height_km) -> int:
    delta = np.diff(np.asarray(height_km, dtype=float))
    scale = max(float(np.max(np.abs(delta), initial=0.0)), 1.0)
    signs = np.sign(delta[np.abs(delta) > 1.0e-9 * scale])
    return int(np.count_nonzero((signs[:-1] > 0.0) & (signs[1:] < 0.0)))


def _load_pharlap_trajectories(radius_km: float) -> list[dict]:
    trajectories = []
    with h5py.File(PREP2_RAW, "r") as source:
        references = source["root_raw/ray_path_data"][()].reshape(-1)
        for reference in references:
            group = source[reference]
            latitude = _dataset(group, "lat")
            longitude = _dataset(group, "lon")
            height = _dataset(group, "height")
            ecef = geocentric_to_ecef(latitude, longitude, height, radius_km)
            path = _dataset(group, "geometric_distance")
            path -= path[0]
            trajectories.append({
                "trajectory_ecef_km": ecef,
                "trajectory_geom_path_km": path,
                "trajectory_geocentric": np.column_stack((latitude, longitude, height)),
                "turning_count": _turning_count(height),
            })
    return trajectories


def _load_ah3d_trajectories() -> list[dict]:
    trajectories = []
    with np.load(PREP3_NPZ, allow_pickle=False) as source:
        indices = sorted(int(key.split("_")[1]) for key in source.files
                         if key.endswith("_ecef_km"))
        for index in indices:
            ecef = np.asarray(source[f"candidate_{index}_ecef_km"], dtype=float)
            geocentric = np.asarray(source[f"candidate_{index}_geocentric"], dtype=float)
            path = np.asarray(source[f"candidate_{index}_geom_path_km"], dtype=float)
            trajectories.append({
                "trajectory_ecef_km": ecef,
                "trajectory_geom_path_km": path,
                "trajectory_geocentric": geocentric,
                "turning_count": _turning_count(geocentric[:, 2]),
            })
    return trajectories


def _bracket_audit(root: dict, scan: list[dict]) -> dict:
    low = float(root["bracket_low_deg"])
    high = float(root["bracket_high_deg"])

    def endpoint(angle):
        candidates = [item for item in scan
                      if item.get("status") == "ground"
                      and item.get("elevation_deg") is not None
                      and abs(item["elevation_deg"] - angle) <= 1.0e-10
                      and item.get("along_track_residual_km") is not None]
        return min(candidates, key=lambda item: abs(item["along_track_residual_km"]))

    left, right = endpoint(low), endpoint(high)
    f_left = float(left["along_track_residual_km"])
    f_right = float(right["along_track_residual_km"])
    slope = (f_right - f_left) / (high - low)
    return {
        "bracket_low_deg": low,
        "bracket_high_deg": high,
        "left_status": left["status"],
        "right_status": right["status"],
        "left_residual_km": f_left,
        "right_residual_km": f_right,
        "residual_slope_km_per_deg": slope,
        "ground_ground_sign_change": bool(f_left * f_right <= 0.0),
        "root_inside_bracket": bool(low <= float(root["elevation_deg"]) <= high),
        "nondegenerate_mapping": bool(np.isfinite(slope) and abs(slope) > 1.0),
    }


def _is_duplicate(left: dict, right: dict) -> tuple[bool, float]:
    rms = _trajectory_rms_km(left["trajectory"], right["trajectory"])
    duplicate = (
        left["status"] == right["status"]
        and left["turning_count"] == right["turning_count"]
        and abs(left["elevation_deg"] - right["elevation_deg"])
        <= DEDUP["launch_deg"]
        and abs(left["bearing_deg"] - right["bearing_deg"])
        <= DEDUP["launch_deg"]
        and rms <= DEDUP["trajectory_rms_km"]
        and abs(left["apex_km"] - right["apex_km"]) <= DEDUP["apex_km"]
        and abs(left["tau_ms"] - right["tau_ms"]) <= DEDUP["tau_ms"]
    )
    return bool(duplicate), rms


def _build_signatures(solver: str, roots: list[dict], trajectories: list[dict],
                      scan: list[dict], bearing_deg: float) -> list[dict]:
    if len(roots) != len(trajectories):
        raise RuntimeError(f"{solver} root/trajectory count mismatch")
    signatures = []
    for root, trajectory in zip(roots, trajectories):
        audit = _bracket_audit(root, scan)
        signatures.append({
            "solver": solver,
            "status": str(root["status"]),
            "turning_count": int(root.get("turning_count", trajectory["turning_count"])),
            "elevation_deg": float(root["elevation_deg"]),
            "bearing_deg": float(bearing_deg),
            "apex_km": float(root["apex_km"]),
            "geometric_path_km": float(root["geometric_path_km"]),
            "group_path_km": float(root["group_path_km"]),
            "tau_ms": float(root["tau_ms"]),
            "landing_lat_deg": float(root["landing_lat_deg"]),
            "landing_lon_deg": float(root["landing_lon_deg"]),
            "along_track_residual_km": float(root["along_track_residual_km"]),
            "cross_track_km": float(root["cross_track_km"]),
            "trajectory_ecef_sha256": _array_hash(trajectory["trajectory_ecef_km"]),
            "trajectory": trajectory,
            "bracket_audit": audit,
        })
    signatures.sort(key=lambda item: item["elevation_deg"])
    if (len(signatures) == 2 and signatures[0]["apex_km"] < signatures[1]["apex_km"]
            and all(item["status"] == "ground" for item in signatures)
            and all(item["turning_count"] == 1 for item in signatures)
            and all(all((audit["ground_ground_sign_change"],
                         audit["root_inside_bracket"],
                         audit["nondegenerate_mapping"]))
                    for audit in (item["bracket_audit"] for item in signatures))):
        signatures[0]["branch_label"] = "low_ray"
        signatures[1]["branch_label"] = "high_ray"
    return signatures


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _range_height(trajectory: dict, radius_km: float, tx_ecef: np.ndarray):
    ecef = np.asarray(trajectory["trajectory_ecef_km"], dtype=float)
    unit = ecef / np.linalg.norm(ecef, axis=1)[:, None]
    tx_unit = tx_ecef / np.linalg.norm(tx_ecef)
    ranges = radius_km * np.arccos(np.clip(unit @ tx_unit, -1.0, 1.0))
    height = np.linalg.norm(ecef, axis=1) - radius_km
    return ranges, height


def _plot(signatures: list[dict], matches: list[dict], scans: dict,
          prep1_input: Path, source_hashes: dict) -> dict:
    data = loadmat(prep1_input, squeeze_me=True, struct_as_record=False)
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    tx = np.array([float(np.asarray(data[name]).item()) for name in
                   ("origin_lat_deg", "origin_lon_deg", "origin_height_km")])
    tx_ecef = geocentric_to_ecef(*tx, radius)
    height_grid = np.asarray(data["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(data["iono_en_grid_m3"], dtype=float)
    profile = ne[ne.shape[0] // 2, ne.shape[1] // 2, :] / 1.0e11
    target = float(np.asarray(data["target_range_km"]).item())

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    colors = {"PHaRLAP": "#0072B2", "AH3D": "#D55E00"}
    for axis, solver in zip(axes[0], ("PHaRLAP", "AH3D")):
        records = [item for item in scans[solver]
                   if item.get("status") == "ground"
                   and item.get("ground_range_km") is not None]
        axis.scatter([item["elevation_deg"] for item in records],
                     [item["ground_range_km"] for item in records],
                     s=9, alpha=0.5, color=colors[solver], label="ground samples")
        for item in [value for value in signatures if value["solver"] == solver]:
            marker = "v" if item["branch_label"] == "low_ray" else "^"
            axis.scatter(item["elevation_deg"], target, marker=marker, s=130,
                         edgecolor="black", color=colors[solver], zorder=4,
                         label=item["branch_label"])
        axis.axhline(target, color="black", linestyle="--", linewidth=1.0,
                     label="RX target range")
        axis.set(title=f"{solver} independent ground map",
                 xlabel="Launch elevation (deg)", ylabel="Terminal ground range (km)")
        axis.grid(alpha=0.25)
        handles, labels = axis.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        axis.legend(unique.values(), unique.keys(), fontsize=8)

    ranges = np.linspace(0.0, target * 1.12, 350)
    mesh = axes[1, 0].pcolormesh(
        ranges, height_grid, np.broadcast_to(profile[:, None],
        (height_grid.size, ranges.size)), cmap="jet", shading="auto", alpha=0.5)
    styles = {("PHaRLAP", "low_ray"): ("#0072B2", "--"),
              ("PHaRLAP", "high_ray"): ("#0072B2", "-"),
              ("AH3D", "low_ray"): ("#D55E00", "--"),
              ("AH3D", "high_ray"): ("#D55E00", "-")}
    for item in signatures:
        ray_range, ray_height = _range_height(item["trajectory"], radius, tx_ecef)
        color, linestyle = styles[(item["solver"], item["branch_label"])]
        axes[1, 0].plot(ray_range, ray_height, color=color, linestyle=linestyle,
                        linewidth=1.7, label=f"{item['solver']} {item['branch_label']}")
    axes[1, 0].axvline(target, color="black", linestyle=":", label="RX")
    axes[1, 0].set(xlabel="Great-circle range from TX (km)",
                   ylabel="Physical height (km)", title="Topology-matched branches")
    axes[1, 0].set_ylim(0.0, 500.0)
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=7)
    fig.colorbar(mesh, ax=axes[1, 0], label="Ne (1e11 m^-3)")

    axes[1, 1].axis("off")
    table_rows = []
    for match in matches:
        table_rows.append([
            match["branch_label"],
            f"{match['pharlap_elevation_deg']:.5f}",
            f"{match['ah3d_elevation_deg']:.5f}",
            f"{match['trajectory_rms_km']:.3f}",
            f"{match['apex_abs_difference_km']:.3f}",
            f"{match['tau_abs_difference_ms']:.4f}",
        ])
    table = axes[1, 1].table(
        cellText=table_rows,
        colLabels=["branch", "PH elev", "AH elev", "traj RMS km", "|d apex| km", "|d tau| ms"],
        cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.6)
    axes[1, 1].set_title(
        "Cross-solver pairing diagnostics only\nnot A4-C1 error acceptance", pad=12)

    for label, axis in zip("ABCD", axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top",
                  ha="left", fontsize=13, fontweight="bold")
    fig.suptitle(
        "A4-C1 Prep-4 | Chapman F2 branch identification and joint dedup\n"
        "one-way | 10 MHz | B=0 O/X degenerate | topology labels | not variational Hessian or C1 validation",
        fontsize=12)
    fig.savefig(FIGURE, dpi=300, facecolor="white", metadata={
        "Title": "A4-C1 Prep-4 branch identification",
        "Description": "Topology labels and solver-specific root retention; not A4-C1 acceptance",
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
    valid = manifest.get("status") == "passed_branch_identification_and_joint_dedup"
    expected_sources = {
        "python": _sha256(Path(__file__)),
        "p2p_dedup_source": _sha256(P2P_SOURCE),
        "prep2_manifest": _sha256(PREP2_MANIFEST),
        "prep3_manifest": _sha256(PREP3_MANIFEST),
    }
    valid = valid and manifest.get("source_sha256") == expected_sources
    for item in manifest.get("frozen_artifacts", []):
        path = Path(item["path"])
        valid = valid and path.is_file() and _sha256(path) == item["sha256"]
    if valid:
        print("A4-C1-Prep-4 frozen: OK")
        return 0
    if allow_rebuild:
        return None
    print("A4-C1-Prep-4 frozen artifacts changed; rebuild refused")
    return 1


def run(*, allow_rebuild: bool = False) -> int:
    existing = _verify_existing(allow_rebuild)
    if existing is not None:
        return existing
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        _verify_manifest(PREP2_MANIFEST, "passed_pharlap_full_domain_scan")
        _verify_manifest(PREP3_MANIFEST, "passed_ah3d_full_domain_scan")
        prep2 = json.loads(PREP2_REPORT.read_text(encoding="utf-8"))
        prep3 = json.loads(PREP3_REPORT.read_text(encoding="utf-8"))
        input2 = prep2["prep1_input"]
        input3 = prep3["prep1_input"]
        if input2 != input3 or _sha256(Path(input2["path"])) != input2["sha256"]:
            raise RuntimeError("Prep-2/Prep-3 frozen input mismatch")
        data = loadmat(input2["path"], squeeze_me=True, struct_as_record=False)
        radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
        bearing = float(np.asarray(data["bearing_deg"]).item())
        scans = {"PHaRLAP": _read_csv(PREP2_SCAN), "AH3D": _read_csv(PREP3_SCAN)}
        ph_roots = list(prep2["summary"]["roots"])
        ah_roots = list(prep3["roots"])
        ph_trajectories = _load_pharlap_trajectories(radius)
        ah_trajectories = _load_ah3d_trajectories()
        signatures = (
            _build_signatures("PHaRLAP", ph_roots, ph_trajectories, scans["PHaRLAP"], bearing)
            + _build_signatures("AH3D", ah_roots, ah_trajectories, scans["AH3D"], bearing)
        )
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_frozen_input", "error": str(error)})
        return 2

    by_solver = {solver: [item for item in signatures if item["solver"] == solver]
                 for solver in ("PHaRLAP", "AH3D")}
    dedup_audit = {}
    unique = []
    for solver, items in by_solver.items():
        pair_duplicate, pair_rms = _is_duplicate(items[0], items[1])
        dedup_audit[solver] = {
            "candidate_count": len(items), "unique_count": 1 if pair_duplicate else len(items),
            "low_high_are_duplicate": pair_duplicate,
            "low_high_trajectory_rms_km": pair_rms,
        }
        unique.extend(items[:1] if pair_duplicate else items)

    matches = []
    for label in ("low_ray", "high_ray"):
        ph = next(item for item in by_solver["PHaRLAP"] if item.get("branch_label") == label)
        ah = next(item for item in by_solver["AH3D"] if item.get("branch_label") == label)
        matches.append({
            "branch_label": label,
            "pharlap_elevation_deg": ph["elevation_deg"],
            "ah3d_elevation_deg": ah["elevation_deg"],
            "elevation_abs_difference_deg": abs(ph["elevation_deg"] - ah["elevation_deg"]),
            "trajectory_rms_km": _trajectory_rms_km(ph["trajectory"], ah["trajectory"]),
            "apex_abs_difference_km": abs(ph["apex_km"] - ah["apex_km"]),
            "group_path_abs_difference_km": abs(ph["group_path_km"] - ah["group_path_km"]),
            "tau_abs_difference_ms": abs(ph["tau_ms"] - ah["tau_ms"]),
            "same_status": ph["status"] == ah["status"] == "ground",
            "same_turning_count": ph["turning_count"] == ah["turning_count"] == 1,
            "pairing_semantics": "same_independently_identified_topological_branch_diagnostic_only",
        })
    cross_matrix = np.array([
        [_trajectory_rms_km(ph["trajectory"], ah["trajectory"])
         for ah in by_solver["AH3D"]]
        for ph in by_solver["PHaRLAP"]])
    matched_cost = float(cross_matrix[0, 0] + cross_matrix[1, 1])
    swapped_cost = float(cross_matrix[0, 1] + cross_matrix[1, 0])

    signature_rows = []
    trajectory_payload = {}
    for item in unique:
        signature_rows.append({key: item[key] for key in (
            "solver", "branch_label", "status", "turning_count", "elevation_deg",
            "bearing_deg", "apex_km", "geometric_path_km", "group_path_km",
            "tau_ms", "landing_lat_deg", "landing_lon_deg",
            "along_track_residual_km", "cross_track_km", "trajectory_ecef_sha256")})
        prefix = f"{item['solver'].lower()}_{item['branch_label']}"
        trajectory_payload[f"{prefix}_ecef_km"] = item["trajectory"]["trajectory_ecef_km"]
        trajectory_payload[f"{prefix}_geom_path_km"] = item["trajectory"]["trajectory_geom_path_km"]
    _write_csv(SIGNATURE_CSV, signature_rows, list(signature_rows[0]))
    _write_csv(MATCH_CSV, matches, list(matches[0]))
    np.savez_compressed(TRAJECTORY_NPZ, **trajectory_payload)
    source_hashes = {
        "prep1_input": input2["sha256"],
        "prep2_report": _sha256(PREP2_REPORT),
        "prep3_report": _sha256(PREP3_REPORT),
        "signatures": _sha256(SIGNATURE_CSV),
        "matching": _sha256(MATCH_CSV),
        "trajectories": _sha256(TRAJECTORY_NPZ),
    }
    figure = _plot(signatures, matches, scans, Path(input2["path"]), source_hashes)

    labels = {solver: [item.get("branch_label") for item in items]
              for solver, items in by_solver.items()}
    checks = {
        "prep2_frozen": True,
        "prep3_frozen": True,
        "same_prep1_input": input2 == input3,
        "all_discovered_roots_retained": len(signature_rows) == len(ph_roots) + len(ah_roots),
        "exactly_two_roots_per_solver_for_current_case": all(
            len(items) == 2 for items in by_solver.values()),
        "all_roots_ground": all(item["status"] == "ground" for item in signatures),
        "all_roots_single_turning": all(item["turning_count"] == 1 for item in signatures),
        "all_brackets_continuous_ground_sign_changes": all(
            all((item["bracket_audit"]["ground_ground_sign_change"],
                 item["bracket_audit"]["root_inside_bracket"],
                 item["bracket_audit"]["nondegenerate_mapping"]))
            for item in signatures),
        "independent_low_high_topology": all(
            value == ["low_ray", "high_ray"] for value in labels.values()),
        "joint_signature_dedup_keeps_two_per_solver": all(
            audit["unique_count"] == 2 and not audit["low_high_are_duplicate"]
            for audit in dedup_audit.values()),
        "same_label_pairing_preferred_over_swapped": matched_cost < swapped_cost,
        "cross_solver_pairs_same_status_and_turning": all(
            item["same_status"] and item["same_turning_count"] for item in matches),
        "trajectory_signatures_complete": all(
            len(item["trajectory_ecef_sha256"]) == 64 for item in signatures),
        "direct_variational_hessian_not_claimed": True,
        "c1_error_acceptance_not_applied": True,
        "visualization": all(figure["checks"].values()),
    }
    status = ("passed_branch_identification_and_joint_dedup"
              if all(checks.values()) else "failed")
    report = {
        "schema_version": "a4-c1-prep4-validation-v1",
        "status": status,
        "case_semantics": (
            "topology_consistent_low_high_labels_not_direct_variational_hessian_"
            "or_a4_c1_external_error_acceptance"),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "theory_semantics": {
            "branch_identification": (
                "lower_launch_and_lower_apex_continuous_ground_branch=low_ray; "
                "higher_launch_and_higher_apex_continuous_ground_branch=high_ray"),
            "literature_context": (
                "Coleman (2011): high-ray optical-path minimum and low-ray saddle; "
                "Nosikov et al. (2017, 2020): systematic multi-branch identification"),
            "scope_limit": (
                "No optical-path Hessian is evaluated in Prep-4; labels are topology-"
                "consistent and literature-guided, pending any direct variational audit"),
        },
        "dedup_thresholds": {**DEDUP,
            "source": "models.point_to_point_ah_3d_ecef.find_seeded_rays_ah3d_ecef"},
        "dedup_audit": dedup_audit,
        "cross_solver_trajectory_rms_matrix_km": cross_matrix.tolist(),
        "same_label_pairing_cost_km": matched_cost,
        "swapped_pairing_cost_km": swapped_cost,
        "signatures": [{key: value for key, value in item.items()
                        if key != "trajectory"} for item in signatures],
        "matches": matches,
        "extra_root_inventory": {"PHaRLAP": [], "AH3D": []},
        "hashes": {
            "python": _sha256(Path(__file__)),
            "p2p_dedup_source": _sha256(P2P_SOURCE),
            **source_hashes,
        },
        "artifacts": {
            "signature_csv": str(SIGNATURE_CSV),
            "matching_csv": str(MATCH_CSV),
            "trajectory_npz": str(TRAJECTORY_NPZ),
            "figure": figure,
        },
        "limitations": [
            "Prep-4 does not apply A4-C1 cross-solver numerical error thresholds",
            "Low/high labels are topological, not a direct optical-path Hessian test",
            "B=0 Chapman branches do not validate O/X splitting or complex 3-D mode completeness",
        ],
    }
    _atomic_json(REPORT, report)
    artifacts = [SIGNATURE_CSV, MATCH_CSV, TRAJECTORY_NPZ, FIGURE, REPORT]
    manifest = {
        "schema_version": "a4-c1-prep4-freeze-v1",
        "status": status,
        "case_semantics": report["case_semantics"],
        "source_sha256": {
            "python": _sha256(Path(__file__)),
            "p2p_dedup_source": _sha256(P2P_SOURCE),
            "prep2_manifest": _sha256(PREP2_MANIFEST),
            "prep3_manifest": _sha256(PREP3_MANIFEST),
        },
        "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                             for path in artifacts],
    }
    _atomic_json(MANIFEST, manifest)
    print(f"A4-C1-Prep-4: {'OK' if status.startswith('passed') else 'FAIL'} "
          f"roots={len(signature_rows)} report={REPORT}")
    return 0 if status.startswith("passed") else 1


if __name__ == "__main__":
    arguments = sys.argv[1:]
    if any(argument != "--rebuild-frozen" for argument in arguments):
        raise SystemExit("usage: _run_a4_c1_prep4.py [--rebuild-frozen]")
    raise SystemExit(run(allow_rebuild="--rebuild-frozen" in arguments))
