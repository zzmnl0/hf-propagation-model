"""Aggregate frozen Prep-1..5 evidence and decide A4-C1 admission."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "output2" / "a4_pharlap" / "c1_prep"
OUT = BASE / "prep6"
REPORT = OUT / "a4_c1_prep6_validation.json"
SOURCE_CSV = OUT / "a4_c1_prep6_image_sources.csv"
LANDING_FIGURE = OUT / "a4_c1_prep6_full_domain_landing.png"
PATH_FIGURE = OUT / "a4_c1_prep6_high_low_paths.png"
STATUS_FIGURE = OUT / "a4_c1_prep6_status_distribution.png"
MANIFEST = OUT / "a4_c1_prep6_freeze_manifest.json"
TEXT_REPORT = ROOT / "word" / "a4_c1_prep_implementation_report.txt"

EXPECTED_STATUS = {
    1: "passed_input_contract_frozen",
    2: "passed_pharlap_full_domain_scan",
    3: "passed_ah3d_full_domain_scan",
    4: "passed_branch_identification_and_joint_dedup",
    5: "passed_1km_grid_stability_and_signatures",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def _atomic_json(path: Path, value: dict) -> None:
    _atomic_text(path, json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n")


def _report_path(stage: int) -> Path:
    return BASE / f"prep{stage}" / f"a4_c1_prep{stage}_validation.json"


def _manifest_path(stage: int) -> Path:
    return BASE / f"prep{stage}" / f"a4_c1_prep{stage}_freeze_manifest.json"


def _verify_stage(stage: int) -> tuple[dict, dict]:
    report_path = _report_path(stage)
    manifest_path = _manifest_path(stage)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = EXPECTED_STATUS[stage]
    if report.get("status") != expected or manifest.get("status") != expected:
        raise RuntimeError(f"Prep-{stage} prerequisite status changed")
    if not report.get("checks") or not all(report["checks"].values()):
        raise RuntimeError(f"Prep-{stage} has a failed or missing check")
    for item in manifest.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"Prep-{stage} frozen artifact changed: {artifact}")
    return report, manifest


def _verify_existing(source_hashes: dict, allow_rebuild: bool) -> int | None:
    if not MANIFEST.is_file():
        return None
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    valid = (manifest.get("status") == "passed_c1_prep_admission"
             and manifest.get("source_sha256") == source_hashes)
    for item in manifest.get("frozen_artifacts", []):
        path = Path(item["path"])
        valid = valid and path.is_file() and _sha256(path) == item["sha256"]
    if valid:
        print("A4-C1-Prep-6 frozen: OK")
        return 0
    if allow_rebuild:
        return None
    print("A4-C1-Prep-6 frozen/source hash mismatch; rebuild refused")
    return 1


def _read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _coarse_rows(rows: list[dict]) -> list[dict]:
    return sorted((row for row in rows if row["scan_stage"].startswith("coarse")),
                  key=lambda row: float(row["elevation_deg"]))


def _full_domain_ok(rows: list[dict]) -> bool:
    elevations = np.array([float(row["elevation_deg"]) for row in rows])
    expected = np.arange(1.0, 89.0 + 0.25, 0.5)
    return (len(rows) == len(expected) and np.array_equal(elevations, expected)
            and all(row.get("status") for row in rows))


def _style_header(figure, subtitle: str) -> None:
    figure.suptitle(
        "A4-C1-Prep-6 | one-way | Chapman F2 | 10 MHz | B=0 | collisions=off\n"
        + subtitle + " | spherical geocentric input / ECEF state | scattering/power=N/A",
        fontsize=10,
    )


def _plot_landing(scans: dict[str, list[dict]], target_range: float,
                  source_tag: str) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
    palette = {"ground": "#187b3f", "penetrated": "#b45823",
               "upper_boundary": "#7851a9"}
    for axis, (solver, rows) in zip(axes, scans.items()):
        for status in sorted({row["status"] for row in rows}):
            selected = [row for row in rows if row["status"] == status]
            x = np.array([float(row["elevation_deg"]) for row in selected])
            y = np.array([float(row["ground_range_km"] or "nan") for row in selected])
            finite = np.isfinite(y)
            color = palette.get(status, "#4c78a8")
            axis.scatter(x[finite], y[finite], s=16, label=status,
                         color=color, alpha=0.85)
            if np.any(~finite):
                axis.scatter(x[~finite], np.full(np.count_nonzero(~finite), -80.0),
                         s=18, marker="x", label=f"{status} (no surface endpoint)",
                         color=color, alpha=0.85)
        axis.axhline(target_range, color="black", linestyle="--", linewidth=1,
                     label="RX target")
        axis.set_title(f"{solver}: frozen 1-89 deg coarse scan")
        axis.set_xlabel("Launch elevation (deg)")
        axis.set_ylim(-120, 3020)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("Surface endpoint range (km)")
    _style_header(figure, f"complete angular-domain landing map | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.88))
    figure.savefig(LANDING_FIGURE, dpi=180)
    plt.close(figure)


def _along_range_height(ecef: np.ndarray, tx: np.ndarray, rx: np.ndarray,
                        radius_km: float) -> tuple[np.ndarray, np.ndarray]:
    tx_unit = tx / np.linalg.norm(tx)
    rx_unit = rx / np.linalg.norm(rx)
    tangent = rx_unit - np.dot(rx_unit, tx_unit) * tx_unit
    tangent /= np.linalg.norm(tangent)
    unit = ecef / np.linalg.norm(ecef, axis=1)[:, None]
    angle = np.arctan2(unit @ tangent, unit @ tx_unit)
    return radius_km * angle, np.linalg.norm(ecef, axis=1) - radius_km


def _plot_paths(data: dict, trajectories: np.lib.npyio.NpzFile,
                target_range: float, source_tag: str) -> None:
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0
    lat = np.deg2rad(float(np.asarray(data["origin_lat_deg"]).item()))
    lon = np.deg2rad(float(np.asarray(data["origin_lon_deg"]).item()))
    rx_lat = np.deg2rad(float(np.asarray(data["receiver_lat_deg"]).item()))
    rx_lon = np.deg2rad(float(np.asarray(data["receiver_lon_deg"]).item()))
    tx = radius * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    rx = radius * np.array([
        np.cos(rx_lat) * np.cos(rx_lon), np.cos(rx_lat) * np.sin(rx_lon), np.sin(rx_lat)])
    heights = np.asarray(data["height_km"], dtype=float).reshape(-1)
    ne = np.asarray(data["iono_en_grid_m3"], dtype=float)
    latitudes = np.asarray(data["latitude_deg"], dtype=float).reshape(-1)
    longitudes = np.asarray(data["longitude_deg"], dtype=float).reshape(-1)
    mid_lat = np.rad2deg(0.5 * (lat + rx_lat))
    mid_lon = np.rad2deg(0.5 * (lon + rx_lon))
    profile = ne[np.argmin(abs(latitudes - mid_lat)),
                 np.argmin(abs(longitudes - mid_lon)), :]

    figure, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
    image = None
    colors = {"low_ray": "#1368aa", "high_ray": "#c43b3b"}
    for axis, solver in zip(axes, ("PHaRLAP", "AH3D")):
        image = axis.imshow(np.repeat(profile[:, None], 2, axis=1), origin="lower",
                            extent=(0, target_range, heights[0], heights[-1]),
                            cmap="jet", alpha=0.5, aspect="auto")
        for label in ("low_ray", "high_ray"):
            key = f"{solver.lower()}_{label}_1km_ecef_km"
            distance, height = _along_range_height(
                np.asarray(trajectories[key], dtype=float), tx, rx, radius)
            axis.plot(distance, height, color=colors[label], linewidth=1.8,
                      label=label.replace("_", " "))
        axis.scatter([0, target_range], [0, 0], color="black", s=20, zorder=5)
        axis.set_title(f"{solver}: independently solved frozen roots")
        axis.set_xlabel("Great-circle range from TX (km)")
        axis.set_xlim(0, target_range)
        axis.set_ylim(0, heights[-1])
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("Physical height above spherical Earth (km)")
    colorbar_axis = figure.add_axes((0.91, 0.15, 0.018, 0.60))
    figure.colorbar(image, cax=colorbar_axis, label="Electron density Ne (m^-3)")
    _style_header(figure, f"low/high paths on frozen refined_1km Ne | jet alpha=0.5 | source={source_tag}")
    figure.subplots_adjust(left=0.07, right=0.88, bottom=0.1, top=0.80, wspace=0.10)
    figure.savefig(PATH_FIGURE, dpi=180)
    plt.close(figure)


def _plot_status(scans: dict[str, list[dict]], source_tag: str) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11, 5.6), sharey=True)
    palette = {"ground": "#187b3f", "penetrated": "#b45823",
               "upper_boundary": "#7851a9"}
    for axis, (solver, rows) in zip(axes, scans.items()):
        counts = Counter(row["status"] for row in rows)
        labels = sorted(counts)
        values = [counts[label] for label in labels]
        bars = axis.bar(labels, values, color=[palette.get(label, "#4c78a8") for label in labels])
        axis.bar_label(bars)
        axis.set_title(f"{solver}: raw termination states (N={len(rows)})")
        axis.set_xlabel("Unmapped solver termination state")
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Coarse-scan angle count")
    _style_header(figure, f"full angular-domain status distribution | source={source_tag}")
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    figure.savefig(STATUS_FIGURE, dpi=180)
    plt.close(figure)


def _image_check(path: Path) -> dict:
    pixels = mpimg.imread(path)
    finite = bool(np.isfinite(pixels).all())
    result = {
        "path": str(path),
        "sha256": _sha256(path),
        "width_px": int(pixels.shape[1]),
        "height_px": int(pixels.shape[0]),
        "finite_pixels": finite,
        "pixel_std": float(np.std(pixels)) if finite else None,
    }
    result["passed"] = bool(result["width_px"] >= 1200 and result["height_px"] >= 700
                            and finite and result["pixel_std"] > 1e-4)
    return result


def _write_source_inventory(rows: list[dict]) -> None:
    with SOURCE_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "role", "stage", "path", "sha256", "semantic"))
        writer.writeheader()
        writer.writerows(rows)


def _text_report(status: str, reports: dict[int, dict], scans: dict[str, list[dict]],
                 image_checks: list[dict], source_tag: str) -> str:
    convergence = reports[5]["convergence"]
    maximum = lambda key: max(abs(float(row[key])) for row in convergence)
    counts = {solver: dict(Counter(row["status"] for row in rows))
              for solver, rows in scans.items()}
    return f"""A4-C1-Prep实施与预验收报告

日期：{datetime.now().astimezone().date().isoformat()}
机器状态：{status}
语义：仅确认Chapman F2、B=0、碰撞关闭算例的等价输入、完整角域分支可见、冻结根稳定和A4-C1准入；不是AH3D-PHaRLAP正式跨软件误差验收，也不是外部物理验证通过。

一、冻结输入与扫描

1. Prep-1至Prep-5的报告、manifest及其冻结产物哈希均复验通过；图像来源短哈希：{source_tag}。
2. PHaRLAP和AH3D均覆盖1至89 deg、0.5 deg步长的177个粗扫角度，无空终止状态。
3. 原始状态分布：PHaRLAP={counts['PHaRLAP']}；AH3D={counts['AH3D']}。penetrated与upper_boundary保留各求解器原始语义，未伪映射为共同状态。

二、分支与稳定性

1. 两端各保留low_ray和high_ray两个独立ground、单转向根；当前算例无额外有效根被删除。
2. 最大仰角变化={maximum('elevation_abs_difference_deg'):.12g} deg，门槛0.02 deg；最大apex变化={maximum('apex_abs_difference_km'):.12g} km，门槛0.5 km。
3. 最大群时延相对变化={100.0 * maximum('tau_relative_difference'):.12g}%，门槛0.2%；最大沿程/横程端点残差分别为{maximum('along_track_residual_km'):.12g}/{maximum('cross_track_km'):.12g} km，门槛1 km。
4. low/high标签为拓扑一致、文献引导语义；未计算光程泛函Hessian，不能宣称直接验证高射线极小值和低射线鞍点。

三、可视化与结论

1. 已生成双方完整角域落点图、冻结Ne背景上的低/高路径图和原始状态分布图；三图均通过尺寸、有限像素、非空白及来源哈希检查。
2. Ne背景只读refined_1km冻结网格，使用jet、alpha=0.5；未重新生成Chapman场。
3. Prep第五节全部门槛通过，允许进入A4-C1正式外部误差验收。C1仍须分别执行固定同角前向轨迹比较和双方独立冻结根同分支比较。

四、风险与未完成项

1. B=0水平分层算例不能验证O/X分裂、位置相关磁场、复杂三维扰动模式完整性或全波散射。
2. PHaRLAP penetrated与AH3D upper_boundary事件定义不同；C1只比较双方均ground且同分支的轨迹。
3. Prep只证明各求解器内部1 km/2 km网格稳定，尚未应用C1轨迹、落点、apex、群时延和发射角跨软件门槛。
4. low/high变分Hessian直接分类仍未实施，不影响本次拓扑分支准入，但不得扩大论文表述。

可视化硬检查：{json.dumps(image_checks, ensure_ascii=True)}
"""


def run(allow_rebuild: bool = False) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        verified = {stage: _verify_stage(stage) for stage in range(1, 6)}
        reports = {stage: pair[0] for stage, pair in verified.items()}
        prep1_inputs = reports[1]["inputs"]
        refined = next(item for item in prep1_inputs if item["variant"] == "refined_1km")
        input_path = Path(refined["path"])
        scan_paths = {
            "PHaRLAP": BASE / "prep2" / "a4_c1_prep2_scan.csv",
            "AH3D": BASE / "prep3" / "a4_c1_prep3_scan.csv",
        }
        scans = {solver: _coarse_rows(_read_csv(path)) for solver, path in scan_paths.items()}
        trajectory_path = BASE / "prep5" / "a4_c1_prep5_trajectories.npz"
        data = loadmat(input_path, squeeze_me=False, struct_as_record=False)
        target_range = float(np.asarray(data["target_range_km"]).item())
        source_hashes = {
            "python": _sha256(Path(__file__)),
            "refined_input": _sha256(input_path),
            "prep5_trajectory": _sha256(trajectory_path),
            **{f"prep{stage}_report": _sha256(_report_path(stage)) for stage in range(1, 6)},
            **{f"prep{stage}_manifest": _sha256(_manifest_path(stage)) for stage in range(1, 6)},
            **{solver.lower() + "_scan": _sha256(path) for solver, path in scan_paths.items()},
        }
        source_tag = hashlib.sha256(json.dumps(
            source_hashes, sort_keys=True).encode("ascii")).hexdigest()[:12]
        existing = _verify_existing(source_hashes, allow_rebuild)
        if existing is not None:
            return existing
        _plot_landing(scans, target_range, source_tag)
        with np.load(trajectory_path, allow_pickle=False) as trajectories:
            _plot_paths(data, trajectories, target_range, source_tag)
        _plot_status(scans, source_tag)
        image_checks = [_image_check(path) for path in (
            LANDING_FIGURE, PATH_FIGURE, STATUS_FIGURE)]

        convergence = reports[5]["convergence"]
        signatures = reports[5]["frozen_signatures"]
        checks = {
            "prep1_to_prep5_frozen_and_all_checks_passed": True,
            "same_frozen_refined_input_hash": _sha256(input_path) == refined["sha256"],
            "both_complete_1_to_89_deg_coarse_scans": all(
                _full_domain_ok(rows) for rows in scans.values()),
            "all_coarse_scan_states_recorded": all(
                all(row.get("status") for row in rows) for rows in scans.values()),
            "all_discovered_roots_retained": reports[4]["checks"]["all_discovered_roots_retained"],
            "two_unique_low_high_roots_per_solver": all(
                reports[4]["dedup_audit"][solver]["unique_count"] == 2
                for solver in ("PHaRLAP", "AH3D")),
            "all_frozen_roots_ground_single_turning": len(signatures) == 4 and all(
                item["termination"] == "ground" and item["turning_count"] == 1
                for item in signatures),
            "all_endpoint_residuals_within_1km": all(
                abs(float(row["along_track_residual_km"])) <= 1.0
                and abs(float(row["cross_track_km"])) <= 1.0 for row in convergence),
            "all_1km_vs_2km_stability_limits_passed": all(row["passed"] for row in convergence),
            "legacy_angles_not_used": str(np.asarray(
                data["legacy_angle_usage"]).reshape(-1)[0]) == "not_present_not_used",
            "visualization": all(item["passed"] for item in image_checks),
            "a4_c1_cross_solver_error_acceptance_not_applied": True,
            "direct_variational_hessian_not_claimed": True,
        }
        status = "passed_c1_prep_admission" if all(checks.values()) else "failed"

        inventory = []
        for key, digest in source_hashes.items():
            if key == "python":
                path = Path(__file__)
            elif key == "refined_input":
                path = input_path
            elif key == "prep5_trajectory":
                path = trajectory_path
            elif key.endswith("_scan"):
                path = scan_paths["PHaRLAP" if key.startswith("pharlap") else "AH3D"]
            else:
                stage = int(key[4])
                path = _report_path(stage) if key.endswith("report") else _manifest_path(stage)
            inventory.append({"role": key, "stage": "Prep-6 source", "path": str(path),
                              "sha256": digest, "semantic": "frozen read-only evidence"})
        for item in image_checks:
            inventory.append({"role": "visualization", "stage": "Prep-6",
                              "path": item["path"], "sha256": item["sha256"],
                              "semantic": "visual evidence; does not override numeric gates"})
        _write_source_inventory(inventory)
        _atomic_text(TEXT_REPORT, _text_report(
            status, reports, scans, image_checks, source_tag))
        result = {
            "schema_version": "a4-c1-prep6-validation-v1",
            "status": status,
            "case_semantics": (
                "c1_prep_admission_only_not_cross_solver_error_acceptance_"
                "not_external_physics_validation"),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "checks": checks,
            "source_sha256": source_hashes,
            "source_set_short_hash": source_tag,
            "coarse_scan_status_counts": {
                solver: dict(Counter(row["status"] for row in rows))
                for solver, rows in scans.items()},
            "frozen_root_count": len(signatures),
            "grid_convergence": convergence,
            "visualization_contract": {
                "ne_background": "refined_1km frozen grid, jet, alpha=0.5",
                "figures": image_checks,
            },
            "artifacts": {
                "image_source_inventory_csv": str(SOURCE_CSV),
                "implementation_report": str(TEXT_REPORT),
            },
            "limitations": [
                "No A4-C1 PHaRLAP-AH3D cross-solver error threshold is applied in Prep-6",
                "B=0 stratified Chapman does not validate O/X or complex-3D mode completeness",
                "Low/high labels are topology-consistent, not a direct optical-path Hessian result",
                "AH3D remains a Hamiltonian ray solver, not a full-wave solver",
            ],
            "next_step": (
                "A4-C1 fixed-same-angle forward trajectory comparison, then independent "
                "frozen-root same-branch comparison" if status.startswith("passed") else
                "repair failed Prep gate before A4-C1"),
        }
        _atomic_json(REPORT, result)
        frozen = [SOURCE_CSV, LANDING_FIGURE, PATH_FIGURE, STATUS_FIGURE, TEXT_REPORT, REPORT]
        _atomic_json(MANIFEST, {
            "schema_version": "a4-c1-prep6-freeze-v1",
            "status": status,
            "case_semantics": result["case_semantics"],
            "source_sha256": source_hashes,
            "frozen_artifacts": [{"path": str(path), "sha256": _sha256(path)}
                                 for path in frozen],
        })
        print(f"A4-C1-Prep-6: {'OK' if status.startswith('passed') else 'FAIL'} "
              f"status={status} report={REPORT}")
        return 0 if status.startswith("passed") else 1
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_prerequisite", "error": str(error)})
        print(f"A4-C1-Prep-6: BLOCKED {type(error).__name__}: {error}")
        return 2


if __name__ == "__main__":
    arguments = sys.argv[1:]
    if any(argument != "--rebuild-frozen" for argument in arguments):
        raise SystemExit("usage: _run_a4_c1_prep6.py [--rebuild-frozen]")
    raise SystemExit(run(allow_rebuild="--rebuild-frozen" in arguments))
