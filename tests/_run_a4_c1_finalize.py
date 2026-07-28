"""Finalize and freeze the failed A4-C1 Chapman external comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import (
    _atomic_json,
    _image_check,
    _load_pharlap_case,
    _refined_input,
    _sha256,
    _trajectory_metrics,
    _verify_manifest,
)


C1 = ROOT / "output2" / "a4_pharlap" / "c1"
OUT = C1 / "final"
WORD = ROOT / "word"
REPORT_TXT = WORD / "a4_c1_chapman_external_validation_report.txt"
VALIDATION = OUT / "a4_c1_validation.json"
METRICS_CSV = OUT / "a4_c1_metrics.csv"
TRAJECTORY_NPZ = OUT / "a4_c1_trajectory_comparison.npz"
SOURCE_LIST = OUT / "a4_c1_source_inventory.json"
IMAGE_LIST = OUT / "a4_c1_image_source_list.json"
SUMMARY_FIGURE = OUT / "a4_c1_validation_summary.png"
PROVENANCE_FIGURE = OUT / "a4_c1_provenance_coverage.png"
DECISION_FIGURE = OUT / "a4_c1_scope_decision.png"
MANIFEST = OUT / "a4_c1_freeze_manifest.json"

C1_0_MANIFEST = C1 / "admission" / "c1_0_freeze_manifest.json"
C1_1_REPORT = C1 / "same_angle" / "c1_1_validation.json"
C1_1_MANIFEST = C1 / "same_angle" / "c1_1_freeze_manifest.json"
C1_2_REPORT = C1 / "independent_roots" / "c1_2_validation.json"
C1_2_MANIFEST = C1 / "independent_roots" / "c1_2_freeze_manifest.json"
C1_3_REPORT = C1 / "attribution" / "c1_3_validation.json"
C1_3_MANIFEST = C1 / "attribution" / "c1_3_freeze_manifest.json"
C1R5_REPORT = C1 / "c1r" / "c1r_5" / "c1r_5_decision.json"
C1R5_MANIFEST = C1 / "c1r" / "c1r_5" / "c1r_5_freeze_manifest.json"
GATE_CSV = C1 / "attribution" / "c1_acceptance_gate_details.csv"
INDEPENDENT_RESIDUAL_NPZ = (
    C1 / "independent_roots" / "c1_independent_root_residuals.npz"
)

FINAL_STATUS = "failed_c1_chapman_high_low_external_comparison"
EXPECTED_FAILURE_GATE = "trajectory_rms_3d_km"
EXPECTED_C1R_STATUS = "completed_c1r5_software_difference_or_unresolved"


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value: dict) -> None:
    _atomic_json(path, value)


def _metric_rows(c1_1: dict, c1_2: dict) -> list[dict]:
    fields = (
        "case_id", "branch_id", "anchor_solver", "elevation_deg",
        "pharlap_elevation_deg", "ah3d_elevation_deg",
        "launch_elevation_absolute_error_deg",
        "landing_absolute_error_km", "landing_relative_error",
        "apex_height_absolute_error_km", "geometric_path_relative_error",
        "group_path_relative_error", "group_delay_relative_error",
        "trajectory_rms_3d_km", "trajectory_max_3d_km",
        "trajectory_horizontal_rms_km", "trajectory_up_rms_km",
        "pharlap_along_track_residual_km", "pharlap_cross_track_residual_km",
        "ah3d_along_track_residual_km", "ah3d_cross_track_residual_km",
    )
    rows = []
    for stage, report in (("C1-1", c1_1), ("C1-2", c1_2)):
        for metric in report["metrics"]:
            row = {
                "stage": stage,
                "comparison_semantics": metric.get(
                    "comparison_semantics",
                    "same_launch_angle_forward_no_root_solving",
                ),
                **{field: metric.get(field) for field in fields},
                "passed": bool(metric.get(
                    "passed",
                    all(value for value in metric["gates"].values()
                        if value is not None),
                )),
                "gates_json": json.dumps(metric["gates"], sort_keys=True),
            }
            rows.append(row)
    return rows


def _write_metrics(rows: list[dict]) -> None:
    with METRICS_CSV.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_trajectory_bundle() -> None:
    input_path, _ = _refined_input()
    from scipy.io import loadmat
    input_data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    radius_km = float(np.asarray(input_data["earth_radius_m"]).item()) / 1000.0
    ph_summary = json.loads(
        (C1 / "same_angle" / "raw" / "c1_same_angle_pharlap_summary.json")
        .read_text(encoding="utf-8")
    )
    ph_cases = {item["case_id"]: item for item in ph_summary["cases"]}
    arrays = {
        "semantics": np.array(
            "normalized_path_cross_solver_residuals; no registration; no root solving"
        )
    }
    for case_id in (
        "low_pharlap_anchor", "low_ah3d_anchor",
        "high_pharlap_anchor", "high_ah3d_anchor",
    ):
        ph_path = _load_pharlap_case(Path(ph_cases[case_id]["raw_mat"]), radius_km)[0]
        with np.load(
            C1 / "same_angle" / "raw" / f"ah3d_{case_id}_raw.npz",
            allow_pickle=False,
        ) as source:
            ah_path = {
                "trajectory_ecef_km": np.asarray(
                    source["repeat1_trajectory_ecef_km"], dtype=float
                ),
                "trajectory_geom_path_km": np.asarray(
                    source["repeat1_trajectory_geom_path_km"], dtype=float
                ),
            }
        _, residual = _trajectory_metrics(ph_path, ah_path, radius_km)
        for key, values in residual.items():
            arrays[f"same_angle_{case_id}_{key}"] = np.asarray(values)
    with np.load(INDEPENDENT_RESIDUAL_NPZ, allow_pickle=False) as source:
        for key in source.files:
            arrays[f"independent_roots_{key}"] = np.asarray(source[key])
    np.savez_compressed(TRAJECTORY_NPZ, **arrays)


def _prior_files() -> list[Path]:
    files = []
    for directory in (
        C1 / "admission", C1 / "same_angle", C1 / "independent_roots",
        C1 / "attribution", C1 / "c1r",
    ):
        files.extend(path for path in directory.rglob("*") if path.is_file())
    files.extend([
        ROOT / "tests" / "_run_a4_c1_validation.py",
        ROOT / "tests" / "_run_a4_c1_independent_roots.py",
        ROOT / "tests" / "_run_a4_c1_attribution.py",
        ROOT / "tests" / "_run_a4_c1r_isolation.py",
        ROOT / "tests" / "_verify_a4_c1_validation.py",
        ROOT / "tests" / "_verify_a4_c1_independent_roots.py",
        ROOT / "tests" / "_verify_a4_c1_attribution.py",
        ROOT / "tests" / "_verify_a4_c1r_isolation.py",
        ROOT / "tests" / "_run_a4_c1_finalize.py",
        ROOT / "tests" / "_verify_a4_c1_finalize.py",
        ROOT / "tests" / "matlab" / "a4_c1_same_angle_pharlap.m",
        WORD / "a4_c1_chapman_external_validation_implementation_plan.txt",
        WORD / "a4_pharlap_external_validation_implementation_plan.txt",
        WORD / "reference-papers.txt",
        Path(r"D:\Matlab\PHaRLAP\a4_c1_same_angle_pharlap.m"),
        Path(r"D:\pharlap\pharlap_4.7.4\mex\raytrace_3d_sp.mexw64"),
        Path(r"D:\pharlap\pharlap_4.7.4\src\matlab\raytrace_3d_sp.m"),
    ])
    input_path, _ = _refined_input()
    files.append(input_path)
    unique = {str(path.resolve()).lower(): path.resolve() for path in files}
    missing = [str(path) for path in unique.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("required frozen source missing: " + ", ".join(missing))
    return sorted(unique.values(), key=lambda path: str(path).lower())


def _write_source_inventory(paths: list[Path]) -> dict:
    categories = {
        "input": 0, "source": 0, "raw_output": 0,
        "comparison": 0, "visualization": 0, "report": 0,
    }
    items = []
    for path in paths:
        suffix = path.suffix.lower()
        if path.name == "a4_c1_chapman_refined_1km.mat":
            category = "input"
        elif suffix in {".py", ".m", ".mexw64"}:
            category = "source"
        elif "raw" in path.parts or suffix == ".mat":
            category = "raw_output"
        elif suffix == ".png":
            category = "visualization"
        elif suffix == ".txt":
            category = "report"
        else:
            category = "comparison"
        categories[category] += 1
        items.append({
            "path": str(path),
            "sha256": _sha256(path),
            "size_bytes": path.stat().st_size,
            "category": category,
        })
    payload = {
        "schema_version": "a4-c1-4-source-inventory-v1",
        "semantics": "all frozen C1/C1R inputs, sources, raw outputs and evidence",
        "category_counts": categories,
        "items": items,
    }
    _write_json(SOURCE_LIST, payload)
    return payload


def _short_hash(value: dict) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:12]


def _plot_summary(c1_1: dict, c1_2: dict, source_tag: str) -> dict:
    cases = c1_1["metrics"]
    labels = [item["case_id"].replace("_", "\n") for item in cases]
    rms = [item["trajectory_rms_3d_km"] for item in cases]
    maximum = [item["trajectory_max_3d_km"] for item in cases]
    x = np.arange(len(cases))
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.6), constrained_layout=True)
    axes[0].bar(x - 0.18, rms, width=0.36, label="trajectory RMS")
    axes[0].bar(x + 0.18, maximum, width=0.36, label="trajectory max")
    axes[0].axhline(2.0, color="#b22222", linestyle="--", label="RMS limit 2 km")
    axes[0].axhline(5.0, color="black", linestyle=":", label="max limit 5 km")
    axes[0].set_xticks(x, labels, fontsize=8)
    axes[0].set_ylabel("Cross-solver residual (km)")
    axes[0].set_title("C1-1 same-angle forward comparison")
    axes[0].grid(axis="y", alpha=0.2)
    axes[0].legend(fontsize=8)

    independent = c1_2["metrics"]
    independent_labels = [item["branch_id"] for item in independent]
    axes[1].bar(
        independent_labels,
        [item["trajectory_rms_3d_km"] for item in independent],
        color=["#1368aa", "#c43b3b"],
    )
    axes[1].axhline(2.0, color="#b22222", linestyle="--", label="RMS limit 2 km")
    axes[1].set_ylabel("Cross-solver residual (km)")
    axes[1].set_title("C1-2 independent roots (diagnostic only)")
    axes[1].grid(axis="y", alpha=0.2)
    axes[1].legend(fontsize=8)
    figure.suptitle(
        "A4-C1-4 final gate evidence | one-way | Chapman F2 | 10 MHz | "
        "B=0 | collision=0\n"
        "Spherical geocentric input / ECEF state | scattering, PE, power=N/A | "
        f"source={source_tag}"
    )
    figure.savefig(SUMMARY_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(SUMMARY_FIGURE)


def _plot_provenance(inventory: dict, source_tag: str) -> dict:
    counts = inventory["category_counts"]
    labels = list(counts)
    values = [counts[label] for label in labels]
    figure, axis = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    bars = axis.bar(labels, values, color="#4c78a8")
    axis.bar_label(bars)
    axis.set_ylabel("Frozen file count")
    axis.set_title(
        "A4-C1-4 provenance coverage | SHA-256 recorded for every listed file\n"
        f"strict no-recompute restore | source={source_tag}"
    )
    axis.grid(axis="y", alpha=0.2)
    figure.savefig(PROVENANCE_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(PROVENANCE_FIGURE)


def _plot_decision(source_tag: str) -> dict:
    stages = ("C1-0", "C1-1", "C1-2", "C1-3", "C1R", "C1-4", "C2")
    labels = ("PASS", "FAIL", "DIAG", "FAIL", "UNRESOLVED", "FROZEN FAIL", "BLOCKED")
    values = (1, -1, 0, -1, 0, -1, -1)
    colors = ("#187b3f", "#b22222", "#777777", "#b22222",
              "#777777", "#b22222", "#555555")
    figure, axis = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
    y = np.arange(len(stages))
    axis.barh(y, values, color=colors)
    axis.set_yticks(y, stages)
    axis.set_xlim(-1.25, 1.25)
    axis.axvline(0, color="black", linewidth=0.8)
    for index, (value, label) in enumerate(zip(values, labels)):
        x_text = 0.04 if value == 0 else (0.9 if value > 0 else -0.68)
        axis.text(
            x_text, index, label, va="center",
            ha="left" if value == 0 else "center",
            color="white" if value else colors[index], fontweight="bold",
        )
    axis.invert_yaxis()
    axis.set_xticks((-1, 0, 1), ("blocked/failed", "diagnostic", "passed"))
    axis.set_title(
        "A4-C1-4 decision boundary | C1 failure preserved; C2 not admitted\n"
        "B=0 smooth Chapman only; no O/X, irregularity, scattering, PE, "
        f"full-wave or power claim | source={source_tag}"
    )
    axis.grid(axis="x", alpha=0.2)
    figure.savefig(DECISION_FIGURE, dpi=180)
    plt.close(figure)
    return _image_check(DECISION_FIGURE)


def _write_report(c1_1: dict, c1_2: dict, c1r5: dict, source_tag: str) -> None:
    failures = [
        item for item in c1_1["metrics"]
        if not item["gates"][EXPECTED_FAILURE_GATE]
    ]
    low = [item for item in c1_1["metrics"] if item["branch_id"] == "low_ray"]
    high_independent = next(
        item for item in c1_2["metrics"] if item["branch_id"] == "high_ray"
    )
    text = f"""A4-C1 Chapman F2无磁高/低射线PHaRLAP外部验证报告

日期：2026-07-24
机器状态：{FINAL_STATUS}
来源集合短哈希：{source_tag}

一、验证范围

1. 本阶段比较PHaRLAP 4.7.4与AH3D在同一冻结Chapman F2电子密度网格上的单程射线结果：10 MHz、B=0、碰撞关闭、球形地心纬经高输入、ECEF内部状态。
2. C1-1固定相同发射角，直接比较两套求解器的前向轨迹；C1-2比较双方独立冻结根，只作为P2P端点约束诊断。
3. 本结论不验证O/X分裂、复杂三维扰动、散射、Ray-PE、全波或绝对功率。

二、原始输出与统一比较

1. PHaRLAP原始MAT、AH3D原始NPZ/JSON、运行日志、冻结输入和源码均保留，未对轨迹执行配准、平移、旋转或缩放。
2. 正式轨迹按累计几何路径归一化至u=0..1的201点，在ECEF中逐点直接相减；轨迹RMS门为2 km，最大误差门为5 km。
3. 低射线两个公共角的轨迹RMS为{low[0]['trajectory_rms_3d_km']:.6f} km和{low[1]['trajectory_rms_3d_km']:.6f} km，全部正式门通过。
4. 高射线PHaRLAP锚点和AH3D锚点的轨迹RMS分别为{failures[0]['trajectory_rms_3d_km']:.6f} km和{failures[1]['trajectory_rms_3d_km']:.6f} km，超过2 km门；对应最大误差为{failures[0]['trajectory_max_3d_km']:.6f} km和{failures[1]['trajectory_max_3d_km']:.6f} km，仍低于5 km门。
5. 双方独立高射线根经端点约束后的轨迹RMS为{high_independent['trajectory_rms_3d_km']:.6f} km，但该结果包含不同发射角补偿，不能覆盖同角失败。

三、差异归因

1. C1R已排除或排除为主因：冻结输入、坐标/地球模型、ground事件、DOP853积分容差、冻结场线性表示和AH3D多尺度Hamiltonian空间导数。
2. PHaRLAP MEX内部插值、ODE和独立变量参数化不可观察，未定位可复现的AH3D生产缺陷；剩余状态为software_difference_or_unresolved。
3. 未修改生产AH3D，未放宽门槛，也未用C1-2独立根结果覆盖C1-1失败。

四、门槛状态与结论

1. 正式失败项只有C1-1两组high_ray的trajectory_rms_3d_km；其他同角正式门及重复性门通过。
2. A4-C1最终状态为failed，而非blocked或passed。
3. 不允许进入A4-C2正式O/X验收；若后续运行C2，只能在另行确认后标记为诊断性结果。

五、冻结与复验

1. 机器产物位于output2/a4_pharlap/c1/final，包括a4_c1_validation.json、指标CSV、轨迹残差NPZ、来源清单、图像来源清单、三张可视化和freeze manifest。
2. tests/_verify_a4_c1_finalize.py执行严格无重算复验；任一输入、源码、MEX、原始输出、比较结果或图像SHA-256漂移均拒绝恢复。
3. 理论依据及完整论文格式沿用word/reference-papers.txt中的Haselgrove（1955）、Haselgrove与Haselgrove（1960）、Jones与Stephenson（1975）、Budden（1985）、Coleman（2011）、Nosikov等（2017、2020）及Cervera与Harris（2014）；C1-4未新增物理参数或文献结论。

六、风险与下一步

1. 当前证据不能判断AH3D或PHaRLAP哪一个更接近真实电离层；两套软件的接收机落点误差只能评价各自求根命中能力，不能代替独立真值。
2. 优先获取PHaRLAP源级算法说明/中间状态，或引入第三套可审计射线实现；在定位AH3D可复现缺陷前不修改生产求解器。
3. C1未通过，因此A4总体仍未完成，B1及后续阶段不得据此准入。
"""
    REPORT_TXT.write_text(text, encoding="utf-8")


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    _verify_manifest(C1_0_MANIFEST, "passed_c1_0_frozen_admission")
    _verify_manifest(C1_1_MANIFEST, "failed_c1_1_numeric_thresholds")
    _verify_manifest(
        C1_2_MANIFEST, "completed_c1_2_diagnostic_after_same_angle_failure"
    )
    _verify_manifest(C1_3_MANIFEST, "failed_pending_c1r")
    c1r5_manifest = json.loads(C1R5_MANIFEST.read_text(encoding="utf-8"))
    if c1r5_manifest.get("status") != EXPECTED_C1R_STATUS:
        raise RuntimeError("C1R-5 state does not authorize failure freeze")
    for item in c1r5_manifest["frozen_artifacts"]:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"C1R-5 artifact drift: {path}")

    c1_1 = json.loads(C1_1_REPORT.read_text(encoding="utf-8"))
    c1_2 = json.loads(C1_2_REPORT.read_text(encoding="utf-8"))
    c1_3 = json.loads(C1_3_REPORT.read_text(encoding="utf-8"))
    c1r5 = json.loads(C1R5_REPORT.read_text(encoding="utf-8"))
    rows = _metric_rows(c1_1, c1_2)
    _write_metrics(rows)
    _write_trajectory_bundle()
    inventory = _write_source_inventory(_prior_files())
    source_tag = _short_hash(inventory)

    visualizations = [
        _plot_summary(c1_1, c1_2, source_tag),
        _plot_provenance(inventory, source_tag),
        _plot_decision(source_tag),
    ]
    if not all(item["passed"] for item in visualizations):
        raise RuntimeError("C1-4 visualization contract failed")

    image_items = []
    for path in sorted(
        [
            item for directory in (
                C1 / "admission", C1 / "same_angle", C1 / "independent_roots",
                C1 / "attribution", C1 / "c1r", OUT,
            )
            for item in directory.rglob("*.png")
        ],
        key=lambda item: str(item).lower(),
    ):
        check = _image_check(path)
        if not check["passed"]:
            raise RuntimeError(f"invalid frozen image: {path}")
        image_items.append(check)
    _write_json(IMAGE_LIST, {
        "schema_version": "a4-c1-4-image-sources-v1",
        "semantics": "C1/C1R evidence images plus C1-4 final status images",
        "images": image_items,
    })

    failed_gates = [
        row for row in csv.DictReader(GATE_CSV.open(encoding="utf-8"))
        if row["stage"] == "C1-1" and row["status"] == "fail"
    ]
    expected_failures = (
        len(failed_gates) == 2
        and all(
            row["branch_id"] == "high_ray"
            and row["gate"] == EXPECTED_FAILURE_GATE
            for row in failed_gates
        )
    )
    checks = {
        "all_prerequisite_manifests_verified": True,
        "formal_failure_signature_preserved": expected_failures,
        "independent_roots_not_used_to_override_same_angle_failure": True,
        "c1r_no_reproducible_ah3d_defect_isolated":
            c1r5["status"] == EXPECTED_C1R_STATUS,
        "production_solver_unchanged": c1r5["checks"][
            "no_production_source_change_during_c1r"
        ],
        "metrics_csv_complete": len(rows) == 6,
        "trajectory_npz_present": TRAJECTORY_NPZ.is_file(),
        "source_inventory_complete": len(inventory["items"]) >= 1,
        "three_final_visualizations_valid": len(visualizations) == 3
            and all(item["passed"] for item in visualizations),
        "scope_limit_explicit": True,
        "c2_admission_denied": True,
    }
    if not all(checks.values()):
        raise RuntimeError("C1-4 completion checks failed")

    _write_report(c1_1, c1_2, c1r5, source_tag)
    validation = {
        "schema_version": "a4-c1-4-validation-v1",
        "status": FINAL_STATUS,
        "created_utc": _utc(),
        "case_semantics": (
            "formal_failed_freeze_of_one_way_10MHz_B0_collision0_Chapman_"
            "PHaRLAP_AH3D_comparison_not_OX_not_scattering_not_PE_not_full_wave"
        ),
        "checks": checks,
        "formal_failed_gates": failed_gates,
        "c1_passed": False,
        "c2_admitted": False,
        "difference_attribution": c1_3["primary_attribution"],
        "remaining_subcause": "software_implementation_difference_or_unresolved",
        "c1r_status": c1r5["status"],
        "source_set_short_hash": source_tag,
        "outputs": {
            "metrics_csv": str(METRICS_CSV.resolve()),
            "trajectory_npz": str(TRAJECTORY_NPZ.resolve()),
            "source_inventory": str(SOURCE_LIST.resolve()),
            "image_source_list": str(IMAGE_LIST.resolve()),
            "report_txt": str(REPORT_TXT.resolve()),
        },
        "visualizations": visualizations,
        "limitations": [
            "B=0 smooth Chapman only",
            "Does not validate O/X, complex 3-D irregularities, scattering, PE, full-wave, or power",
            "Cannot determine which solver is closer to physical truth without independent truth data",
            "PHaRLAP MEX internal interpolation and ODE implementation remain unobservable",
        ],
        "next_step": (
            "Do not start formal C2. Obtain PHaRLAP source-level/intermediate-state "
            "evidence or use a third auditable ray tracer before changing AH3D."
        ),
    }
    _write_json(VALIDATION, validation)

    artifacts = [
        METRICS_CSV, TRAJECTORY_NPZ, SOURCE_LIST, IMAGE_LIST,
        SUMMARY_FIGURE, PROVENANCE_FIGURE, DECISION_FIGURE,
        REPORT_TXT, VALIDATION,
    ]
    _write_json(MANIFEST, {
        "schema_version": "a4-c1-4-freeze-v1",
        "status": FINAL_STATUS,
        "created_utc": _utc(),
        "restore_policy": (
            "strict_no_recompute; reject any source or artifact hash drift; "
            "explicit rebuild may not change frozen physics or gates"
        ),
        "source_set_short_hash": source_tag,
        "source_inventory_sha256": _sha256(SOURCE_LIST),
        "source_files": inventory["items"],
        "frozen_artifacts": [
            {"path": str(path.resolve()), "sha256": _sha256(path)}
            for path in artifacts
        ],
    })
    print(f"A4-C1-4: OK {FINAL_STATUS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
