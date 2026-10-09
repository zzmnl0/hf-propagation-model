"""Assemble P3.1 evidence without promoting unresolved physical mechanisms."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_pass(physical, strict, reference):
    return bool(physical and (strict or reference))


def build(run, parent, output):
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    frozen = read(run / "baseline_freeze.json")
    delivery = read(parent / "delivery_sha256.json")
    history_checks = [sha(parent / "delivery_sha256.json") == frozen["parent_manifest_sha256"]]
    for row in delivery["records"]:
        history_checks.append(sha(parent / row["path"]) == row["sha256"])
    for row in delivery["external_tracked_sources_and_report"]:
        history_checks.append(sha(Path(row["path"])) == row["sha256"])
    unchanged = all(sha(ROOT / row["path"]) == row["sha256"] for row in frozen["records"])
    write(output / "parent_and_source_verification.json", {
        "parent_entries_checked": len(history_checks) - 1,
        "parent_unchanged": all(history_checks),
        "preexisting_source_files_unchanged": unchanged,
        "source_files_checked": len(frozen["records"]),
    })
    index = read(run / "reference_certification_index.json")
    ref = read(run / index["analysis_directory"] / "summary.json")
    refs = {row["case_id"]: row for row in ref["cases"]}
    independent = read(run / index["independent_verification"])
    independent_passed = independent["all_16_reproduced"] is True
    independently_verified = {row["case_id"]: row["passed"] is True for row in independent["cases"]}
    original = read(parent / "analysis/p3_refined_v1/summary.json")
    rows = []
    for old in original["cases"]:
        cert = refs.get(old["case_id"])
        strict = old["strict_passed"] is True
        reference = bool(cert and cert["reference_certified_passed"] and independently_verified.get(old["case_id"], False))
        if cert:
            physical = all(cert["physical_checks"].values()) and cert["checks"]["initial_state"] and cert["checks"]["scope"]
            assert cert["strict_convergence_passed"] == strict
        else:
            physical = bool(strict and old["topology"]["eligible"] and old["topology"]["coverage_complete"])
        passed = case_pass(physical, strict, reference)
        rows.append({"case_id": old["case_id"], "role": old["role"], "scope": old["scope"],
                     "strict_convergence_passed": strict, "reference_certified_passed": reference,
                     "common_physics_passed": bool(physical), "passed": passed,
                     "acceptance_route": "strict_parent_reused" if strict else "independent_reference" if reference else "none",
                     "execution": "parent_raw_reused" if passed else "blocked_actual_entrance_medium",
                     "reason": "" if passed else "Negative centered-Hermite support; PH propagation rule not uniquely identified."})
    for angle in (27, 33):
        rows.append({"case_id": f"new_holdout_f11_e{angle}", "role": "new_holdout", "scope": "fullpath",
                     "strict_convergence_passed": False, "reference_certified_passed": False,
                     "common_physics_passed": False, "passed": False, "acceptance_route": "none",
                     "execution": "not_executed", "reason": "Actual medium and reference topology registration pending."})
    assert len(rows) == 26 and len(refs) == 16
    boundary = read(run / "boundary_independent_check/verification.json")
    full = read(run / "reference_fullpath_v5/controls.json")
    piecewise = read(run / "reference_fullpath_piecewise/verification.json")
    defaults = read(run / "after/equivalence.json")["passed"] and read(run / "after/power_equivalence.json")["passed"]
    engineering = read(run / "r0_tests/v1/summary.json")["full_suite_passed"]
    boundary_identified = boundary["propagation_rule_identified"] is True
    full_formal = full["formal_smooth_fullpath_verified"] is True
    flags = {
        "production_implementation_complete": False,
        "boundary_output_rule_identified": bool(boundary["passed"]),
        "boundary_convention_identified": boundary_identified,
        "fullpath_reference_implementation_verified": bool(full["passed"] and piecewise["passed"]),
        "fullpath_independent_reconstruction_verified": full_formal,
        "strict_convergence_passed": all(r["strict_convergence_passed"] for r in rows),
        "reference_certified_passed": all(r["reference_certified_passed"] for r in rows),
        "r3_reference_certified_passed": bool(ref["all_16_certified"] and independent_passed),
        "default_regression_passed": bool(defaults),
        "engineering_tests_passed": bool(engineering),
        "entrance_consistency_passed": False,
        "post_repair_fd_grid_regression_passed": False,
        "holdout_passed": all(r["passed"] for r in rows if r["role"] == "new_holdout"),
        "source_verification_passed": bool(all(history_checks) and unchanged and independent_passed),
        "p31_scope_passed": False, "recommend_p4": False,
    }
    summary = {"schema": "p31-delivery-v1", "flags": flags, "cases": rows,
               "case_count": len(rows), "passed_count": sum(r["passed"] for r in rows),
               "strict_passed_count": sum(r["strict_convergence_passed"] for r in rows),
               "reference_passed_count": sum(r["reference_certified_passed"] for r in rows),
               "production_modified": False, "reference_index": index,
               "stage_status": "PARTIAL_BLOCKED_AT_ATTRIBUTION_REVIEW",
               "scope": "6371 km sphere, radial Ne, B=0, collisionless, single-hop forward",
               "overall_3d_magnetoionic_equivalence_passed": False}
    summary["r3_maximum_cross_error"] = {key: max(row["actual"][key] for row in refs.values())
                                          for key in ("rms_m", "max_m", "landing_m", "group_rel")}
    write(output / "summary.json", summary)
    with (output / "case_matrix.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    write(output / "gate_self_checks.json", {
        "wrong_physics_cannot_pass": not case_pass(False, True, True),
        "reference_does_not_require_strict": case_pass(True, False, True),
        "no_route_cannot_pass": not case_pass(True, False, False),
    })
    causes = [
        {"cause": "PH exported density follows clipped centered-slope Hermite", "status": "confirmed_output_only", "evidence": "72 boundary calls; independent native-node reproduction", "repair": "None: exported density alone does not identify propagation."},
        {"cause": "PH propagation uses a pure clipped Hamiltonian", "status": "unresolved", "evidence": "30 native state/output mismatches; 25 return to unit momentum at next node", "repair": "No production clipping or reset compensation."},
        {"cause": "Raw polynomial without state reset explains every boundary state", "status": "excluded_as_complete_explanation", "evidence": "Negative-support interior nodes repeatedly return to unit momentum", "repair": "Do not implement a raw-density-only compatibility rule."},
        {"cause": "R3 strict failure proves large propagation error", "status": "excluded_in_certified_scope" if flags["r3_reference_certified_passed"] else "unresolved", "evidence": "Independent 80/120/160 dps reference and separate quadrature scan", "repair": "Separate reference route; retain strict failures."},
        {"cause": "All R3 nonmonotonicity is floating point roundoff", "status": "unresolved", "evidence": "Accuracy certification does not uniquely explain nonmonotonicity", "repair": "No integrator change based on this assertion."},
        {"cause": "Fullpath independent reference unavailable", "status": "implementation_gap_closed_on_legal_controls", "evidence": "Turning, ground, group, early termination and piecewise controls", "repair": "Connect formal smooth cases only after actual medium identified."},
    ]
    write(output / "cause_matrix.json", causes)
    text = f"""P3.1 剩余误差闭合交付报告

状态：部分完成，停在入口传播规则归因评审节点；recommend_p4=false。
正式矩阵：{summary['passed_count']}/26。严格通道 {summary['strict_passed_count']} 例；独立参考通道 {summary['reference_passed_count']} 例。
原 P3 的 2/24 严格通过记录保持不变。新增通道不表示严格单调收敛，也不证明全部非单调变化都来自舍入。

一、入口归因
训练48次、8/15 MHz留出24次PH运行均保留原smooth_s5输入。导出Ne符合截零后的中心斜率Hermite，最大节点差1.4841556549072266e-5 m^-3。
独立复核发现30个状态/输出不一致节点，其中12个来自留出；25个在下一原生节点恢复单位波矢模长。
输出规则的证据充分，传播规则仍不唯一：可能涉及零密度分支的状态投影、跨界积分或输出状态时序。
纯未截零且从不重置的模型不能解释全部状态；纯截零Hamiltonian也尚未得到唯一确认。
精确20 km端点已作位置和原生群程PCHIP诊断，未采用跨重置的状态导数。
第三档clip候选最大端点差2.84652e-5 m，raw候选6.91465e-5 m。如此小的轨迹差不能替代机制辨识。
PH端点重建缺少独立误差包络，formal_20km_endpoint_acceptance=false。
0.5/0.25/0.125 km方向差分已记录，但依赖输出重建，未作为PH内部梯度的独立收敛证明，也不称为完整三维梯度。
最小复现：boundary/minimal_repro/run_minimal_repro.m；保留逐点状态及输入MAT供PH内部规则核实。

二、R3独立认证
精确转换冻结双精度输入与系数，使用80/120/160位精度和独立求积阶数扫描；不重新拟合斜率或常数。
分别记录AH—参考、PH—参考和AH—PH，查询原生节点、100 km实际右支撑及正式比较点。
传播包络包含位置误差和h/4乘方向导数误差；纯重建、统计采样和参考数值误差分别记账。
参考稳定性是数值估计，不是严格区间证明。详见reference_certification_index.json指向的结果及独立复算。
首次全量独立复算有一例预算复现失败：独立验收器把100 km以后的支撑区间曲率也纳入正式区间界限。已将共同曲率计算限定到[0,100]，保留右侧支撑用于同一插值；新增域外高曲率反例，未修改物理阈值或v4原预算。失败记录保留在reference_independent_final_v1。
16例跨程序最大指标：{json.dumps(summary['r3_maximum_cross_error'], ensure_ascii=False)}。

三、全程参考
合法仿射和非线性分段三次介质的独立径向积分与Cartesian控制通过。
包括上行、转向、下行、自然落点、群程和提前终止反例；参考实现验证为真。
原smooth_s5全程实际介质仍未明确，因此formal全程重建标志保持false，不能以控制例替代原四角。

四、生产与留出
本轮没有生产补丁。原严格拒绝负密度、默认PCHIP路径及已有显式兼容配置保持原状。
六项生产工程测试通过；默认四例26项逐值相等检查和4项功率检查通过。
FD/网格没有新增生产设置，保留P3原记录；入口修复后的受影响回归仍待执行，不标为新通过。
11 MHz、27/33度历史暴露审计未发现既往运行，但实际介质与拓扑未冻结，尚未执行该留出；预登记备选也未运行。

五、下一步
使用现有最小MATLAB复现核实PH零密度分支、状态投影及输出采样时序；须同时解释上/下行、幅度和留出状态。
若仍无法区分内部机制，继续保留未解决，不增加经验密度缩放或按角度调参。
唯一归因后提交共同根因补丁或显式兼容设计，再运行原四角、12 MHz指定回归与预登记新留出。
只有26例规定通道、全程自然终点/群程、FD/网格及来源复核全部满足，才重新评估P4准入。

限制：本轮只覆盖6371 km球面、径向Ne、B=0、无碰撞、单跳前向。三维横向、O/X、IGRF和P2P均未执行。
"""
    (output / "report.txt").write_text(text, encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = sorted(refs)
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    for ax, key, label in zip(axes, ("rms_m", "group_rel"), ("RMS position (m)", "Group relative difference")):
        ax.semilogy(range(len(names)), [max(refs[n]["actual"][key], 1e-20) for n in names], "o-", label="AH - PH")
        ax.semilogy(range(len(names)), [refs[n]["budget"][key] for n in names], ".-", label="Independent numerical budget")
        ax.set_ylabel(label); ax.grid(alpha=.3); ax.legend()
        ax.set_xticks(range(len(names)), names, rotation=65, ha="right", fontsize=7)
    fig.savefig(output / "r3_reference_budgets.png", dpi=160); plt.close(fig)
    print(json.dumps({"passed": summary["passed_count"], "total": 26, "recommend_p4": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.run.resolve(), args.parent.resolve(), args.output.resolve())
