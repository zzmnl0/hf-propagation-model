"""Aggregate bounded A3Rnew evidence and decide limited A4 readiness."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "output2" / "a3rnew_p5"
RESULT = OUT / "a3rnew_p5_acceptance.json"
REPORT = ROOT / "word" / "a3rnew_p5_accuracy_acceptance_report.txt"

EVIDENCE = {
    "p0c_preparation": ROOT / "output2/a3r_benchmark/a3r_p0c_reference_preparation.json",
    "p0c_closure": ROOT / "output2/a3r_benchmark/a3r_p0c_closure.json",
    "p2_6l": ROOT / "output2/a3rnew_p2l/p2_6l_reference_acceptance.json",
    "p3": ROOT / "output2/a3rnew_p3/a3rnew_p3_validation.json",
}

FOCUSED_CHECKS = (
    "tests/_verify_AH3D_physics.py",
    "tests/_verify_ah3d_frequency_continuity.py",
    "tests/_verify_ah3d_convergence.py",
    "tests/_verify_p2p_multiscale_reference.py",
    "tests/_verify_p0c_high_low_reference.py",
    "tests/_verify_a3rnew_p2l_reference.py",
    "tests/_verify_a3rnew_p3_validation.py",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _run_check(relative_path: str) -> dict:
    started = time.perf_counter()
    try:
        run = subprocess.run(
            [sys.executable, relative_path], cwd=ROOT, text=True,
            capture_output=True, timeout=600, check=False,
        )
        output = (run.stdout + "\n" + run.stderr).strip().splitlines()
        return {
            "status": "passed" if run.returncode == 0 else "failed",
            "returncode": run.returncode,
            "wall_s": round(time.perf_counter() - started, 6),
            "output_tail": output[-4:],
        }
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout_failed", "returncode": None,
            "wall_s": round(time.perf_counter() - started, 6),
            "output_tail": ["focused check exceeded 600 s"],
        }


def _max_integration_deltas(closure: dict) -> dict:
    landing, height, delay = [], [], []
    for repeat in closure["repeats"]:
        for item in repeat["independent"]:
            research, reference = item["research"], item["reference"]
            landing.append(((research["x_land_km"] - reference["x_land_km"]) ** 2
                            + (research["y_land_km"] - reference["y_land_km"]) ** 2) ** 0.5)
            height.append(abs(research["h_reflect_km"] - reference["h_reflect_km"]))
            delay.append(abs(research["tau_ms"] - reference["tau_ms"])
                         / max(abs(reference["tau_ms"]), 1e-15))
    return {
        "landing_km": max(landing),
        "reflection_height_km": max(height),
        "group_delay_relative": max(delay),
    }


def _evidence_integrity(p0c: dict, p2: dict, p3: dict) -> dict:
    from legacy.models.p2p_mode_audit_3d import _code_hash, _hash

    p3_prerequisites = p3["provenance"]["prerequisite_hashes"]
    p3_code = p3["provenance"]["code_hashes"]
    audit_path = next((ROOT / "output2/a3rnew_p2l/audit").rglob("audit.json"))
    audit = _load(audit_path)
    frozen = audit["frozen"]
    repeat_provider_hashes = {item["provider_hash"] for item in p0c["repeats"]}
    repeat_config_hashes = {item["provider_config_hash"] for item in p0c["repeats"]}
    return {
        "p3_prerequisite_hashes_match": all(
            _sha256(ROOT / path) == expected for path, expected in p3_prerequisites.items()),
        "p3_code_hashes_match": all(
            _sha256(ROOT / path) == expected for path, expected in p3_code.items()),
        "p2_checkpoint_frozen_hash_match": _hash(frozen) == audit["frozen_hash"],
        "p2_checkpoint_code_hash_match": _code_hash() == frozen["code_hash"],
        "provider_hash_chain_match": (
            len(repeat_provider_hashes) == 1
            and p2["provider_hash"] == p2["reference_provider_hash"]
            == frozen["provider_hash"] == p3["coleman_2011"]["frozen_inputs"]["provider_hash"]
            == next(iter(repeat_provider_hashes))
        ),
        "provider_config_hash_chain_match": (
            len(repeat_config_hashes) == 1
            and frozen["config_hash"] == next(iter(repeat_config_hashes))
        ),
        "p2_checkpoint_complete": audit["completed"] is True,
    }


def _write_report(result: dict) -> None:
    metrics = result["metrics"]
    checks = result["focused_checks"]
    risks = "\n".join(f"{i}. {risk}" for i, risk in enumerate(result["residual_risks"], 1))
    check_lines = "\n".join(
        f"{i}. {name}: {record['status']} ({record['wall_s']:.3f} s)"
        for i, (name, record) in enumerate(checks.items(), 1)
    )
    text = f"""A3Rnew-P5new 精度验收与A4有限准入报告
日期：{result['date']}

一、验收语义

P5new仅汇总P0C、P1、P2-6L和P3的冻结证据并执行短时回归，不重复高密度角域、复杂场O/X完整性或仰角精细分层计算。AH Hamiltonian方程、Appleton--Hartree色散、DOP853积分、D=0投影和多尺度中心差分Jacobian均未简化。结论为a4_limited_scope_ready，不是复杂三维场全模式完整性通过。

二、分项状态

1. trajectory_accuracy_status={result['statuses']['trajectory_accuracy_status']}
2. p2p_correctness_status={result['statuses']['p2p_correctness_status']}
3. completeness_status={result['statuses']['completeness_status']}
4. literature_status={result['statuses']['literature_status']}
5. performance_status={result['statuses']['performance_status']}
6. external_validation_status={result['statuses']['external_validation_status']}

三、关键精度证据

1. 最大Hamiltonian残差：{metrics['hamiltonian_residual_max']:.6e}，阈值不大于1e-6。
2. research/reference最大落点差：{metrics['integration_deltas_max']['landing_km']:.6f} km，阈值不大于1 km。
3. 最大反射高度差：{metrics['integration_deltas_max']['reflection_height_km']:.6f} km，阈值不大于0.5 km。
4. 最大群时延相对差：{metrics['integration_deltas_max']['group_delay_relative']:.6e}，阈值不大于1e-3。
5. P2P最大命中误差：{metrics['p2p_landing_error_max_km']:.6f} km，阈值小于5 km。
6. Coleman高/低射线召回率：{metrics['p0c_mode_recall']:.1%}；正式根数={metrics['p0c_formal_mode_count']}，未对其余风险候选作额外正式根判定。

四、轻量回归

{check_lines}

五、A4准入

gate_status={result['gate_status']}。A4仅允许先开展平滑场、可识别单轨迹/O-X分支的PHaRLAP对照；不得解释为复杂场模式完整性验证。

六、风险与未完成项

{risks}

七、理论依据与文献语义

P5new未新增物理参数或方程，沿用Haselgrove Hamiltonian射线理论、Appleton--Hartree磁离子色散、Jones--Stephenson三维射线追踪，以及Coleman（2011）和Nosikov等（2020）的高/低射线识别语义。完整引用统一维护在word/reference-papers.txt；Coleman结果为literature_guided_parameterized_reproduction，Nosikov结果为partial_reproduction_with_assumptions，不宣称严格复现。
"""
    REPORT.write_text(text, encoding="utf-8")


def main() -> int:
    missing = [str(path.relative_to(ROOT)) for path in EVIDENCE.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing frozen evidence: " + ", ".join(missing))

    prep = _load(EVIDENCE["p0c_preparation"])
    p0c = _load(EVIDENCE["p0c_closure"])
    p2 = _load(EVIDENCE["p2_6l"])
    p3 = _load(EVIDENCE["p3"])
    focused = {path: _run_check(path) for path in FOCUSED_CHECKS}
    integrity = _evidence_integrity(p0c, p2, p3)
    integration = _max_integration_deltas(p0c)
    thresholds = p3["thresholds"]
    modes = p2["p0c_modes"]
    metrics = {
        "hamiltonian_residual_max": p3["coleman_2011"]["maximum_errors"]["hamiltonian_residual"],
        "integration_deltas_max": integration,
        "p2p_landing_error_max_km": max(mode["p2p_landing_error_km"] for mode in modes),
        "p0c_mode_recall": sum(item["matched"] for item in p2["p0c_matches"])
        / len(prep["independent_roots_deg"]),
        "p0c_formal_mode_count": len(modes),
        "p0c_unresolved_candidate_count": len(p2["p0c_unresolved_not_evaluated"]),
    }
    acceptance_checks = {
        "evidence_integrity": all(integrity.values()),
        "focused_checks": all(item["status"] == "passed" for item in focused.values()),
        "p0c_passed": p0c["acceptance"]["p0_overall_status"] == "passed",
        "p1_multiscale_in_production": all(
            mode["p2p_jacobian_effective"] == "central_reference_multiscale" for mode in modes),
        "p2_limited_passed": p2["status"] == "reference_case_limited_pass",
        "p3_passed": p3["status"] == "passed" and all(p3["checks"].values()),
        "trajectory_thresholds": (
            metrics["hamiltonian_residual_max"] <= thresholds["hamiltonian_residual_le"]
            and integration["landing_km"] <= thresholds["integration_landing_delta_km_le"]
            and integration["reflection_height_km"] <= thresholds["integration_height_delta_km_le"]
            and integration["group_delay_relative"] <= thresholds["integration_delay_relative_le"]
        ),
        "p2p_thresholds": metrics["p2p_landing_error_max_km"] < thresholds["landing_error_km_lt"],
        "bounded_completeness_semantics": (
            p2["not_run"]["high_density_frozen_reference_recall"] == "deferred_as_P2-6H"
            and p2["not_run"]["complex_scene_ox_p2p_completeness"] == "not_run"
        ),
    }
    ready = all(acceptance_checks.values())
    result = {
        "schema_version": "a3rnew-p5-v1",
        "stage": "A3Rnew-P5new",
        "date": date.today().isoformat(),
        "status": "a4_limited_scope_ready" if ready else "a4_not_ready",
        "gate_status": "a4_limited_scope_ready" if ready else "a4_not_ready",
        "semantics": "accuracy evidence aggregation; not complex-scene full-mode completeness",
        "equation_policy": "production equations and numerical tolerances unchanged",
        "statuses": {
            "trajectory_accuracy_status": "internal_reference_pass" if ready else "internal_reference_failed",
            "p2p_correctness_status": "reference_case_pass" if ready else "reference_case_failed",
            "completeness_status": "reference_case_limited_pass" if ready else "reference_case_not_passed",
            "literature_status": "literature_guided_and_partial",
            "performance_status": "recorded_not_gated",
            "external_validation_status": "not_run_until_A4",
        },
        "acceptance_checks": acceptance_checks,
        "evidence_integrity": integrity,
        "evidence_hashes": {
            name: _sha256(path) for name, path in EVIDENCE.items()
        },
        "metrics": metrics,
        "thresholds": thresholds,
        "focused_checks": focused,
        "performance_evidence": {
            "source": "stored P0C/P2-6L diagnostics plus current focused-check wall time",
            "p0c_repeat_count": len(p0c["repeats"]),
            "p2_6l_profile_wall_s": p2["p0c_diagnostics"]["solver"]["profile_wall_s"],
            "focused_checks_wall_s": sum(item["wall_s"] for item in focused.values()),
            "speedup_gate": "not_applicable",
        },
        "not_run": {
            **p2["not_run"],
            "p2_6h": "deferred_full_completeness_audit",
            "fine_elevation_layering": "not_run_by_design",
            "performance_benchmark_rerun": "not_run_by_design",
            "pharlap_external_validation": "not_run_until_A4",
        },
        "residual_risks": [
            "P2-6H复杂场高密度完整性审计继续延期",
            "baseline/TID的O/X P2P完整性、连续两级稳定性及continuation等价性尚未验证",
            "变分Jacobian精度闭环继续延期，且不用于生产验收",
            "near-cutoff不可微、隧穿和模式转换超出射线理论有效范围",
            "AH3D仍是有限差分Hamiltonian射线追踪，不是全波求解器",
            "local_2p5d PE不是完整3D PE/全波，尚未闭合路径间相干耦合",
            "绝对功率尚未使用实测数据校准",
            "P5new内部证据不能替代A4 PHaRLAP外部验证",
        ],
        "a4_scope": "smooth field and identifiable single trajectory/O-X branch PHaRLAP comparison only",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _write_report(result)
    print("A3Rnew-P5new: {}".format(result["gate_status"]))
    return 0 if ready else 1


if __name__ == "__main__":
    raise SystemExit(main())
