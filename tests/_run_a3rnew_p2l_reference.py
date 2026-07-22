"""P2-6L lightweight full-domain reference-case acceptance.

This is deliberately a limited reference gate.  It does not validate mode
completeness for baseline or TID providers; those remain recorded as P2-6H.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(TESTS))

import config as cfg
from _benchmark_p2p_ah3d import CASES, _build_case, _mode_record
from _verify_AH3D import run as verify_ah3d
from _verify_a3rnew_p2_audit import main as verify_audit
from _verify_p0c_high_low_reference import main as verify_p0c_history
from _verify_p2p_multiscale_reference import run as verify_multiscale
from models import point_to_point_ah_3d as p2p
from legacy.models.p2p_mode_audit_3d import adaptive_landing_audit
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D
from models.scenario_provider_factory import ScenarioProviderFactory


OUT = ROOT / "output2" / "a3rnew_p2l"
P0C = ROOT / "output2" / "a3r_benchmark" / "a3r_p0c_reference_preparation.json"


def _jsonable(value):
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(_jsonable(value), indent=2, ensure_ascii=True),
                         encoding="utf-8")
    temporary.replace(path)


def _formal_ok(mode: dict) -> bool:
    reference = mode.get("p2p_jacobian_reference", {})
    return bool(
        mode.get("branch_status") == "ground"
        and mode.get("convergence_status") == "passed"
        and mode.get("p2p_jacobian_effective") == "central_reference_multiscale"
        and bool(reference.get("stable"))
        and float(mode.get("p2p_landing_error_km", np.inf)) < 5.0
        and mode.get("convergence_reference_profile") == "reference"
    )


def _match_roots(modes: list[dict], roots: list[dict]) -> tuple[list[dict], list[dict]]:
    available = set(range(len(modes)))
    matches = []
    for root in roots:
        ranked = sorted(available, key=lambda index: abs(
            float(modes[index]["beta_deg_tx"]) - float(root["beta_deg"])))
        index = ranked[0] if ranked else None
        mode = modes[index] if index is not None else None
        beta_error = (abs(float(mode["beta_deg_tx"]) - float(root["beta_deg"]))
                      if mode is not None else np.inf)
        height_error = (abs(float(mode["h_reflect_km"]) - float(root["h_reflect_km"]))
                        if mode is not None else np.inf)
        matched = bool(mode is not None and beta_error <= 0.05 and height_error <= 1.0)
        if matched:
            available.remove(index)
        matches.append({"ray_label": root["ray_label"], "matched": matched,
                        "beta_error_deg": beta_error, "height_error_km": height_error,
                        "mode": _mode_record(mode) if mode is not None else None})
    return matches, [_mode_record(modes[index]) for index in sorted(available)]


def _stats() -> dict:
    return {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {}, "profile_accepted_steps": {},
        "reference_ray_calls": 0, "reference_nfev": 0,
        "reference_accepted_steps": 0, "jacobian_rejected_local_roots": {},
    }


def _reference_cell(seed: dict, root: dict) -> bool:
    (beta0, beta1), (az0, az1) = seed["angle_bounds"]
    beta = float(root["beta_deg"])
    return beta0 <= beta <= beta1 and az0 - 1.0e-10 <= 0.0 <= az1 + 1.0e-10


def _select_reference_seeds(seeds: list[dict], roots: list[dict]) -> list[dict]:
    selected = []
    for root in roots:
        candidates = [seed for seed in seeds if _reference_cell(seed, root)]
        if not candidates:
            raise RuntimeError("P0C reference root is outside all audited risk cells")
        selected.append(max(candidates, key=lambda seed: seed["audit_level"]))
    return selected


def _reference_seed(seed: dict, root: dict) -> dict:
    (beta0, beta1), _ = seed["angle_bounds"]
    beta = float(root["beta_deg"])
    return {
        **seed,
        "beta_deg": beta,
        "az_deg": 0.0,
        "provenance": "independent_p0c_reference",
        "angle_bounds": ((beta0, beta1), (-1.0, 1.0)),
    }


def run() -> int:
    prep = json.loads(P0C.read_text(encoding="utf-8"))
    if prep.get("status") != "p0c_1_passed":
        raise RuntimeError("P0C reference preparation is not passed")
    spec = next(item for item in CASES if item["name"] == "coleman_chapman_high_low")
    provider, n3d, setup_s = _build_case(ScenarioProviderFactory(cfg), spec)
    roots = prep["root_rays"]
    params = {
        **cfg.P2P_3D,
        "audit_mode": "adaptive_2d",
        "audit_beta_bounds_deg": (2.0, 88.0),
        "audit_az_bounds_deg": (-180.0, 180.0),
        "audit_beta_step_deg": 15.0,
        "audit_az_step_deg": 60.0,
        "audit_max_refinement": 1,
        "audit_target_near_km": 150.0,
        "audit_fold_condition": 1.0e4,
        "audit_checkpoint_dir": str(OUT / "audit"),
        "audit_resume": True,
        "audit_progress": True,
        "candidate_nearest": 1,
        "jacobian": "central_reference_multiscale",
    }
    target = (float(prep["selected_target_distance_km"]), 0.0)
    n_ah = RefractiveIndexAH3D(
        n3d, float(spec["freq_MHz"]), "O",
        geomag={**cfg.GEOMAG, "fH_MHz": 0.0}, integration_profile="smoke")
    trace_cache, stats = {}, _stats()
    seeds, audit = adaptive_landing_audit(
        cfg.TX_POS, target, n_ah, float(spec["freq_MHz"]), params, 0.0,
        lambda beta, az: p2p._trace_cached(
            cfg.TX_POS, n_ah, float(spec["freq_MHz"]), beta, az, 0.0, False,
            "smoke", trace_cache, stats))
    selected_cells = _select_reference_seeds(seeds, roots)
    selected = [_reference_seed(cell, root)
                for cell, root in zip(selected_cells, roots)]
    modes = []
    selected_records = []
    for root, seed in zip(roots, selected):
        local = p2p._shoot_ah_3d(
            cfg.TX_POS, target, n_ah, float(spec["freq_MHz"]),
            seed["beta_deg"], seed["az_deg"], max_iter=int(params["max_iter"]),
            tol_km=float(params["tol_km"]), d_beta=float(params["d_beta"]),
            d_az=float(params["d_az"]), jacobian=params["jacobian"], t=0.0,
            jacobian_rtol=float(params["jacobian_rtol"]),
            jacobian_scales=tuple(params["jacobian_scales"]),
            jacobian_min_scale=float(params["jacobian_min_scale"]),
            jacobian_scale_rtol=float(params["jacobian_scale_rtol"]),
            jacobian_condition_max=float(params["jacobian_condition_max"]),
            trace_cache=trace_cache, trace_stats=stats, defer_formal=True,
            seed_info=seed, angle_bounds=seed["angle_bounds"])
        mode = (p2p._formalize_ah_3d(
            cfg.TX_POS, target, n_ah, float(spec["freq_MHz"]), local,
            float(params["tol_km"]), 0.0, "research", trace_cache, stats)
            if local is not None else None)
        selected_records.append({"ray_label": root["ray_label"], "seed": seed,
                                 "formal_passed": mode is not None})
        if mode is not None:
            modes.append(mode)
    diagnostics = {"audit": audit, "solver": stats}
    matches, unexpected = _match_roots(modes, roots)
    reference_checks = {
        "p0c_history": verify_p0c_history() == 0,
        "p0c_lightweight_recall": all(item["matched"] for item in matches),
        "p0c_unselected_risk_candidates_recorded": bool(seeds),
        "p0c_all_formal_modes_strict": bool(modes) and all(_formal_ok(mode) for mode in modes),
        "ah3d_ox_direct_reference": verify_ah3d() == 0,
        "audit_checkpoint_periodic": verify_audit() == 0,
        "multiscale_jacobian": verify_multiscale() == 0,
    }
    passed = all(reference_checks.values())
    report = {
        "schema": "a3rnew-p2-6l-v1",
        "status": "reference_case_limited_pass" if passed else "reference_case_limited_fail",
        "semantics": "P2-6L reference-only acceptance; not complex-scene mode completeness",
        "configuration": {
            "label": "lightweight_full_domain_preacceptance",
            "wave_modes": ["O"],
            "p0c_geomag": "B=0",
            "full_beta_domain_deg": [2.0, 88.0],
            "full_az_domain_deg": [-180.0, 180.0],
            "beta_step_deg": 15.0,
            "az_step_deg": 60.0,
            "max_refinement": 1,
            "jacobian": "central_reference_multiscale",
            "formal_profiles": ["research", "reference"],
            "candidate_nearest": 1,
        },
        "provider_hash": provider.get_metadata().get("provider_hash"),
        "reference_provider_hash": prep["provider_metadata"].get("provider_hash"),
        "setup_wall_s": setup_s,
        "reference_checks": reference_checks,
        "p0c_matches": matches,
        "p0c_unexpected_evaluated_formal_roots": unexpected,
        "p0c_reference_selected_cells": selected_records,
        "p0c_unresolved_not_evaluated": [
            seed for seed in seeds if seed not in selected_cells
        ],
        "p0c_modes": [_mode_record(mode) for mode in modes],
        "p0c_diagnostics": diagnostics,
        "not_run": {
            "baseline_mode_completeness": "not_verified; prior zero-mode rejection chain remains unresolved",
            "tid_mode_completeness": "not_run",
            "complex_scene_ox_p2p_completeness": "not_run",
            "two_level_mode_set_stability": "not_run",
            "high_density_frozen_reference_recall": "deferred_as_P2-6H",
            "continuation_equivalence": "not_run",
        },
        "deferred_gate": "P2-6H deferred_full_completeness_audit",
    }
    _write(OUT / "p2_6l_reference_acceptance.json", report)
    print("P2-6L status={} modes={} output={}".format(
        report["status"], len(modes), OUT / "p2_6l_reference_acceptance.json"))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(run())
