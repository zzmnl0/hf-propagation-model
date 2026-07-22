"""A3R-P0 reproducible single-process AH3D P2P benchmark.

The frozen set is a moderate-density numerical reference, not physical truth.
Run from ``code`` with the pytorch_cpu environment::

    python tests/_benchmark_p2p_ah3d.py
"""

from __future__ import annotations

import argparse
import cProfile
import hashlib
import json
import os
import platform
import random
import statistics
import sys
import time
import pstats
from pathlib import Path

# Set numerical-library threads before importing NumPy/SciPy.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import scipy
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
import legacy.models.point_to_point_ah_3d as p2p_ah
from models.dispersion_ah import appleton_hartree_roots
from legacy.models.ne_provider import BaseNeProvider3D
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D
from models.scenario_provider_factory import ScenarioProviderFactory


SCHEMA_VERSION = "a3r-p0c-v1"
BENCHMARK_FILENAME = "a3r_p0r_benchmark.json"
REFERENCE_SEED_COUNT = 8  # Reduced from a dense scan to bound formal runtime.
CASES = (
    {"name": "vacuum_uniform_degenerate", "scenario": "vacuum", "range_km": 250.0,
     "freq_MHz": 30.0, "expected_semantics": "no_reflecting_mode_allowed"},
    {"name": "smooth_layer_single_mode", "scenario": "baseline", "range_km": 600.0,
     "freq_MHz": 10.0, "candidate_nearest": 1,
     "expected_semantics": "one_or_more_converged_modes"},
    {"name": "coleman_chapman_high_low", "scenario": "coleman_chapman",
     "range_km": 1100.0, "freq_MHz": 10.0, "candidate_nearest": 3,
     "expected_semantics": "literature_guided_same_branch_high_low_reference"},
    {"name": "tid_3d", "scenario": "tid", "range_km": cfg.RX_RANGE,
     "freq_MHz": cfg.FREQ_MHZ, "candidate_nearest": 1,
     "expected_semantics": "parameterized_3d_tid"},
)


class _VacuumProvider(BaseNeProvider3D):
    """Exact zero-density 3-D provider for the non-reflecting limit."""

    name = "vacuum_3d"

    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        return 0.0

    def get_Ne_grid(self, x_km, z_km, t: float = 0.0, y_km: float = 0.0):
        return np.zeros((len(x_km), len(z_km)), dtype=float)

    def get_metadata(self) -> dict:
        return {
            "name": self.name, "dimensions": "3d", "units": "m^-3",
            "config_hash": _stable_hash({"model": self.name}),
            "provider_hash": _stable_hash({"model": self.name, "Ne": 0.0}),
        }


class _ColemanChapmanProvider(BaseNeProvider3D):
    """Coleman (2011) uniform Chapman F2 benchmark field."""

    name = "coleman_2011_chapman_f2_3d"
    foF2_MHz = 7.0
    hmF2_km = 350.0
    ymF2_km = 100.0
    horizontal_gradient_MHz_per_deg = 0.0

    @property
    def peak_Ne_m3(self) -> float:
        return (self.foF2_MHz * 1.0e6 / 8.98) ** 2

    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        z_km = float(np.asarray(pos_km, dtype=float).ravel()[-1])
        u = (z_km - self.hmF2_km) / self.ymF2_km
        return float(self.peak_Ne_m3 * np.exp(0.5 * (1.0 - u - np.exp(-u))))

    def get_Ne_grid(self, x_km, z_km, t: float = 0.0, y_km: float = 0.0):
        z = np.asarray(z_km, dtype=float)
        u = (z - self.hmF2_km) / self.ymF2_km
        profile = self.peak_Ne_m3 * np.exp(0.5 * (1.0 - u - np.exp(-u)))
        return np.repeat(profile[None, :], len(x_km), axis=0)

    def get_metadata(self) -> dict:
        frozen = {
            "model": self.name, "foF2_MHz": self.foF2_MHz,
            "hmF2_km": self.hmF2_km, "ymF2_km": self.ymF2_km,
            "horizontal_gradient_MHz_per_deg": self.horizontal_gradient_MHz_per_deg,
            "collision_frequency_Hz": 0.0, "geomagnetic_field": "B=0",
            "density_formula": "Ne=Nmax*exp(0.5*(1-u-exp(-u)))",
            "u": "(z-hmF2)/ymF2", "source": "Coleman (2011), RS5016",
        }
        return {
            **frozen, "dimensions": "3d_uniform_xy", "units": "m^-3",
            "peak_Ne_m3": self.peak_Ne_m3,
            "config_hash": _stable_hash(frozen),
            "provider_hash": _stable_hash({**frozen, "peak_Ne_m3": self.peak_Ne_m3}),
        }


def _stable_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _atomic_write(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True),
                         encoding="utf-8")
    temporary.replace(path)


def _machine_metadata() -> dict:
    return {
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "processes": 1,
        "threads_per_process": 1,
        "thread_environment": {name: os.environ.get(name) for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS")},
        "timer": "time.perf_counter_ns",
        "memory_semantics": "not_measured_during_timing_to_avoid_tracemalloc_distortion",
    }


def _mode_record(mode: dict) -> dict:
    keys = (
        "wave_mode", "branch_id", "beta_deg_tx", "az_deg_tx", "x_land_km",
        "y_land_km", "h_reflect_km", "geom_path_km", "group_path_km",
        "tau_ms", "p2p_landing_error_km", "convergence_status",
        "branch_status", "p2p_jacobian_effective", "p2p_jacobian_fallback_count",
    )
    out = {}
    for key in keys:
        value = mode.get(key)
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, float) and not np.isfinite(value):
            value = None
        out[key] = value
    return out


def _mode_signature(modes: list[dict]) -> list[tuple]:
    return sorted((
        str(m.get("wave_mode")),
        round(float(m.get("beta_deg_tx", 0.0)), 6),
        round(float(m.get("az_deg_tx", 0.0)), 6),
        round(float(m.get("tau_ms", 0.0)), 6),
    ) for m in modes)


def _build_case(factory, spec: dict):
    started = time.perf_counter_ns()
    provider = (_VacuumProvider() if spec["scenario"] == "vacuum" else
                _ColemanChapmanProvider() if spec["scenario"] == "coleman_chapman" else
                factory.create(spec["scenario"], seed=0, t_s=0.0, validation=True))
    ne_slice = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z, t=0.0, y_km=0.0)
    n3d = RefractiveIndex3D(
        ne_slice, cfg.BG_X, cfg.BG_Z, spec["freq_MHz"], earth_flat=True,
        ne_provider=provider)
    elapsed = (time.perf_counter_ns() - started) * 1e-9
    return provider, n3d, elapsed


def _run_once(factory, spec: dict, n_init: int, wave_mode: str = "both") -> dict:
    provider, n3d, setup_s = _build_case(factory, spec)
    params = {
        **cfg.P2P_3D,
        "seed_source": "internal_scan",
        "seed_betas_deg": [],
        "n_init": int(n_init),
        "jacobian": "central",
        "candidate_nearest": int(spec.get("candidate_nearest", 1)),
    }
    rx = (float(spec["range_km"]), 0.0)
    counters = {"ray_calls": 0, "ray_rhs_calls": 0, "accepted_steps": 0}
    if spec["scenario"] == "vacuum":
        cpu0 = time.process_time_ns()
        wall0 = time.perf_counter_ns()
        dip = np.radians(float(cfg.GEOMAG["dip_deg"]))
        decl = np.radians(float(cfg.GEOMAG["decl_deg"]))
        b_hat = np.array([np.cos(dip) * np.cos(decl),
                          np.cos(dip) * np.sin(decl), -np.sin(dip)])
        roots = appleton_hartree_roots(
            2.0 * np.pi * float(spec["freq_MHz"]) * 1.0e6,
            0.0, np.array([1.0, 0.0, 0.0]), b_hat,
            float(cfg.GEOMAG["fH_MHz"]))
        wall_s = (time.perf_counter_ns() - wall0) * 1e-9
        cpu_s = (time.process_time_ns() - cpu0) * 1e-9
        metadata = provider.get_metadata()
        return {
            "benchmark_operation": "analytic_O_X_vacuum_dispersion_limit",
            "setup_wall_s": setup_s,
            "algorithmic_wall_s": wall_s,
            "algorithmic_cpu_s": cpu_s,
            "peak_python_allocated_bytes": None,
            "provider_hash": metadata.get("provider_hash"),
            "provider_config_hash": metadata.get("config_hash"),
            "mode_count": 0,
            "modes": [],
            "mode_signature": [],
            "fallback_count": 0,
            "solver_nfev": 0,
            "solver_njev": 0,
            "candidate_audit_rays": 0,
            "vacuum_n2": {branch: float(np.real(value))
                           for branch, value in roots.items()},
            "vacuum_n2_max_abs_error": float(max(
                abs(np.real(value) - 1.0) for value in roots.values())),
            **counters,
        }
    original_trace = p2p_ah.trace_single_ray_ah_3d

    def counted_trace(*args, **kwargs):
        call_no = counters["ray_calls"] + 1
        branch = getattr(args[3], "wave_mode", "?") if len(args) > 3 else "?"
        profile = getattr(args[3], "integration_profile", "?") if len(args) > 3 else "?"
        if call_no == 1 or call_no % 10 == 0 or profile != "smoke":
            print("  ray_call={} branch={} profile={} beta={:.4f} az={:.4f} sensitivity={}".format(
                call_no, branch, profile, float(args[1]), float(args[2]),
                bool(kwargs.get("with_sensitivity", False))), flush=True)
        ray = original_trace(*args, **kwargs)
        counters["ray_calls"] += 1
        counters["ray_rhs_calls"] += int(ray.get("nfev", 0))
        counters["accepted_steps"] += int(ray.get("accepted_steps", 0))
        return ray

    cpu0 = time.process_time_ns()
    wall0 = time.perf_counter_ns()
    p2p_ah.trace_single_ray_ah_3d = counted_trace
    p2p_diagnostics = {}
    try:
        modes = p2p_ah.find_all_rays_ah_3d(
            cfg.TX_POS, rx, n3d, float(spec["freq_MHz"]), params,
            t=0.0, wave_mode=wave_mode, geomag=cfg.GEOMAG,
            rays_fan=None, continuation_modes=None,
            diagnostics=p2p_diagnostics)
    finally:
        p2p_ah.trace_single_ray_ah_3d = original_trace
    wall_s = (time.perf_counter_ns() - wall0) * 1e-9
    cpu_s = (time.process_time_ns() - cpu0) * 1e-9
    metadata = provider.get_metadata()
    return {
        "benchmark_operation": "provider_ready_candidate_discovery_plus_P2P",
        "setup_wall_s": setup_s,
        "algorithmic_wall_s": wall_s,
        "algorithmic_cpu_s": cpu_s,
        "peak_python_allocated_bytes": None,
        "provider_hash": metadata.get("provider_hash"),
        "provider_config_hash": metadata.get("config_hash"),
        "scenario": spec["scenario"],
        "provider_components": metadata.get("components", []),
        "mode_count": len(modes),
        "modes": [_mode_record(mode) for mode in modes],
        "mode_signature": _mode_signature(modes),
        "fallback_count": int(sum(
            m.get("p2p_jacobian_fallback_count", 0) for m in modes)),
        "solver_nfev": int(sum(m.get("p2p_nfev", 0) for m in modes)),
        "solver_njev": int(sum(m.get("p2p_njev", 0) for m in modes)),
        "candidate_audit_rays": int(sum(
            m.get("p2p_candidate_audit_rays", 0) for m in modes)),
        **counters,
        "actual_ray_calls_including_reference": int(
            p2p_diagnostics.get("cache_misses", 0)
            + p2p_diagnostics.get("reference_ray_calls", 0)),
        "p0r_diagnostics": p2p_diagnostics,
    }


def _summarize_runs(runs: list[dict]) -> dict:
    signatures = [[tuple(item) for item in run["mode_signature"]] for run in runs]
    fields = ("setup_wall_s", "algorithmic_wall_s", "algorithmic_cpu_s")
    timing = {}
    for field in fields:
        values = [float(run[field]) for run in runs]
        timing[field] = {
            "median": statistics.median(values),
            "min": min(values),
            "max": max(values),
        }
    return {
        "repeat_mode_set_consistent": all(sig == signatures[0] for sig in signatures),
        "timing": timing,
        "runs": runs,
    }


def _profile_representative(output_dir: Path, case_name: str,
                            branch: str, n_init: int) -> int:
    spec = next(item for item in CASES if item["name"] == case_name)
    factory = ScenarioProviderFactory(cfg)
    _, n3d, _ = _build_case(factory, spec)
    params = {**cfg.P2P_3D, "seed_source": "internal_scan",
              "seed_betas_deg": [], "n_init": n_init, "jacobian": "central",
              "candidate_nearest": int(spec.get("candidate_nearest", 1))}
    rx = (float(spec["range_km"]), 0.0)
    n_ah = RefractiveIndexAH3D(
        n3d, float(spec["freq_MHz"]), branch, geomag=cfg.GEOMAG,
        integration_profile="smoke")
    betas, _ = p2p_ah._seed_betas_3d(rx, params, None)
    trace_cache = {}
    trace_stats = {
        "cache_hits": 0, "cache_misses": 0, "avoided_ray_calls": 0,
        "profile_wall_s": {}, "profile_nfev": {},
        "profile_accepted_steps": {}, "reference_ray_calls": 0,
        "reference_nfev": 0, "reference_accepted_steps": 0,
        "print_progress": True,
    }
    candidate_diag = {}
    fan_started = time.perf_counter()
    seeds, audit_count = p2p_ah._candidate_seeds(
        cfg.TX_POS, rx, n_ah, float(spec["freq_MHz"]), betas, None, 0.0,
        nearest=int(spec.get("candidate_nearest", 1)),
        trace_cache=trace_cache, trace_stats=trace_stats,
        diagnostics=candidate_diag)
    candidate_fan_wall_s = time.perf_counter() - fan_started
    seed = next((item for item in seeds if item["provenance"] == "crossing"),
                seeds[0] if seeds else None)
    if seed is None:
        raise RuntimeError(f"no representative candidate for {case_name}/{branch}")

    profiler = cProfile.Profile()
    profiler.enable()
    solved = p2p_ah._shoot_ah_3d(
        cfg.TX_POS, rx, n_ah, float(spec["freq_MHz"]),
        seed["beta_deg"], seed["az_deg"],
        max_iter=int(params.get("max_iter", 12)),
        tol_km=float(params.get("tol_km", 5.0)),
        d_beta=float(params.get("d_beta", 0.05)),
        d_az=float(params.get("d_az", 0.10)), jacobian="central", t=0.0,
        trace_cache=trace_cache, trace_stats=trace_stats, seed_info=seed)
    profiler.disable()
    stats = pstats.Stats(profiler).strip_dirs().sort_stats("cumulative")
    entries = []
    for (filename, line, function), values in list(stats.stats.items()):
        cc, nc, total, cumulative, _ = values
        entries.append({"function": f"{filename}:{line}:{function}",
                        "primitive_calls": cc, "calls": nc,
                        "total_s": total, "cumulative_s": cumulative})
    entries.sort(key=lambda item: item["cumulative_s"], reverse=True)

    def cumulative_named(*names):
        return max((item["cumulative_s"] for item in entries
                    if any(name in item["function"] for name in names)), default=0.0)

    record = {
        "case": case_name, "branch": branch, "n_init": n_init,
        "representative_seed": seed, "candidate_audit_rays": audit_count,
        "candidate_fan_wall_s": candidate_fan_wall_s,
        "candidate_diagnostics": candidate_diag,
        "trace_stats": trace_stats, "formal_mode_passed": solved is not None,
        "root_cumulative_s": cumulative_named(":root", "appleton_hartree_roots"),
        "provider_query_cumulative_s": cumulative_named(":Ne_at", ":get_Ne"),
        "top20_cumulative": entries[:20],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "a3r_p0r_profile.json"
    report = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    report[f"{case_name}:{branch}"] = record
    _atomic_write(path, report)
    print(f"A3R-P0R profile output={path}", flush=True)
    return 0 if solved is not None else 1


def _summarize_checkpoints(output_dir: Path) -> int:
    cases = {}
    hashes = {}
    for spec in CASES:
        path = output_dir / f"a3r_p0r_{spec['name']}.json"
        if not path.exists():
            continue
        report = json.loads(path.read_text(encoding="utf-8"))
        case = report.get("cases", {}).get(spec["name"])
        if case is not None:
            cases[spec["name"]] = case
            hashes[spec["name"]] = report.get("benchmark_config_hash")

    required_repeats = {spec["name"]: (1 if spec["scenario"] == "vacuum" else 3)
                        for spec in CASES}
    repeats_complete = all(
        len(cases.get(name, {}).get("runs", [])) >= count
        for name, count in required_repeats.items())
    consistent = all(case.get("repeat_mode_set_consistent", False)
                     for case in cases.values()) and len(cases) == len(CASES)
    formal_modes = [mode for name, case in cases.items()
                    if name != "vacuum_uniform_degenerate"
                    for run in case.get("runs", [])
                    for mode in run.get("modes", [])]
    formal_passed = bool(formal_modes) and all(
        mode.get("convergence_status") == "passed"
        and mode.get("branch_status") == "ground"
                    and mode.get("p2p_jacobian_effective")
                    == "central_reference_multiscale"
        for mode in formal_modes)
    within_hard_limit = all(
        float(run.get("algorithmic_wall_s", float("inf"))) <= 9000.0
        for name, case in cases.items() if name != "vacuum_uniform_degenerate"
        for run in case.get("runs", []))
    warnings = [
        {"case": name, "repeat": i + 1,
         "algorithmic_wall_s": run["algorithmic_wall_s"]}
        for name, case in cases.items() if name != "vacuum_uniform_degenerate"
        for i, run in enumerate(case.get("runs", []))
        if float(run.get("algorithmic_wall_s", 0.0)) > 7200.0
    ]
    runtime_case_names = (
        "vacuum_uniform_degenerate",
        "smooth_layer_single_mode",
        "tid_3d",
    )
    runtime_repeats_complete = all(
        len(cases.get(name, {}).get("runs", [])) >= required_repeats[name]
        for name in runtime_case_names)
    runtime_closure_passed = all((
        runtime_repeats_complete, consistent, formal_passed, within_hard_limit,
    ))
    high_modes = (cases.get("coleman_chapman_high_low", {}).get("runs", [{}])[0]
                  .get("modes", []))
    high_low_demonstrated = any(
        left.get("wave_mode") == right.get("wave_mode")
        and (abs(float(left["h_reflect_km"]) - float(right["h_reflect_km"]))
             >= float(cfg.P2P_3D["clust_h_km"])
             or abs(float(left["tau_ms"]) - float(right["tau_ms"]))
             >= float(cfg.P2P_3D["clust_tau_ms"]))
        for i, left in enumerate(high_modes) for right in high_modes[i + 1:])
    acceptance = {
        "runtime_closure_status": (
            "passed" if runtime_closure_passed else "runtime_blocked"),
        "repeats_complete": repeats_complete,
        "all_repeat_mode_sets_consistent": consistent,
        "all_formal_modes_passed": formal_passed,
        "within_9000_s_hard_limit": within_hard_limit,
        "performance_warnings": warnings,
        "high_low_physics_demonstrated": high_low_demonstrated,
        "high_low_case_status": ("passed" if high_low_demonstrated
                                  else "case_physics_not_demonstrated"),
    }
    acceptance["p0_reference_complete"] = all((
        repeats_complete, consistent, formal_passed, within_hard_limit,
        high_low_demonstrated,
    ))
    acceptance["p0_overall_status"] = (
        "passed" if acceptance["p0_reference_complete"]
        else "case_physics_not_demonstrated" if not high_low_demonstrated
        else "runtime_blocked")
    acceptance["p1_gate"] = (
        "open" if acceptance["p0_reference_complete"] else "closed")
    aggregate = {
        "schema_version": SCHEMA_VERSION,
        "benchmark_config_hash_semantics": "per_case_frozen_hash",
        "case_benchmark_config_hashes": hashes,
        "cases": cases,
        "acceptance": acceptance,
    }
    _atomic_write(output_dir / BENCHMARK_FILENAME, aggregate)
    print(f"A3R-P0R summary={output_dir / BENCHMARK_FILENAME}", flush=True)
    return 0 if acceptance["p0_reference_complete"] else 1


def _audit_high_low_distance(output_dir: Path, n_init: int) -> int:
    spec = next(item for item in CASES if item["name"] == "coleman_chapman_high_low")
    factory = ScenarioProviderFactory(cfg)
    provider, n3d, setup_s = _build_case(factory, spec)
    params = {**cfg.P2P_3D, "seed_source": "internal_scan",
              "seed_betas_deg": [], "n_init": n_init}
    betas, _ = p2p_ah._seed_betas_3d((float(spec["range_km"]), 0.0), params, None)
    fans = {}
    models = {}
    fan_wall_s = {}
    for branch in ("O", "X"):
        model = RefractiveIndexAH3D(
            n3d, float(spec["freq_MHz"]), branch, geomag=cfg.GEOMAG,
            integration_profile="smoke")
        models[branch] = model
        started = time.perf_counter()
        fans[branch] = [p2p_ah.trace_single_ray_ah_3d(
            cfg.TX_POS, float(beta), 0.0, model, float(spec["freq_MHz"]), t=0.0)
            for beta in betas]
        fan_wall_s[branch] = time.perf_counter() - started

    metadata = provider.get_metadata()
    report = {
        "scenario": spec["scenario"], "frequency_MHz": spec["freq_MHz"],
        "n_init": n_init, "distance_domain_km": [600.0, 1600.0],
        "distance_step_km": 100.0, "setup_wall_s": setup_s,
        "candidate_fan_wall_s": fan_wall_s,
        "fan_samples": {
            branch: [{"beta_deg": float(ray["beta_deg_tx"]),
                      "x_land_km": float(ray["x_land_km"]),
                      "branch_status": ray.get("branch_status")}
                     for ray in fans[branch]]
            for branch in ("O", "X")
        },
        "provider_hash": metadata.get("provider_hash"),
        "provider_config_hash": metadata.get("config_hash"),
        "audits": [], "candidate_distances_km": [],
        "selected_distance_km": None,
    }
    path = output_dir / "a3r_p0r_high_low_distance_audit.json"
    for distance in np.arange(600.0, 1600.0 + 1e-9, 100.0):
        branches = {}
        selected = False
        for branch in ("O", "X"):
            diag = {}
            seeds, _ = p2p_ah._candidate_seeds(
                cfg.TX_POS, (float(distance), 0.0), models[branch],
                float(spec["freq_MHz"]), betas, fans[branch], 0.0,
                nearest=int(spec.get("candidate_nearest", 3)), diagnostics=diag)
            physical = [seed for seed in seeds
                        if seed["provenance"] in ("crossing", "fold")]
            branches[branch] = {
                "physical_candidate_count": len(physical),
                "physical_candidates": physical,
                "candidate_diagnostics": diag,
            }
            selected = selected or len(physical) >= 2
        report["audits"].append({"distance_km": float(distance),
                                 "branches": branches})
        if selected and report["selected_distance_km"] is None:
            report["selected_distance_km"] = float(distance)
        if selected:
            report["candidate_distances_km"].append(float(distance))
        output_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write(path, report)
    print(f"A3R-P0R high-low audit={path}", flush=True)
    print(f"selected_distance_km={report['selected_distance_km']}", flush=True)
    return 0 if report["selected_distance_km"] is not None else 1


def _prepare_p0c_reference(output_dir: Path, n_init: int) -> int:
    """Freeze Coleman inputs and independently locate a stable two-root interval."""
    spec = next(item for item in CASES if item["name"] == "coleman_chapman_high_low")
    provider, n3d, setup_s = _build_case(ScenarioProviderFactory(cfg), spec)
    geomag_b0 = {**cfg.GEOMAG, "fH_MHz": 0.0}
    model = RefractiveIndexAH3D(
        n3d, float(spec["freq_MHz"]), "O", geomag=geomag_b0,
        integration_profile="smoke")
    cache = {}

    def trace(beta: float) -> dict:
        key = float(beta)
        if key not in cache:
            cache[key] = p2p_ah.trace_single_ray_ah_3d(
                cfg.TX_POS, key, 0.0, model, float(spec["freq_MHz"]), t=0.0)
        return cache[key]

    def landing(beta: float) -> float:
        ray = trace(beta)
        if ray.get("branch_status") != "ground":
            raise ValueError("non-ground ray inside bounded root interval")
        return float(ray["x_land_km"])

    betas = np.linspace(2.0, 88.0, int(n_init))
    coarse = [trace(float(beta)) for beta in betas]
    valid = [(float(ray["beta_deg_tx"]), float(ray["x_land_km"]))
             for ray in coarse if ray.get("branch_status") == "ground"
             and np.isfinite(ray.get("x_land_km", np.nan))]
    folds = []
    for i in range(1, len(valid) - 1):
        left_slope = (valid[i][1] - valid[i - 1][1]) / (valid[i][0] - valid[i - 1][0])
        right_slope = (valid[i + 1][1] - valid[i][1]) / (valid[i + 1][0] - valid[i][0])
        if left_slope * right_slope < 0.0:
            folds.append((i, left_slope, right_slope))

    result = {
        "stage": "P0C-0/P0C-1", "status": "case_physics_not_demonstrated",
        "reference_semantics": "independent_literature_guided_numerical_root_reference",
        "literature_parameter_audit": {
            "Coleman_2011": {
                "available": ["foF2", "hmF2", "ymF2", "frequency",
                              "representative_low_and_high_launch_angles",
                              "Haselgrove_phase_path_and_apex_tables"],
                "not_exactly_specified": ["the approximately 1100 km receiver range"],
                "use": "freeze field; independently select receiver inside computed two-root interval",
            },
            "Nosikov_2017": {
                "available": ["IRI-2007", "8 MHz", "2014-06-22 12:00 UT",
                              "Kaliningrad-Stockholm", "two high and two low rays"],
                "use": "external literature trend only; not frozen numerical truth",
            },
            "Nosikov_2020": {
                "status": "partial_reproduction_with_assumptions",
                "reason": "simplified-case endpoint and supplementary inputs are incomplete",
            },
        },
        "provider_metadata": provider.get_metadata(), "setup_wall_s": setup_s,
        "frequency_MHz": float(spec["freq_MHz"]), "wave_semantics": "isotropic_B0",
        "branch_semantics": "magnetoionic_degenerate_single_count",
        "coarse_elevation_domain_deg": [2.0, 88.0], "coarse_sample_count": int(n_init),
        "coarse_ground_samples": [{"beta_deg": b, "x_land_km": x} for b, x in valid],
        "fold_candidates": [], "selected_target_distance_km": None,
        "independent_root_brackets_deg": [], "independent_roots_deg": [],
        "p2_full_azimuth_audit_required": True,
        "magnetoionic_caustic_validation_required": True,
    }
    if folds:
        i, _, _ = max(folds, key=lambda item: valid[item[0]][1])
        lo = valid[max(0, i - 1)][0]
        hi = valid[min(len(valid) - 1, i + 1)][0]
        refined_betas = np.linspace(lo, hi, 9)
        refined = sorted((float(beta), landing(float(beta))) for beta in refined_betas)
        slopes = np.diff([x for _, x in refined]) / np.diff([b for b, _ in refined])
        changes = np.where(slopes[:-1] * slopes[1:] < 0.0)[0]
        if changes.size:
            j = int(changes[np.argmax([refined[k + 1][1] for k in changes])]) + 1
            left = refined[:j + 1]
            right = refined[j:]
            overlap_lo = max(min(x for _, x in left), min(x for _, x in right))
            overlap_hi = min(max(x for _, x in left), max(x for _, x in right))
            if overlap_hi > overlap_lo:
                target = overlap_lo + 0.5 * (overlap_hi - overlap_lo)

                def bracket(samples):
                    for a, b in zip(samples[:-1], samples[1:]):
                        if (a[1] - target) * (b[1] - target) <= 0.0:
                            return (a[0], b[0])
                    return None

                brackets = [bracket(left), bracket(right)]
                if all(brackets):
                    roots = [float(brentq(lambda beta: landing(beta) - target,
                                          a, b, xtol=1.0e-5, rtol=1.0e-10))
                             for a, b in brackets]
                    result.update({
                        "status": "p0c_1_passed",
                        "selected_target_distance_km": float(target),
                        "target_fold_margin_fraction": 0.5,
                        "independent_root_brackets_deg": [list(v) for v in brackets],
                        "independent_roots_deg": roots,
                        "root_rays": [{
                            "ray_label": label, "beta_deg": root,
                            "x_land_km": landing(root),
                            "h_reflect_km": float(trace(root).get("h_reflect_km", np.nan)),
                            "branch_status": trace(root).get("branch_status"),
                        } for label, root in zip(("low_ray", "high_ray"), roots)],
                    })
            result["fold_candidates"].append({
                "coarse_beta_deg": valid[i][0], "refined_samples": [
                    {"beta_deg": b, "x_land_km": x} for b, x in refined]})

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "a3r_p0c_reference_preparation.json"
    result["reference_config_hash"] = _stable_hash({
        "provider": provider.get_metadata().get("config_hash"),
        "frequency_MHz": result["frequency_MHz"],
        "elevation_domain_deg": result["coarse_elevation_domain_deg"],
        "sample_count": result["coarse_sample_count"],
        "target_distance_km": result["selected_target_distance_km"],
    })
    _atomic_write(path, result)
    print(f"A3R-P0C preparation={path}", flush=True)
    print(f"status={result['status']}", flush=True)
    return 0 if result["status"] == "p0c_1_passed" else 1


def _p0c_trace_record(ray: dict, label: str, profile: str) -> dict:
    keys = ("beta_deg_tx", "az_deg_tx", "x_land_km", "y_land_km",
            "h_reflect_km", "geom_path_km", "group_path_km", "tau_ms",
            "branch_status", "convergence_status", "hamiltonian_residual_max",
            "nfev", "accepted_steps", "integration_profile")
    record = {"ray_label": label, "profile": profile}
    for key in keys:
        value = ray.get(key)
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, float) and not np.isfinite(value):
            value = None
        record[key] = value
    return record


def _close_p0c(output_dir: Path, repeats: int, n_init: int) -> int:
    """Run the independent P0C roots and match them to production P2P."""
    prep_path = output_dir / "a3r_p0c_reference_preparation.json"
    prep = json.loads(prep_path.read_text(encoding="utf-8"))
    if prep.get("status") != "p0c_1_passed":
        raise RuntimeError("P0C-1 preparation is not passed")
    spec = next(item for item in CASES if item["name"] == "coleman_chapman_high_low")
    target = float(prep["selected_target_distance_km"])
    roots = prep["independent_roots_deg"]
    config = {"schema_version": SCHEMA_VERSION, "prep_hash": prep["reference_config_hash"],
              "target_distance_km": target, "roots_deg": roots, "repeats": int(repeats),
              "n_init": int(n_init), "wave_mode": "O", "geomag": "B=0",
              "formal_profiles": ["research", "reference"], "tol_km": 5.0,
              "tightened": {"root_xtol_deg": 5e-6, "trace_profile": "reference"}}
    config_hash = _stable_hash(config)
    path = output_dir / "a3r_p0c_closure.json"
    report = {"schema_version": SCHEMA_VERSION, "stage": "P0C-2/P0C-3",
              "configuration": config, "configuration_hash": config_hash,
              "repeats": []}
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        if old.get("configuration_hash") != config_hash:
            raise RuntimeError("existing P0C closure checkpoint has a different config hash")
        report = old
    provider, n3d, setup_s = _build_case(ScenarioProviderFactory(cfg), spec)
    model = RefractiveIndexAH3D(n3d, float(spec["freq_MHz"]), "O",
                                geomag={**cfg.GEOMAG, "fH_MHz": 0.0},
                                integration_profile="research")
    if "root_tolerance_convergence" not in report:
        cache = {}

        def landing(beta: float) -> float:
            key = float(beta)
            if key not in cache:
                model.integration_profile = "smoke"
                cache[key] = p2p_ah.trace_single_ray_ah_3d(
                    cfg.TX_POS, key, 0.0, model, float(spec["freq_MHz"]), t=0.0)
            return float(cache[key]["x_land_km"])

        root_checks = []
        for label, bracket in zip(("low_ray", "high_ray"),
                                  prep["independent_root_brackets_deg"]):
            coarse = float(brentq(lambda beta: landing(beta) - target, *bracket,
                                  xtol=1.0e-5, rtol=1.0e-10))
            tight = float(brentq(lambda beta: landing(beta) - target, *bracket,
                                 xtol=5.0e-6, rtol=1.0e-12))
            root_checks.append({"ray_label": label, "coarse_root_deg": coarse,
                                "tight_root_deg": tight,
                                "root_change_deg": abs(tight - coarse),
                                "tight_landing_residual_km": abs(landing(tight) - target)})
        report["root_tolerance_convergence"] = root_checks
        _atomic_write(path, report)
    for repeat in range(len(report["repeats"]), repeats):
        print("P0C-2 repeat={}/{} start".format(repeat + 1, repeats), flush=True)
        independent = []
        for label, beta in zip(("low_ray", "high_ray"), roots):
            profiles = {}
            for profile in ("research", "reference"):
                model.integration_profile = profile
                ray = p2p_ah.trace_single_ray_ah_3d(
                    cfg.TX_POS, float(beta), 0.0, model, float(spec["freq_MHz"]), t=0.0)
                profiles[profile] = _p0c_trace_record(ray, label, profile)
            research = profiles["research"]
            reference = profiles["reference"]
            independent.append({"ray_label": label, "beta_seed_deg": float(beta),
                                "research": research, "reference": reference,
                                "converged": research.get("convergence_status") == "passed"
                                and reference.get("convergence_status") == "passed"})
        params = {**cfg.P2P_3D, "n_init": int(n_init), "jacobian": "central",
                  "seed_source": "internal_scan", "seed_betas_deg": [],
                  "candidate_nearest": 3}
        diagnostics = {}
        production = p2p_ah.find_all_rays_ah_3d(
            cfg.TX_POS, (target, 0.0), n3d, float(spec["freq_MHz"]), params,
            t=0.0, wave_mode="O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0},
            diagnostics=diagnostics)
        matches = []
        for ref in independent:
            rr = ref["reference"]
            candidates = sorted(production, key=lambda item: (
                abs(float(item["beta_deg_tx"]) - float(ref["beta_seed_deg"])),
                abs(float(item["h_reflect_km"]) - float(rr["h_reflect_km"]))))
            match = candidates[0] if candidates else None
            matches.append({"ray_label": ref["ray_label"], "matched": match is not None,
                            "beta_error_deg": (abs(float(match["beta_deg_tx"]) - float(ref["beta_seed_deg"]))
                                                if match else None),
                            "height_error_km": (abs(float(match["h_reflect_km"]) - float(rr["h_reflect_km"]))
                                                 if match else None),
                            "tau_error_ms": (abs(float(match["tau_ms"]) - float(rr["tau_ms"]))
                                              if match else None),
                            "production": _mode_record(match) if match else None})
        repeat_result = {"repeat": repeat + 1, "setup_wall_s": setup_s,
                         "independent": independent, "production": [_mode_record(m) for m in production],
                         "matches": matches, "diagnostics": diagnostics,
                         "provider_hash": provider.get_metadata().get("provider_hash"),
                         "provider_config_hash": provider.get_metadata().get("config_hash")}
        report["repeats"].append(repeat_result)
        _atomic_write(path, report)
        print("P0C-2 repeat={}/{} done modes={}".format(
            repeat + 1, repeats, len(production)), flush=True)
    all_passed = len(report["repeats"]) >= repeats and all(
        item.get("converged") for run in report["repeats"] for item in run["independent"])
    signatures = [tuple((m["ray_label"], round(float(m["beta_error_deg"]), 6))
                        for m in run["matches"] if m["matched"])
                  for run in report["repeats"]]
    repeat_consistent = len(signatures) == repeats and len(set(signatures)) == 1
    one_to_one = all(len(run["matches"]) == 2 and all(m["matched"] for m in run["matches"])
                     and len({m["production"]["beta_deg_tx"] for m in run["matches"]}) == 2
                     for run in report["repeats"])
    root_tightened = all(item["root_change_deg"] <= 1.0e-5
                         and item["tight_landing_residual_km"] < 0.05
                         for item in report["root_tolerance_convergence"])
    integration_converged = all(
        abs(float(item["research"]["x_land_km"]) - float(item["reference"]["x_land_km"])) < 0.1
        and abs(float(item["research"]["h_reflect_km"]) - float(item["reference"]["h_reflect_km"])) < 0.1
        and abs(float(item["research"]["tau_ms"]) - float(item["reference"]["tau_ms"])) < 1.0e-4
        for run in report["repeats"] for item in run["independent"])
    prior_path = output_dir / BENCHMARK_FILENAME
    prior = json.loads(prior_path.read_text(encoding="utf-8")) if prior_path.exists() else {}
    prior_runtime_passed = prior.get("acceptance", {}).get("runtime_closure_status") == "passed"
    complete = bool(all_passed and repeat_consistent and one_to_one
                    and root_tightened and integration_converged and prior_runtime_passed)
    report["acceptance"] = {
        "p0c_2_independent_research_reference_passed": all_passed,
        "p0c_2_integration_tightening_converged": integration_converged,
        "p0c_2_root_tolerance_tightening_converged": root_tightened,
        "p0c_2_repeat_mode_signature_consistent": repeat_consistent,
        "p0c_2_production_one_to_one_match": one_to_one,
        "p0c_3_prior_vacuum_single_mode_tid_complete": prior_runtime_passed,
        "p0c_3_p0_reference_complete": complete,
        "p0_overall_status": "passed" if complete else "case_physics_not_demonstrated",
        "p1_gate": "open" if complete else "closed",
        "p2_full_azimuth_audit_required": True,
        "magnetoionic_caustic_validation_required": True,
    }
    _atomic_write(path, report)
    print("P0C-3 status={}".format(report["acceptance"]["p0_overall_status"]), flush=True)
    return 0 if report["acceptance"]["p0c_3_p0_reference_complete"] else 1


def run(output_dir: Path, repeats: int, n_init: int, cases: set[str] | None,
        range_km: float | None = None, wave_mode: str = "both") -> int:
    random.seed(0)
    np.random.seed(0)
    selected = [spec for spec in CASES if not cases or spec["name"] in cases]
    if range_km is not None:
        if len(selected) != 1 or selected[0]["name"] != "coleman_chapman_high_low":
            raise ValueError("--range-km is only valid for the high-low case")
        selected = [{**selected[0], "range_km": float(range_km),
                     "case_range_override": True}]
    if not selected:
        raise ValueError("no benchmark cases selected")
    config_snapshot = {
        "schema_version": SCHEMA_VERSION,
        "seed": 0,
        "reference_semantics": "moderate_density_central_numerical_reference",
        "density_rationale": "8 elevation samples bound AH3D formal runtime",
        "wave_modes": ["O", "X"] if wave_mode == "both" else [wave_mode],
        "continuation": False,
        "integration_final_profile": "research",
        "jacobian": "central",
        "launch_elevation_domain_deg": [2.0, 88.0],
        "candidate_scan_elevation_domain_deg": [cfg.RT_3D["beta_min"], cfg.RT_3D["beta_max"]],
        "candidate_scan_azimuth_domain_deg": [0.0, 0.0],
        "candidate_scan_azimuth_limitation": "P0_records_current_solver; P2 must add full azimuth audit",
        "n_init": int(n_init),
        "p2p": {**{key: cfg.P2P_3D[key] for key in sorted(cfg.P2P_3D)},
                "n_init": int(n_init), "jacobian": "central"},
        "cases": selected,
    }
    config_hash = _stable_hash(config_snapshot)
    report = {
        "schema_version": SCHEMA_VERSION,
        "benchmark_config_hash": config_hash,
        "configuration": config_snapshot,
        "machine": _machine_metadata(),
        "algorithmic_timing_semantics": "provider_ready_then_candidate_discovery_plus_P2P",
        "throughput_speedup": None,
        "risks": [
            "moderate density can miss narrow caustic branches",
            "current candidate audit is elevation-domain only at azimuth zero",
            "central numerical reference is not physical truth",
            "memory measurement is excluded because tracemalloc distorts AH3D timing",
            "provider/dependency hash changes invalidate speed comparison",
        ],
        "cases": {},
    }
    factory = ScenarioProviderFactory(cfg)
    output_dir.mkdir(parents=True, exist_ok=True)
    range_suffix = (f"_{int(range_km)}km" if range_km is not None else "")
    mode_suffix = (f"_{wave_mode}" if wave_mode != "both" else "")
    checkpoint = output_dir / (
        f"a3r_p0r_{selected[0]['name']}{range_suffix}{mode_suffix}.json"
        if len(selected) == 1 else BENCHMARK_FILENAME)
    if checkpoint.exists():
        existing = json.loads(checkpoint.read_text(encoding="utf-8"))
        if existing.get("benchmark_config_hash") != config_hash:
            raise RuntimeError("existing P0R checkpoint has a different config hash")
        report = existing
    for spec in selected:
        print("A3R-P0 case={} repeats={} n_init={}".format(
            spec["name"], repeats, n_init), flush=True)
        existing_case = report["cases"].get(spec["name"], {})
        runs = list(existing_case.get("runs", []))
        for repeat in range(len(runs), repeats):
            print(" repeat={}/{} start".format(repeat + 1, repeats), flush=True)
            result = _run_once(factory, spec, n_init, wave_mode)
            runs.append(result)
            print(" repeat={}/{} done wall_s={:.3f} modes={} rays={}".format(
                repeat + 1, repeats, result["algorithmic_wall_s"],
                result["mode_count"], result["ray_calls"]), flush=True)
            report["cases"][spec["name"]] = {
                "specification": spec,
                **_summarize_runs(runs),
            }
            _atomic_write(checkpoint, report)
        if not runs:
            raise RuntimeError(f"no completed repeats for {spec['name']}")
        report["cases"][spec["name"]] = {
            "specification": spec,
            **_summarize_runs(runs),
        }
    report["acceptance"] = {
        "all_repeat_mode_sets_consistent": all(
            value["repeat_mode_set_consistent"] for value in report["cases"].values()),
        "config_hash_required_for_speed_comparison": True,
        "selected_run_complete": all(spec["name"] in report["cases"] for spec in selected),
        "p0_reference_complete": len(report["cases"]) == len(CASES),
    }
    _atomic_write(checkpoint, report)
    print("A3R-P0R output={}".format(checkpoint), flush=True)
    return 0 if report["acceptance"]["all_repeat_mode_sets_consistent"] else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="A3R-P0 AH3D P2P benchmark")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "output2" / "a3r_benchmark")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--n-init", type=int, default=REFERENCE_SEED_COUNT)
    parser.add_argument("--case", action="append", choices=[c["name"] for c in CASES])
    parser.add_argument("--profile-case", choices=[c["name"] for c in CASES])
    parser.add_argument("--profile-branch", choices=("O", "X"), default="O")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--audit-high-low-distance", action="store_true")
    parser.add_argument("--prepare-p0c", action="store_true")
    parser.add_argument("--close-p0c", action="store_true")
    parser.add_argument("--range-km", type=float)
    parser.add_argument("--wave-mode", choices=("O", "X", "both"), default="both")
    args = parser.parse_args()
    if args.repeats < 1 or args.n_init < 3:
        parser.error("repeats must be >=1 and n-init must be >=3")
    if args.profile_case:
        return _profile_representative(
            args.output_dir, args.profile_case, args.profile_branch, args.n_init)
    if args.summarize_only:
        return _summarize_checkpoints(args.output_dir)
    if args.prepare_p0c:
        return _prepare_p0c_reference(args.output_dir, args.n_init)
    if args.close_p0c:
        return _close_p0c(args.output_dir, args.repeats, args.n_init)
    if args.audit_high_low_distance:
        return _audit_high_low_distance(args.output_dir, args.n_init)
    return run(args.output_dir, args.repeats, args.n_init,
               None if args.case is None else set(args.case), args.range_km,
               args.wave_mode)


if __name__ == "__main__":
    raise SystemExit(main())
