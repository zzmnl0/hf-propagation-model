"""Run staged A4-C1R high-ray numerical-isolation diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
import models.ray_tracer_ah_3d_ecef as ahcore
import tests._run_a4_c1_prep3 as prep3
from models.ne_provider_geocentric import ChapmanGeocentricNeProvider3D
from tests._run_a4_c1_validation import (
    _image_check,
    _load_pharlap_case,
    _range_height,
    _refined_input,
    _sha256,
    _trajectory_metrics,
    _verify_manifest,
)
from utils import geocentric_to_ecef


BASE = ROOT / "output2" / "a4_pharlap" / "c1"
OUT = BASE / "c1r"
C1_0_MANIFEST = BASE / "admission" / "c1_0_freeze_manifest.json"
C1_1_REPORT = BASE / "same_angle" / "c1_1_validation.json"
C1_1_MANIFEST = BASE / "same_angle" / "c1_1_freeze_manifest.json"
C1_1_RAW = BASE / "same_angle" / "raw"
C1_2_MANIFEST = BASE / "independent_roots" / "c1_2_freeze_manifest.json"
C1_3_MANIFEST = BASE / "attribution" / "c1_3_freeze_manifest.json"
AH_SOLVER = ROOT / "models" / "ray_tracer_ah_3d_ecef.py"
AH_PROVIDER = ROOT / "models" / "ne_provider_geocentric.py"
CONFIG = ROOT / "config.py"
PHARLAP_HELP = Path(r"D:\pharlap\pharlap_4.7.4\src\matlab\raytrace_3d_sp.m")
PHARLAP_MEX = Path(r"D:\pharlap\pharlap_4.7.4\mex\raytrace_3d_sp.mexw64")
ANGLE_DEG = 39.54254824570964
EXPECTED = {
    "trajectory_rms_3d_km": 2.410089444517647,
    "trajectory_max_3d_km": 3.9247869073626145,
    "landing_absolute_error_km": 3.9247869686352432,
}
REFINED_PROFILE = {
    "rtol": 3.0e-10,
    "atol_r_km": 1.0e-8,
    "atol_k": 1.0e-12,
    "ds_min_km": 0.001,
    "ds_max_km": 0.25,
}
INTEGRATION_LIMITS = {"trajectory_rms_3d_km": 0.1, "landing_delta_km": 0.2}


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        value, indent=2, ensure_ascii=True, allow_nan=False,
        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist(),
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _sources(input_path: Path) -> dict:
    return {
        "runner": _sha256(Path(__file__)),
        "refined_input": _sha256(input_path),
        "ah_solver": _sha256(AH_SOLVER),
        "ah_provider": _sha256(AH_PROVIDER),
        "config": _sha256(CONFIG),
        "pharlap_help": _sha256(PHARLAP_HELP),
        "pharlap_mex": _sha256(PHARLAP_MEX),
        "c1_0_manifest": _sha256(C1_0_MANIFEST),
        "c1_1_manifest": _sha256(C1_1_MANIFEST),
        "c1_2_manifest": _sha256(C1_2_MANIFEST),
        "c1_3_manifest": _sha256(C1_3_MANIFEST),
    }


def _manifest(stage: int, status: str, sources: dict, artifacts: list[Path],
              checks: dict) -> dict:
    stage_dir = OUT / f"c1r_{stage}"
    result = {
        "schema_version": f"a4-c1r-{stage}-freeze-v1",
        "status": status,
        "created_utc": _utc(),
        "semantics": (
            "A4-C1R validation-only high-ray numerical isolation; one-way; "
            "10 MHz; B=0; collision=0; no root solving; not C1 acceptance"
        ),
        "source_sha256": sources,
        "checks": checks,
        "frozen_artifacts": [
            {"path": str(path.resolve()), "sha256": _sha256(path)}
            for path in artifacts
        ],
    }
    _atomic_json(stage_dir / f"c1r_{stage}_freeze_manifest.json", result)
    return result


def _require_stage(stage: int, statuses: set[str]) -> dict:
    path = OUT / f"c1r_{stage}" / f"c1r_{stage}_freeze_manifest.json"
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("status") not in statuses:
        raise RuntimeError(f"C1R-{stage} status blocks the next stage")
    for item in result.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"C1R-{stage} frozen artifact drift: {artifact}")
    return result


def _short_source(sources: dict) -> str:
    return hashlib.sha256(json.dumps(
        sources, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()[:12]


def _path(ray: dict) -> dict:
    keys = (
        "trajectory_ecef_km",
        "trajectory_geocentric_lat_lon_height",
        "trajectory_physical_height_km",
        "trajectory_geom_path_km",
    )
    return {key: np.asarray(ray[key], dtype=float) for key in keys}


def _frozen_ah() -> tuple[dict, dict]:
    summary_path = C1_1_RAW / "ah3d_high_pharlap_anchor_summary.json"
    npz_path = C1_1_RAW / "ah3d_high_pharlap_anchor_raw.npz"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    with np.load(npz_path, allow_pickle=False) as data:
        arrays = {
            key[len("repeat1_"):]: np.asarray(data[key], dtype=float)
            for key in data.files if key.startswith("repeat1_")
        }
    return summary["runs"][0], arrays


def _pharlap() -> tuple[dict, dict]:
    summary = json.loads(
        (C1_1_RAW / "c1_same_angle_pharlap_summary.json").read_text(encoding="utf-8"))
    case = next(item for item in summary["cases"]
                if item["case_id"] == "high_pharlap_anchor")
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    return case["runs"][0], _load_pharlap_case(
        Path(case["raw_mat"]), context["radius"])[0]


def _trace_cache(stage_dir: Path, label: str, context: dict, provider,
                 profile: str, derivative: str, sources: dict) -> tuple[dict, dict]:
    raw = stage_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    summary_path = raw / f"{label}_summary.json"
    npz_path = raw / f"{label}_raw.npz"
    signature = {
        **sources, "label": label, "profile": profile, "derivative": derivative,
        "provider_hash": provider.get_metadata().get("provider_hash",
                                                     provider.get_metadata().get("grid_hash")),
    }
    if summary_path.is_file() and npz_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if (summary.get("source_sha256") == signature
                and summary.get("raw_npz_sha256") == _sha256(npz_path)):
            with np.load(npz_path, allow_pickle=False) as data:
                ray = {key: np.asarray(data[key], dtype=float) for key in data.files}
            ray.update(summary["ray_metadata"])
            return summary["record"], ray

    original_profile = prep3.PROFILE
    original_gradient = ahcore.spatial_gradient_ah3d_ecef
    local_context = {**context, "provider": provider}
    try:
        prep3.PROFILE = profile
        if derivative == "analytic":
            ahcore.spatial_gradient_ah3d_ecef = _analytic_spatial_gradient
        record, ray = prep3._trace(prep3._model(local_context), local_context, ANGLE_DEG)
    finally:
        ahcore.spatial_gradient_ah3d_ecef = original_gradient
        prep3.PROFILE = original_profile

    arrays = {key: np.asarray(value) for key, value in ray.items()
              if isinstance(value, np.ndarray)}
    metadata = {key: value for key, value in ray.items() if key not in arrays}
    np.savez_compressed(npz_path, **arrays)
    _atomic_json(summary_path, {
        "schema_version": "a4-c1r-single-high-ray-v1",
        "status": "completed_validation_only_trace",
        "semantics": (
            f"same high PHaRLAP anchor; {profile}; {derivative}; "
            "diagnostic only; no root solving"
        ),
        "record": record,
        "ray_metadata": metadata,
        "source_sha256": signature,
        "raw_npz_sha256": _sha256(npz_path),
    })
    return record, ray


def _analytic_provider(input_path: Path, radius: float):
    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    return ChapmanGeocentricNeProvider3D(
        radius,
        foF2_MHz=float(np.asarray(data["foF2_MHz"]).item()),
        hmF2_km=float(np.asarray(data["hmF2_km"]).item()),
        ymF2_km=float(np.asarray(data["ymF2_km"]).item()),
    )


def _analytic_spatial_gradient(state, model, *, t: float = 0.0):
    if not isinstance(model.provider, ChapmanGeocentricNeProvider3D):
        raise TypeError("analytic C1R gradient requires Chapman diagnostic Provider")
    position = np.asarray(state, dtype=float)[:3]
    radius = float(np.linalg.norm(position))
    height = radius - model.earth_radius_km
    provider = model.provider
    u = (height - provider.hmF2_km) / provider.ymF2_km
    ne = provider.get_Ne_ecef(position, model.earth_radius_km, t=t)
    q = np.exp(np.clip(-u, -700.0, 700.0))
    dne_dh = ne * 0.5 * (-1.0 + q) / provider.ymF2_km
    grad_ne = dne_dh * position / radius
    f_hz = model.omega / (2.0 * np.pi)
    gradient = model.k0**2 * cfg.K_FP**2 / f_hz**2 * grad_ne
    model.fd_evaluation_count += 1
    model.fd_selected_platform_counts[(0.0,)] = (
        model.fd_selected_platform_counts.get((0.0,), 0) + 1)
    return gradient, {
        "steps_km": [], "selected_steps_km": [],
        "relative_spread": 0.0, "absolute_spread": 0.0,
        "stability_ratio": 0.0, "estimates": gradient[None, :],
        "semantics": "validation_only_analytic_B0_Chapman_D_r",
    }


def _surface_range(record: dict) -> float:
    value = record.get("ground_range_km")
    if value is None or not np.isfinite(float(value)):
        raise RuntimeError("diagnostic ray did not reach ground")
    return float(value)


def _background(axis, context: dict):
    profile = np.asarray(context["source_ne"], dtype=float)[0, 0, :]
    height = np.asarray(context["source_height"], dtype=float)
    x_max = float(context["target_range_km"]) * 1.08
    image = np.tile(profile[:, None], (1, 300))
    axis.imshow(image, origin="lower", aspect="auto", cmap="jet", alpha=0.5,
                extent=(0.0, x_max, height[0], height[-1]))
    axis.set_xlim(0.0, x_max)
    axis.set_ylim(0.0, min(500.0, float(height[-1])))
    axis.set_xlabel("Ground range (km)")
    axis.set_ylabel("Physical height (km)")


def _plot_paths(axis, context: dict, paths: list[tuple[str, dict, str]]):
    _background(axis, context)
    for label, path, style in paths:
        distance, height = _range_height(path, context["tx"], context["radius"])
        axis.plot(distance, height, style, linewidth=1.8, label=label)
    axis.legend(fontsize=8)


def _image(path: Path) -> dict:
    check = _image_check(path)
    if not check["passed"]:
        raise RuntimeError(f"invalid C1R visualization: {path}")
    return check


def run_c1r0() -> str:
    stage_dir = OUT / "c1r_0"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    prerequisite = {
        "c1_0": _verify_manifest(C1_0_MANIFEST, "passed_c1_0_frozen_admission")["status"],
        "c1_1": _verify_manifest(C1_1_MANIFEST, "failed_c1_1_numeric_thresholds")["status"],
        "c1_2": _verify_manifest(
            C1_2_MANIFEST, "completed_c1_2_diagnostic_after_same_angle_failure")["status"],
        "c1_3": _verify_manifest(C1_3_MANIFEST, "failed_pending_c1r")["status"],
    }
    c1 = json.loads(C1_1_REPORT.read_text(encoding="utf-8"))
    metric = next(item for item in c1["metrics"]
                  if item["case_id"] == "high_pharlap_anchor")
    observed = {key: float(metric[key]) for key in EXPECTED}
    signature_ok = all(np.isclose(
        observed[key], EXPECTED[key], rtol=0.0, atol=1.0e-12) for key in EXPECTED)
    checks = {
        "all_prerequisite_manifests_frozen": len(prerequisite) == 4,
        "failure_signature_exact": signature_ok,
        "same_angle_failed_only_trajectory_rms": (
            not metric["gates"]["trajectory_rms_3d_km"]
            and all(value for key, value in metric["gates"].items()
                    if key != "trajectory_rms_3d_km")),
    }
    report_path = stage_dir / "c1r_0_failure_signature.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-0-v1",
        "status": "passed_c1r0_failure_signature_frozen" if all(checks.values())
                  else "blocked_c1r0_source_or_signature_drift",
        "semantics": "frozen failure signature; no ray recomputation",
        "angle_deg": ANGLE_DEG,
        "expected": EXPECTED,
        "observed": observed,
        "prerequisite_status": prerequisite,
        "checks": checks,
        "source_sha256": sources,
    })
    figure_path = stage_dir / "c1r_0_failure_signature_status.png"
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
    labels = ("RMS", "Maximum", "Landing")
    expected = list(EXPECTED.values())
    actual = list(observed.values())
    x = np.arange(3)
    axes[0].bar(x - 0.18, expected, 0.36, label="Frozen expected", color="#4c78a8")
    axes[0].bar(x + 0.18, actual, 0.36, label="Observed", color="#f58518")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Distance (km)")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].barh(
        np.arange(len(checks)), [int(value) for value in checks.values()],
        color=["#187b3f" if value else "#b22222" for value in checks.values()])
    axes[1].set_yticks(np.arange(len(checks)),
                       [key.replace("_", " ") for key in checks])
    axes[1].set_xlim(0, 1.05)
    axes[1].set_xticks((0, 1), ("FAIL", "PASS"))
    figure.suptitle(
        "A4-C1R-0 frozen high-ray failure signature | one-way | 10 MHz | "
        "B=0 | collision=0\n"
        f"same PHaRLAP launch angle | no recomputation | source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    status = json.loads(report_path.read_text(encoding="utf-8"))["status"]
    _manifest(0, status, sources, [report_path, figure_path], checks)
    if not all(checks.values()):
        raise RuntimeError("C1R-0 failed")
    return status


def run_c1r1() -> str:
    _require_stage(0, {"passed_c1r0_failure_signature_frozen"})
    stage_dir = OUT / "c1r_1"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    context = prep3._load_context(input_path)
    cfg.UPDATE2["AH3D_INTEGRATION"]["c1r_refined"] = dict(REFINED_PROFILE)
    reference_record, reference_ray = _trace_cache(
        stage_dir, "reference_fd", context, context["provider"], "reference", "fd", sources)
    refined_record, refined_ray = _trace_cache(
        stage_dir, "refined_fd", context, context["provider"], "c1r_refined", "fd", sources)
    pair, residual = _trajectory_metrics(
        _path(reference_ray), _path(refined_ray), context["radius"])
    landing_delta = abs(_surface_range(reference_record) - _surface_range(refined_record))
    ph_record, ph_path = _pharlap()
    ph_reference, _ = _trajectory_metrics(ph_path, _path(reference_ray), context["radius"])
    ph_refined, _ = _trajectory_metrics(ph_path, _path(refined_ray), context["radius"])
    excluded = (
        pair["trajectory_rms_3d_km"] < INTEGRATION_LIMITS["trajectory_rms_3d_km"]
        and landing_delta < INTEGRATION_LIMITS["landing_delta_km"])
    status = ("passed_c1r1_integration_not_primary"
              if excluded else "stopped_c1r1_integration_candidate")
    checks = {
        "same_angle": True,
        "same_formal_provider": reference_ray["provider_hash"] == refined_ray["provider_hash"],
        "both_ground_single_turning": (
            reference_record["status"] == refined_record["status"] == "ground"
            and reference_record["turning_count"] == refined_record["turning_count"] == 1),
        "both_constraints_valid": (
            reference_record["constraint_ok"] and refined_record["constraint_ok"]),
        "integration_delta_below_diagnostic_limits": excluded,
    }
    report_path = stage_dir / "c1r_1_integration_isolation.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-1-v1",
        "status": status,
        "semantics": (
            "reference-vs-validation-refined DOP853 A/B on identical frozen linear "
            "Provider; diagnostic convergence, not C1 acceptance"
        ),
        "profiles": {
            "reference": cfg.UPDATE2["AH3D_INTEGRATION"]["reference"],
            "c1r_refined": REFINED_PROFILE,
        },
        "pair_metrics": pair,
        "landing_delta_km": landing_delta,
        "pharlap_comparison": {
            "reference_trajectory_rms_km": ph_reference["trajectory_rms_3d_km"],
            "refined_trajectory_rms_km": ph_refined["trajectory_rms_3d_km"],
            "reference_landing_error_km": abs(
                _surface_range(reference_record) - float(ph_record["ground_range_km"])),
            "refined_landing_error_km": abs(
                _surface_range(refined_record) - float(ph_record["ground_range_km"])),
        },
        "limits": INTEGRATION_LIMITS,
        "records": {"reference": reference_record, "refined": refined_record},
        "checks": checks,
        "source_sha256": sources,
    })
    figure_path = stage_dir / "c1r_1_integration_isolation.png"
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    _plot_paths(axes[0, 0], context, [
        ("DOP853 reference", _path(reference_ray), "k-"),
        ("DOP853 C1R refined", _path(refined_ray), "r--"),
    ])
    axes[0, 0].set_title("Same frozen linear Ne Provider")
    axes[0, 1].plot(residual["u"], residual["total"], color="black", label="3-D")
    axes[0, 1].plot(residual["u"], residual["horizontal"], color="#1368aa",
                    label="horizontal")
    axes[0, 1].plot(residual["u"], np.abs(residual["up"]), color="#c43b3b",
                    label="|up|")
    axes[0, 1].axhline(INTEGRATION_LIMITS["trajectory_rms_3d_km"],
                       color="gray", linestyle=":", label="RMS diagnostic limit")
    axes[0, 1].set(xlabel="Normalized geometric path", ylabel="A/B separation (km)",
                   title="Reference-refined trajectory separation")
    axes[0, 1].legend(fontsize=8)
    values = [pair["trajectory_rms_3d_km"], landing_delta]
    limits = [INTEGRATION_LIMITS["trajectory_rms_3d_km"],
              INTEGRATION_LIMITS["landing_delta_km"]]
    x = np.arange(2)
    axes[1, 0].bar(x - 0.18, values, 0.36, label="Observed", color="#4c78a8")
    axes[1, 0].bar(x + 0.18, limits, 0.36, label="Diagnostic limit", color="#f58518")
    axes[1, 0].set_xticks(x, ("Trajectory RMS", "Landing"))
    axes[1, 0].set_ylabel("Distance (km)")
    axes[1, 0].legend()
    labels = ("Accepted steps", "RHS evaluations", "Hamiltonian residual x1e10")
    for index, (name, record) in enumerate((
            ("reference", reference_record), ("refined", refined_record))):
        axes[1, 1].plot(
            np.arange(3),
            [record["accepted_steps"], record["nfev"],
             1e10 * record["hamiltonian_residual_max"]],
            marker="o", label=name)
    axes[1, 1].set_xticks(np.arange(3), labels, rotation=15)
    axes[1, 1].set_title("Solver diagnostics (mixed scales)")
    axes[1, 1].legend()
    figure.suptitle(
        "A4-C1R-1 integration isolation | high ray | same launch angle | one-way | "
        "10 MHz | B=0\n"
        f"DOP853 full Hamiltonian | frozen linear Ne | source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    _manifest(1, status, sources, [report_path, figure_path] + sorted(
        (stage_dir / "raw").glob("*")), checks)
    if not excluded:
        raise RuntimeError("C1R-1 identified integration as a candidate; stop")
    return status


def run_c1r2() -> str:
    _require_stage(1, {"passed_c1r1_integration_not_primary"})
    stage_dir = OUT / "c1r_2"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    context = prep3._load_context(input_path)
    analytic = _analytic_provider(input_path, context["radius"])
    analytic_record, analytic_ray = _trace_cache(
        stage_dir, "analytic_provider_fd", context, analytic, "research", "fd", sources)
    frozen_record, frozen_path = _frozen_ah()
    pair, residual = _trajectory_metrics(
        frozen_path, _path(analytic_ray), context["radius"])
    landing_delta = abs(_surface_range(frozen_record) - _surface_range(analytic_record))
    ph_record, ph_path = _pharlap()
    formal_ph, _ = _trajectory_metrics(ph_path, frozen_path, context["radius"])
    analytic_ph, _ = _trajectory_metrics(ph_path, _path(analytic_ray), context["radius"])

    heights = np.asarray(context["source_height"], dtype=float)
    linear_ne = np.array([
        context["provider"].get_Ne_geocentric(
            context["tx"][0], context["tx"][1], height) for height in heights])
    analytic_ne = np.array([
        analytic.get_Ne_geocentric(
            context["tx"][0], context["tx"][1], height) for height in heights])
    ecef = geocentric_to_ecef(
        np.full_like(heights, context["tx"][0]),
        np.full_like(heights, context["tx"][1]), heights, context["radius"])
    linear_grad = np.array([
        np.linalg.norm(context["provider"].get_grad_Ne_ecef(
            point, context["radius"], step_km=0.5)) for point in ecef])
    u = (heights - analytic.hmF2_km) / analytic.ymF2_km
    analytic_grad = np.abs(
        analytic_ne * 0.5 * (-1.0 + np.exp(-u)) / analytic.ymF2_km)
    mask = analytic_ne > 0.01 * np.max(analytic_ne)
    ne_relative_max = float(np.max(
        np.abs(linear_ne[mask] - analytic_ne[mask]) / analytic_ne[mask]))
    grad_normalized_max = float(np.max(
        np.abs(linear_grad[mask] - analytic_grad[mask]))
        / np.max(analytic_grad[mask]))
    small = (
        pair["trajectory_rms_3d_km"] < 0.25 * EXPECTED["trajectory_rms_3d_km"]
        and landing_delta < 0.25 * EXPECTED["landing_absolute_error_km"])
    status = ("passed_c1r2_frozen_field_not_primary"
              if small else "completed_c1r2_field_effect_unresolved")
    checks = {
        "formal_path_is_frozen_c1_1": True,
        "analytic_provider_diagnostic_only": True,
        "same_equation_profile_and_angle": True,
        "analytic_trace_ground_single_turning": (
            analytic_record["status"] == "ground"
            and analytic_record["turning_count"] == 1),
        "field_representation_effect_below_quarter_failure": small,
    }
    report_path = stage_dir / "c1r_2_field_representation_isolation.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-2-v1",
        "status": status,
        "semantics": (
            "frozen refined_1km linear Provider vs same-parameter analytic Chapman "
            "Provider; different field representation; diagnostic only"
        ),
        "field_metrics": {
            "significant_region_ne_relative_error_max": ne_relative_max,
            "gradient_error_normalized_to_peak_gradient": grad_normalized_max,
        },
        "trajectory_pair_metrics": pair,
        "landing_delta_km": landing_delta,
        "pharlap_comparison": {
            "formal_linear_trajectory_rms_km": formal_ph["trajectory_rms_3d_km"],
            "analytic_trajectory_rms_km": analytic_ph["trajectory_rms_3d_km"],
            "formal_linear_landing_error_km": abs(
                _surface_range(frozen_record) - float(ph_record["ground_range_km"])),
            "analytic_landing_error_km": abs(
                _surface_range(analytic_record) - float(ph_record["ground_range_km"])),
        },
        "record": analytic_record,
        "checks": checks,
        "source_sha256": sources,
    })
    figure_path = stage_dir / "c1r_2_field_representation_isolation.png"
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    axes[0, 0].plot(linear_ne, heights, label="Frozen 1-km linear")
    axes[0, 0].plot(analytic_ne, heights, "--", label="Analytic Chapman")
    axes[0, 0].set(xlabel="Ne (m^-3)", ylabel="Height (km)",
                   title="Electron-density representation")
    axes[0, 0].legend()
    axes[0, 1].plot(linear_grad, heights, label="Frozen Provider FD gradient")
    axes[0, 1].plot(analytic_grad, heights, "--", label="Analytic radial gradient")
    axes[0, 1].set(xlabel="|grad Ne| (m^-3 km^-1)", ylabel="Height (km)",
                   title="Vertical-gradient representation")
    axes[0, 1].legend(fontsize=8)
    _plot_paths(axes[1, 0], context, [
        ("Frozen linear / research", frozen_path, "k-"),
        ("Analytic Chapman / research", _path(analytic_ray), "r--"),
    ])
    axes[1, 0].set_title("Provider A/B trajectories")
    axes[1, 1].plot(residual["u"], residual["total"], color="black", label="3-D")
    axes[1, 1].plot(residual["u"], residual["horizontal"], label="horizontal")
    axes[1, 1].plot(residual["u"], np.abs(residual["up"]), label="|up|")
    axes[1, 1].set(xlabel="Normalized geometric path",
                   ylabel="Linear-analytic separation (km)",
                   title="Field-representation trajectory effect")
    axes[1, 1].legend()
    figure.suptitle(
        "A4-C1R-2 frozen-field interpolation isolation | high ray | one-way | "
        "10 MHz | B=0\n"
        f"same AH3D equation/profile/angle | analytic branch diagnostic only | "
        f"source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    _manifest(2, status, sources, [report_path, figure_path] + sorted(
        (stage_dir / "raw").glob("*")), checks)
    return status


def _gradient_samples(context: dict, provider, ray: dict) -> dict:
    positions = np.asarray(ray["trajectory_ecef_km"], dtype=float)
    wavevectors = np.asarray(ray["wavevector_ecef_km^-1"], dtype=float)
    indices = np.unique(np.linspace(0, len(positions) - 1, 41).astype(int))
    model = prep3._model({**context, "provider": provider})
    fd = []
    exact = []
    original = ahcore.spatial_gradient_ah3d_ecef
    for index in indices:
        state = np.r_[positions[index], wavevectors[index]]
        fd.append(original(state, model)[0])
        exact.append(_analytic_spatial_gradient(state, model)[0])
    fd = np.asarray(fd)
    exact = np.asarray(exact)
    error = np.linalg.norm(fd - exact, axis=1)
    relative = error / np.maximum(np.linalg.norm(exact, axis=1), model.k0**2 * 1e-12)
    return {
        "u": indices / max(len(positions) - 1, 1),
        "fd": fd,
        "analytic": exact,
        "absolute_error": error,
        "relative_error": relative,
    }


def run_c1r3() -> str:
    _require_stage(2, {
        "passed_c1r2_frozen_field_not_primary",
        "completed_c1r2_field_effect_unresolved",
    })
    stage_dir = OUT / "c1r_3"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    context = prep3._load_context(input_path)
    analytic = _analytic_provider(input_path, context["radius"])
    c1r2_dir = OUT / "c1r_2"
    fd_record, fd_ray = _trace_cache(
        c1r2_dir, "analytic_provider_fd", context, analytic, "research", "fd", sources)
    exact_record, exact_ray = _trace_cache(
        stage_dir, "analytic_provider_analytic_gradient", context, analytic,
        "research", "analytic", sources)
    pair, residual = _trajectory_metrics(
        _path(fd_ray), _path(exact_ray), context["radius"])
    landing_delta = abs(_surface_range(fd_record) - _surface_range(exact_record))
    samples = _gradient_samples(context, analytic, fd_ray)
    ph_record, ph_path = _pharlap()
    fd_ph, _ = _trajectory_metrics(ph_path, _path(fd_ray), context["radius"])
    exact_ph, _ = _trajectory_metrics(ph_path, _path(exact_ray), context["radius"])
    formal_record, formal_path = _frozen_ah()
    formal_ph, _ = _trajectory_metrics(ph_path, formal_path, context["radius"])
    rms_closure = (
        fd_ph["trajectory_rms_3d_km"] - exact_ph["trajectory_rms_3d_km"])
    landing_closure = (
        abs(_surface_range(fd_record) - float(ph_record["ground_range_km"]))
        - abs(_surface_range(exact_record) - float(ph_record["ground_range_km"])))
    small = (
        pair["trajectory_rms_3d_km"] < INTEGRATION_LIMITS["trajectory_rms_3d_km"]
        and landing_delta < INTEGRATION_LIMITS["landing_delta_km"])
    supported = (
        rms_closure > 0.5 * formal_ph["trajectory_rms_3d_km"]
        and landing_closure > 0.5 * EXPECTED["landing_absolute_error_km"])
    status = (
        "passed_c1r3_spatial_derivative_not_primary" if small
        else "completed_c1r3_spatial_derivative_candidate" if supported
        else "completed_c1r3_spatial_derivative_unresolved"
    )
    checks = {
        "same_analytic_provider_profile_and_angle": True,
        "both_ground_single_turning": (
            fd_record["status"] == exact_record["status"] == "ground"
            and fd_record["turning_count"] == exact_record["turning_count"] == 1),
        "analytic_gradient_validation_only": True,
        "finite_gradient_samples": bool(
            np.isfinite(samples["relative_error"]).all()),
        "derivative_effect_below_diagnostic_limits": small,
    }
    report_path = stage_dir / "c1r_3_spatial_derivative_isolation.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-3-v1",
        "status": status,
        "semantics": (
            "same analytic Chapman Provider with production multiscale centered "
            "D_r versus validation-only analytic B=0 radial D_r"
        ),
        "gradient_metrics": {
            "relative_error_max": float(np.max(samples["relative_error"])),
            "relative_error_median": float(np.median(samples["relative_error"])),
            "absolute_error_max": float(np.max(samples["absolute_error"])),
        },
        "trajectory_pair_metrics": pair,
        "landing_delta_km": landing_delta,
        "pharlap_comparison": {
            "formal_linear_fd_rms_km": formal_ph["trajectory_rms_3d_km"],
            "analytic_provider_fd_rms_km": fd_ph["trajectory_rms_3d_km"],
            "analytic_provider_exact_gradient_rms_km": exact_ph["trajectory_rms_3d_km"],
            "rms_closure_km": rms_closure,
            "landing_closure_km": landing_closure,
        },
        "records": {"finite_difference": fd_record, "analytic": exact_record},
        "checks": checks,
        "source_sha256": sources,
    })
    sample_path = stage_dir / "c1r_3_gradient_samples.npz"
    np.savez_compressed(sample_path, **samples)
    figure_path = stage_dir / "c1r_3_spatial_derivative_isolation.png"
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    _plot_paths(axes[0, 0], context, [
        ("Multiscale centered D_r", _path(fd_ray), "k-"),
        ("Analytic Chapman D_r", _path(exact_ray), "r--"),
    ])
    axes[0, 0].set_title("Same analytic Chapman field")
    axes[0, 1].semilogy(
        samples["u"], np.maximum(samples["relative_error"], 1e-18), marker="o")
    axes[0, 1].set(xlabel="Normalized sampled path",
                   ylabel="Relative |D_r(FD)-D_r(exact)|",
                   title="Hamiltonian spatial-gradient agreement")
    axes[0, 1].grid(alpha=0.25)
    axes[1, 0].plot(residual["u"], residual["total"], color="black", label="3-D")
    axes[1, 0].plot(residual["u"], residual["horizontal"], label="horizontal")
    axes[1, 0].plot(residual["u"], np.abs(residual["up"]), label="|up|")
    axes[1, 0].set(xlabel="Normalized geometric path",
                   ylabel="FD-analytic trajectory separation (km)",
                   title="Derivative-only trajectory effect")
    axes[1, 0].legend()
    labels = ("Formal linear FD", "Analytic field FD", "Analytic field exact D_r")
    values = (
        formal_ph["trajectory_rms_3d_km"],
        fd_ph["trajectory_rms_3d_km"],
        exact_ph["trajectory_rms_3d_km"],
    )
    axes[1, 1].bar(np.arange(3), values, color=("#4c78a8", "#f58518", "#54a24b"))
    axes[1, 1].axhline(2.0, color="#b22222", linestyle=":", label="C1 RMS gate")
    axes[1, 1].set_xticks(np.arange(3), labels, rotation=12)
    axes[1, 1].set_ylabel("RMS to frozen PHaRLAP trajectory (km)")
    axes[1, 1].set_title("Attribution context; not a formal same-field pass")
    axes[1, 1].legend()
    figure.suptitle(
        "A4-C1R-3 Hamiltonian derivative isolation | high ray | one-way | "
        "10 MHz | B=0\n"
        f"analytic D_r is validation-only | source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    _manifest(3, status, sources, [report_path, sample_path, figure_path] + sorted(
        (stage_dir / "raw").glob("*")), checks)
    return status


def run_c1r4() -> str:
    stage3 = _require_stage(3, {
        "passed_c1r3_spatial_derivative_not_primary",
        "completed_c1r3_spatial_derivative_candidate",
        "completed_c1r3_spatial_derivative_unresolved",
    })
    stage_dir = OUT / "c1r_4"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    context = prep3._load_context(input_path)
    _, ph_path = _pharlap()
    _, formal_path = _frozen_ah()
    _, residual = _trajectory_metrics(ph_path, formal_path, context["radius"])
    checkpoints = np.array((0.0, 0.25, 0.5, 0.75, 1.0))
    sampled = np.interp(checkpoints, residual["u"], residual["total"])
    help_text = PHARLAP_HELP.read_text(encoding="utf-8", errors="replace")
    visible = {
        "haselgrove_hamiltonian_documented": "Haselgrove" in help_text,
        "spherical_earth_user_radius_documented": (
            "Spherical Earth coordinate system" in help_text
            and "radius of Earth" in help_text),
        "grid_smoothness_required": "gradients must be smooth" in help_text,
        "ode_tolerance_and_step_controls_documented": (
            "ODE solver tolerence" in help_text
            and "minimum step size" in help_text
            and "maximum step size" in help_text),
        "internal_interpolation_formula_visible": False,
        "internal_ode_algorithm_visible": False,
        "source_level_mex_implementation_visible": False,
    }
    stage1 = json.loads(
        (OUT / "c1r_1" / "c1r_1_integration_isolation.json").read_text(encoding="utf-8"))
    stage2 = json.loads(
        (OUT / "c1r_2" / "c1r_2_field_representation_isolation.json").read_text(
            encoding="utf-8"))
    stage3_report = json.loads(
        (OUT / "c1r_3" / "c1r_3_spatial_derivative_isolation.json").read_text(
            encoding="utf-8"))
    categories = {
        "input": "excluded",
        "coordinate_earth_model": "excluded",
        "ground_event": "excluded_as_primary",
        "DOP853_integration": (
            "excluded_as_primary" if stage1["status"].startswith("passed")
            else "supported"),
        "frozen_field_representation": (
            "excluded_as_primary" if stage2["status"].startswith("passed")
            else "unresolved"),
        "Hamiltonian_spatial_derivative": (
            "excluded_as_primary" if stage3_report["status"].startswith("passed")
            else "supported" if stage3_report["status"].endswith("candidate")
            else "unresolved"),
        "PHaRLAP_internal_implementation": "unobservable_external_implementation",
    }
    status = "completed_c1r4_remaining_implementation_attribution"
    checks = {
        "residual_present_before_ground_event": bool(sampled[1] > 0.5),
        "pharlap_public_help_hashed": PHARLAP_HELP.is_file(),
        "pharlap_mex_hashed": PHARLAP_MEX.is_file(),
        "unobservable_internals_not_inferred": all(
            not visible[key] for key in (
                "internal_interpolation_formula_visible",
                "internal_ode_algorithm_visible",
                "source_level_mex_implementation_visible")),
        "no_repeated_terminal_event_experiment": True,
    }
    report_path = stage_dir / "c1r_4_remaining_attribution.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-4-v1",
        "status": status,
        "semantics": (
            "reuse frozen path residual and audit only PHaRLAP public help/MEX "
            "visibility; no PHaRLAP rerun and no inference about hidden internals"
        ),
        "path_residual_checkpoints": {
            f"u={u:.2f}": float(value) for u, value in zip(checkpoints, sampled)
        },
        "pharlap_public_interface": visible,
        "categories": categories,
        "checks": checks,
        "source_sha256": sources,
        "stage3_manifest_status": stage3["status"],
    })
    figure_path = stage_dir / "c1r_4_remaining_attribution.png"
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.8), constrained_layout=True)
    axes[0].plot(residual["u"], residual["total"], color="black")
    axes[0].scatter(checkpoints, sampled, color="#d62728", zorder=3)
    axes[0].axvline(0.5, color="gray", linestyle=":")
    axes[0].set(xlabel="Normalized geometric path", ylabel="AH3D-PHaRLAP 3-D residual (km)",
                title="Residual accumulates before ground event")
    labels = list(visible)
    values = [1 if visible[key] else 0 for key in labels]
    axes[1].barh(np.arange(len(labels)), values,
                 color=["#4c78a8" if value else "#bab0ac" for value in values])
    axes[1].set_yticks(np.arange(len(labels)), [key.replace("_", " ") for key in labels])
    axes[1].set_xlim(0, 1.05)
    axes[1].set_xticks((0, 1), ("NOT VISIBLE", "DOCUMENTED"))
    axes[1].set_title("PHaRLAP public-interface observability")
    figure.suptitle(
        "A4-C1R-4 event/external-implementation attribution | high ray | one-way | "
        "10 MHz | B=0\n"
        f"frozen evidence only | PHaRLAP MEX internals unobservable | "
        f"source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    _manifest(4, status, sources, [report_path, figure_path], checks)
    return status


def run_c1r5() -> str:
    _require_stage(4, {"completed_c1r4_remaining_implementation_attribution"})
    stage_dir = OUT / "c1r_5"
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    sources = _sources(input_path)
    stage4 = json.loads(
        (OUT / "c1r_4" / "c1r_4_remaining_attribution.json").read_text(
            encoding="utf-8"))
    categories = stage4["categories"]
    defect_candidates = [
        key for key, value in categories.items() if value == "supported"
        and key in {"DOP853_integration", "frozen_field_representation",
                    "Hamiltonian_spatial_derivative"}
    ]
    if defect_candidates:
        status = "completed_c1r5_ah3d_defect_candidate_requires_repair_plan"
        decision = (
            "Stop before production changes; prepare a minimal repair plan and "
            "rerun C0R-light plus C1-0..C1-3 after separate confirmation."
        )
    else:
        status = "completed_c1r5_software_difference_or_unresolved"
        decision = (
            "No reproducible AH3D defect isolated; keep failed_pending_c1r, "
            "do not modify production AH3D, and keep C1-4/C2 blocked."
        )
    checks = {
        "all_hypotheses_have_state": all(
            value in {"supported", "excluded", "excluded_as_primary", "unresolved",
                      "unobservable_external_implementation"}
            for value in categories.values()),
        "no_production_source_change_during_c1r": all(
            sources[key] == json.loads(
                (OUT / "c1r_0" / "c1r_0_failure_signature.json").read_text(
                    encoding="utf-8"))["source_sha256"][key]
            for key in ("ah_solver", "ah_provider", "config")),
        "formal_c1_failure_not_overwritten": (
            json.loads(C1_1_REPORT.read_text(encoding="utf-8"))["status"]
            == "failed_c1_1_numeric_thresholds"),
        "decision_boundary_explicit": True,
    }
    report_path = stage_dir / "c1r_5_decision.json"
    _atomic_json(report_path, {
        "schema_version": "a4-c1r-5-v1",
        "status": status,
        "semantics": (
            "C1R hypothesis matrix and decision only; not A4-C1 pass and not "
            "production repair authorization"
        ),
        "hypotheses": categories,
        "ah3d_defect_candidates": defect_candidates,
        "decision": decision,
        "checks": checks,
        "source_sha256": sources,
    })
    figure_path = stage_dir / "c1r_5_hypothesis_matrix.png"
    labels = list(categories)
    states = ("supported", "excluded", "unresolved")
    matrix = np.zeros((len(labels), len(states)), dtype=int)
    for row, label in enumerate(labels):
        value = categories[label]
        column = (0 if value == "supported" else
                  1 if value in {"excluded", "excluded_as_primary"} else 2)
        matrix[row, column] = 1
    from matplotlib.colors import ListedColormap
    figure, axis = plt.subplots(figsize=(12, 6.8), constrained_layout=True)
    axis.imshow(matrix, aspect="auto",
                cmap=ListedColormap(("#ffffff", "#4c78a8")), vmin=0, vmax=1)
    axis.set_xticks(np.arange(3), ("SUPPORTED", "EXCLUDED", "UNRESOLVED/UNOBSERVABLE"))
    axis.set_yticks(np.arange(len(labels)), [label.replace("_", " ") for label in labels])
    for row in range(len(labels)):
        for column in range(3):
            axis.text(column, row, "X" if matrix[row, column] else "",
                      ha="center", va="center", color="white", fontweight="bold")
    axis.set_title(
        "A4-C1R-5 hypothesis evidence matrix | high ray | one-way | 10 MHz | B=0\n"
        f"formal C1 remains failed | no threshold relaxation | source={_short_source(sources)}")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    _image(figure_path)
    _manifest(5, status, sources, [report_path, figure_path], checks)
    if not all(checks.values()):
        raise RuntimeError("C1R-5 machine closure failed")
    return status


RUNNERS = {
    0: run_c1r0,
    1: run_c1r1,
    2: run_c1r2,
    3: run_c1r3,
    4: run_c1r4,
    5: run_c1r5,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("all", "0", "1", "2", "3", "4", "5"),
                        default="all")
    args = parser.parse_args()
    stages = range(6) if args.stage == "all" else (int(args.stage),)
    for stage in stages:
        print(f"C1R-{stage}: {RUNNERS[stage]()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
