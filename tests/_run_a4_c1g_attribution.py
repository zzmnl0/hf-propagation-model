"""Run A4-C1G spherical Chapman high-ray attribution."""

from __future__ import annotations

import argparse
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
from numpy.polynomial.legendre import leggauss
from scipy.integrate import IntegrationWarning, cumulative_trapezoid, quad
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from tests import _run_a4_c1_prep3 as prep3
from tests import _run_a4_c1q_attribution as c1q
from tests._run_a4_c1_validation import (
    CASES, _image_check, _refined_input, _resample, _sha256,
)
from utils import geocentric_to_ecef


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1g"
PLAN = ROOT / "word" / "a4_c1g_spherical_high_ray_equation_and_parameter_attribution_plan.txt"
REPORT = ROOT / "word" / "a4_c1g_spherical_high_ray_equation_and_parameter_attribution_report.txt"
FINAL_C1 = ROOT / "output2" / "a4_pharlap" / "c1" / "final" / "a4_c1_validation.json"
C1Q_FINAL = ROOT / "output2" / "a4_pharlap" / "c1" / "c1q" / "c1q_6" / "c1q_6_freeze_manifest.json"
C1WR_FINAL = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr" / "a4_c1wr_manifest.json"
P2PF_FINAL = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr_p2pf" / "a4_c1wr_p2pf_manifest.json"
CASE_META = {case_id: (branch, anchor, angle)
             for case_id, branch, anchor, angle in CASES}


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, value) -> None:
    def clean(item):
        if isinstance(item, dict):
            return {str(key): clean(entry) for key, entry in item.items()}
        if isinstance(item, (list, tuple)):
            return [clean(entry) for entry in item]
        if isinstance(item, np.ndarray):
            return clean(item.tolist())
        if isinstance(item, np.generic):
            return clean(item.item())
        if isinstance(item, float) and not np.isfinite(item):
            return None
        return item

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        clean(value), indent=2, ensure_ascii=True, allow_nan=False,
    ) + "\n", encoding="utf-8")
    temporary.replace(path)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _hash_dict(value: dict) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"),
    ).encode("ascii")).hexdigest()


def _manifest(stage: int, status: str, checks: dict, artifacts: list[Path],
              sources: dict, semantics: str) -> None:
    directory = OUT / f"c1g_{stage}"
    _write_json(directory / f"c1g_{stage}_freeze_manifest.json", {
        "schema_version": f"a4-c1g-{stage}-freeze-v1",
        "status": status,
        "created_utc": _utc(),
        "semantics": semantics,
        "checks": checks,
        "source_sha256": sources,
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in artifacts
        ],
    })


def _load_manifest(path: Path, expected: set[str]) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") not in expected:
        raise RuntimeError(f"unexpected upstream status: {path}")
    for item in value.get("frozen_artifacts", []):
        artifact = Path(item["path"])
        if not artifact.is_file() or _sha256(artifact) != item["sha256"]:
            raise RuntimeError(f"upstream artifact drift: {artifact}")
    return value


def _sources(input_path: Path) -> dict:
    return {
        "runner": _sha256(Path(__file__)),
        "plan": _sha256(PLAN),
        "refined_input": _sha256(input_path),
        "c1_final": _sha256(FINAL_C1),
        "c1q_final_manifest": _sha256(C1Q_FINAL),
        "c1wr_final_manifest": _sha256(C1WR_FINAL),
        "p2pf_final_manifest": _sha256(P2PF_FINAL),
        "pharlap_help": _sha256(
            Path(r"D:\pharlap\pharlap_4.7.4\src\matlab\raytrace_3d_sp.m")),
        "pharlap_mex": _sha256(
            Path(r"D:\pharlap\pharlap_4.7.4\mex\raytrace_3d_sp.mexw64")),
        "ah_solver": _sha256(ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
        "ah_dispersion": _sha256(ROOT / "models" / "dispersion_ah.py"),
        "ah_provider": _sha256(ROOT / "models" / "ne_provider_geocentric.py"),
        "config": _sha256(ROOT / "config.py"),
    }


def run_c1g0() -> str:
    directory = OUT / "c1g_0"
    directory.mkdir(parents=True, exist_ok=True)
    input_path, refined = _refined_input()
    final = json.loads(FINAL_C1.read_text(encoding="utf-8"))
    q = _load_manifest(C1Q_FINAL, {
        "completed_c1q6_pharlap_effective_dynamics_not_recoverable",
        "pharlap_effective_dynamics_not_recoverable",
    })
    wr = _load_manifest(C1WR_FINAL, {"identified_and_repaired_grid_derivative"})
    p2pf = _load_manifest(P2PF_FINAL, {
        "passed_current_ecef_chapman_p2p_root_freeze",
    })
    angles = [item[3] for item in CASES]
    checks = {
        "original_c1_failure_frozen":
            final["status"] == "failed_c1_chapman_high_low_external_comparison",
        "c1q_unresolved_frozen": "pharlap_effective_dynamics_not_recoverable"
            in q["status"],
        "c1wr_grid_derivative_repair_frozen":
            wr["status"] == "identified_and_repaired_grid_derivative",
        "p2pf_current_roots_frozen":
            p2pf["status"] == "passed_current_ecef_chapman_p2p_root_freeze",
        "four_angles_exact": angles == [
            21.94290326504532, 21.94562872092976,
            39.54254824570964, 39.58701138438118,
        ],
        "refined_input_hash_frozen": _sha256(input_path) == refined["sha256"],
        "plan_state_valid": any(
            state in PLAN.read_text(encoding="utf-8")
            for state in (
                "planned_waiting_user_confirmation",
                "completed_identified_spherical_equation_or_invariant_defect",
            )),
    }
    if not all(checks.values()):
        raise RuntimeError(f"C1G-0 admission failed: {checks}")
    payload = {
        "status": "passed_c1g0_frozen_admission",
        "created_utc": _utc(),
        "scope": "B0_spherical_Chapman_same_angle_attribution",
        "angles_deg": angles,
        "checks": checks,
        "upstream": {
            "c1": final["status"], "c1q": q["status"],
            "c1wr": wr["status"], "p2pf": p2pf["status"],
        },
    }
    artifact = directory / "c1g_0_admission.json"
    _write_json(artifact, payload)
    _manifest(0, payload["status"], checks, [artifact], _sources(input_path),
              "freeze_only; no ray recomputation; original C1 failure retained")
    return payload["status"]


def _vertical_profiles(context: dict):
    provider = context["provider"]
    heights = np.asarray(provider.height_km, dtype=float)
    values = np.asarray(provider.Ne_m3, dtype=float)
    vertical = values[values.shape[0] // 2, values.shape[1] // 2]
    cubic = CubicSpline(heights, vertical, bc_type="not-a-knot",
                        extrapolate=False)
    peak = (7.0e6 / 8.98) ** 2

    def linear(height):
        return np.interp(height, heights, vertical, left=0.0, right=0.0)

    def cubic_grid(height):
        value = cubic(height)
        return np.where(np.isfinite(value), np.maximum(value, 0.0), 0.0)

    def analytic(height):
        u = (np.asarray(height, dtype=float) - 350.0) / 100.0
        exponent = 0.5 * (1.0 - u - np.exp(np.clip(-u, -700.0, 700.0)))
        return peak * np.exp(np.clip(exponent, -745.0, 700.0))

    return {"linear": linear, "cubic": cubic_grid, "analytic": analytic}


def _exact_solution(ne_function, radius: float, frequency_MHz: float,
                    *, elevation_deg: float | None = None,
                    impact_km: float | None = None) -> dict:
    f_hz = frequency_MHz * 1e6

    def n2_r(r):
        ne = np.asarray(ne_function(np.asarray(r) - radius), dtype=float)
        return 1.0 - cfg.K_FP**2 * ne / f_hz**2

    n0 = float(np.sqrt(max(float(n2_r(radius)), 0.0)))
    if impact_km is None:
        if elevation_deg is None:
            raise ValueError("elevation_deg or impact_km is required")
        impact_km = n0 * radius * np.cos(np.radians(elevation_deg))
    impact_km = float(impact_km)

    def turning_function(r):
        return np.sqrt(max(float(n2_r(r)), 0.0)) * r - impact_km

    samples = radius + np.linspace(0.0, 462.0, 18481)
    values = np.array([turning_function(value) for value in samples])
    crossings = np.flatnonzero(values[:-1] * values[1:] <= 0.0)
    crossings = crossings[samples[crossings] > radius + 1e-8]
    if crossings.size == 0:
        raise RuntimeError("no ascending turning point for spherical quadrature")
    index = int(crossings[0])
    turning = brentq(turning_function, samples[index], samples[index + 1],
                     xtol=1e-12, rtol=1e-14)
    span = turning - radius
    t_max = np.sqrt(span)

    def radicand(r):
        return float(n2_r(r)) * r * r - impact_km * impact_km

    derivative_step = 1e-4
    derivative = (
        radicand(turning + derivative_step)
        - radicand(turning - derivative_step)
    ) / (2.0 * derivative_step)
    if derivative >= 0.0:
        raise RuntimeError("turning radicand derivative must be negative")

    def components(t):
        t = np.asarray(t, dtype=float)
        r = turning - t * t
        n2 = np.asarray(n2_r(r), dtype=float)
        n = np.sqrt(np.maximum(n2, 0.0))
        rad = np.maximum(n2 * r * r - impact_km * impact_km, 0.0)
        denominator = np.sqrt(rad)
        numerators = np.stack((
            np.full_like(r, impact_km) / r,
            n * r,
            r,
            n2 * r,
        ), axis=-1)
        scale = np.divide(
            2.0 * t, denominator,
            out=np.full_like(t, 2.0 / np.sqrt(-derivative)),
            where=denominator > 1e-13,
        )
        return numerators * scale[..., None]

    quad_values = []
    quad_errors = []
    for component in range(4):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            value, error = quad(
                lambda t: float(components(t)[component]), 0.0, t_max,
                epsabs=1e-10, epsrel=1e-11, limit=500,
            )
        quad_values.append(value)
        quad_errors.append(error)

    nodes, weights = leggauss(1024)
    t_nodes = 0.5 * t_max * (nodes + 1.0)
    gauss_values = 0.5 * t_max * np.sum(
        weights[:, None] * components(t_nodes), axis=0)
    quad_values = np.asarray(quad_values)
    uncertainty = np.abs(quad_values - gauss_values)

    t_ascending = np.linspace(0.0, t_max, 6001)
    transformed = components(t_ascending)
    cumulative = cumulative_trapezoid(
        transformed, t_ascending, axis=0, initial=0.0)
    half = cumulative[-1]
    t_up = t_ascending[::-1]
    up_integrals = cumulative[::-1]
    t_down = t_ascending[1:]
    down_integrals = cumulative[1:]
    r_up = turning - t_up**2
    r_down = turning - t_down**2
    theta = np.r_[half[0] - up_integrals[:, 0],
                  half[0] + down_integrals[:, 0]]
    geom = np.r_[half[1] - up_integrals[:, 1],
                 half[1] + down_integrals[:, 1]]
    group = np.r_[half[2] - up_integrals[:, 2],
                  half[2] + down_integrals[:, 2]]
    phase = np.r_[half[3] - up_integrals[:, 3],
                  half[3] + down_integrals[:, 3]]
    radial = np.r_[r_up, r_down]
    tx = geocentric_to_ecef(*context_global["tx"], radius)
    tx_hat = tx / np.linalg.norm(tx)
    rx = geocentric_to_ecef(*context_global["rx"], radius)
    normal = np.cross(tx_hat, rx / np.linalg.norm(rx))
    normal /= np.linalg.norm(normal)
    tangent = np.cross(normal, tx_hat)
    position = radial[:, None] * (
        np.cos(theta)[:, None] * tx_hat
        + np.sin(theta)[:, None] * tangent
    )
    return {
        "impact_km": impact_km,
        "n0": n0,
        "turning_radius_km": turning,
        "apex_km": turning - radius,
        "ground_range_km": 2.0 * radius * quad_values[0],
        "geom_path_km": 2.0 * quad_values[1],
        "group_path_km": 2.0 * quad_values[2],
        "phase_path_km": 2.0 * quad_values[3],
        "quadrature_uncertainty": {
            "half_angle_rad": uncertainty[0],
            "geom_path_km": 2.0 * uncertainty[1],
            "group_path_km": 2.0 * uncertainty[2],
            "phase_path_km": 2.0 * uncertainty[3],
            "quad_reported_max": max(quad_errors),
        },
        "position_ecef_km": position,
        "trajectory_geom_path_km": geom,
        "trajectory_group_path_km": group,
        "trajectory_phase_path_km": phase,
        "trajectory_height_km": radial - radius,
    }


def _solver_cases(context: dict) -> dict:
    pharlap = c1q._ph_cases(context["radius"])
    result = {}
    for case_id in CASE_META:
        result[(case_id, "PHaRLAP")] = pharlap[case_id]
        result[(case_id, "AH3D")] = c1q._ah_case(case_id, context)
    return result


def _solver_metrics(case: dict, context: dict, ne_function) -> dict:
    position = np.asarray(case["position_ecef_km"], dtype=float)
    p = np.asarray(case["phase_index_vector"], dtype=float)
    radius = context["radius"]
    radial = np.linalg.norm(position, axis=1)
    angular_momentum = np.cross(position, p)
    impact = np.linalg.norm(angular_momentum, axis=1)
    unit_normal = angular_momentum / impact[:, None]
    reference_normal = unit_normal[0]
    plane_angle = np.degrees(np.arccos(np.clip(
        unit_normal @ reference_normal, -1.0, 1.0)))
    ne = ne_function(radial - radius)
    n2 = 1.0 - cfg.K_FP**2 * ne / (context["frequency_MHz"] * 1e6)**2
    dispersion = np.sum(p * p, axis=1) - n2
    radial_hat = position[0] / radial[0]
    tangential = p[0] - np.dot(p[0], radial_hat) * radial_hat
    elevation = np.degrees(np.arctan2(
        np.dot(p[0], radial_hat), np.linalg.norm(tangential)))
    start_hat = position[0] / radial[0]
    end_hat = position[-1] / radial[-1]
    ground_range = radius * np.arccos(np.clip(
        np.dot(start_hat, end_hat), -1.0, 1.0))
    return {
        "initial_elevation_from_state_deg": float(elevation),
        "initial_impact_km": float(impact[0]),
        "impact_relative_drift_max": float(
            np.max(np.abs(impact / impact[0] - 1.0))),
        "plane_normal_drift_max_deg": float(np.max(plane_angle)),
        "dispersion_abs_max": float(np.max(np.abs(dispersion))),
        "ground_range_km": float(ground_range),
        "apex_km": float(np.max(radial - radius)),
        "geom_path_km": float(case["geom_path_km"][-1]),
        "group_path_km": float(case["group_path_km"][-1]),
        "impact_series": impact,
        "plane_angle_series_deg": plane_angle,
        "dispersion_series": dispersion,
    }


def _trajectory_error(case: dict, exact: dict) -> tuple[float, float]:
    count = 401
    u = np.linspace(0.0, 1.0, count)
    solver_geom = np.asarray(case["geom_path_km"], dtype=float)
    solver_position = _resample(
        solver_geom, case["position_ecef_km"], u * solver_geom[-1])
    exact_geom = np.asarray(exact["trajectory_geom_path_km"], dtype=float)
    exact_position = _resample(
        exact_geom, exact["position_ecef_km"], u * exact_geom[-1])
    difference = np.linalg.norm(solver_position - exact_position, axis=1)
    return float(np.sqrt(np.mean(difference**2))), float(np.max(difference))


def _plot_trajectories(rows: list[dict], exact_paths: dict,
                       solver_cases: dict, context: dict, path: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=False)
    for axis, (case_id, branch, anchor, _) in zip(axes.flat, CASES):
        for solver, color in (("PHaRLAP", "#111111"), ("AH3D", "#1769aa")):
            case = solver_cases[(case_id, solver)]
            xyz = c1q._xyz(case["position_ecef_km"], context)
            axis.plot(xyz[:, 0], xyz[:, 2], color=color, lw=1.4,
                      label=solver)
        exact = exact_paths[(case_id, "linear", "nominal", "shared")]
        xyz = c1q._xyz(exact["position_ecef_km"], context)
        axis.plot(xyz[:, 0], xyz[:, 2], color="#c7362f", ls="--", lw=1.3,
                  label="spherical integral")
        axis.set_title(f"{case_id} | {branch} | {anchor} angle")
        axis.set(xlabel="Great-circle range (km)", ylabel="Height (km)")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle(
        "A4-C1G-1 | B=0 spherical Chapman | same nominal angle | "
        "no root solving")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _plot_invariants(rows: list[dict], series: dict, path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 5))
    for (case_id, solver), values in series.items():
        if "high" not in case_id:
            continue
        u = np.linspace(0.0, 1.0, len(values["impact"]))
        label = f"{case_id}/{solver}"
        axes[0].plot(u, values["impact"] / values["impact"][0] - 1.0,
                     label=label)
        axes[1].plot(u, values["plane"], label=label)
        axes[2].plot(u, values["dispersion"], label=label)
    axes[0].set(xlabel="Normalized point index", ylabel="b/b0 - 1",
                title="Bouguer invariant drift")
    axes[1].set(xlabel="Normalized point index", ylabel="Angle (deg)",
                title="Propagation-plane drift")
    axes[2].set(xlabel="Normalized point index", ylabel="|p|^2 - n^2",
                title="Dispersion residual")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    figure.suptitle("A4-C1G-1 spherical invariants | high-ray objects")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _plot_matrix(rows: list[dict], path: Path) -> None:
    selected = [row for row in rows
                if row["profile"] == "linear"
                and row["launch_reference"] == "nominal"]
    labels = [f"{row['case_id']}\n{row['solver']}" for row in selected]
    values = np.array([
        [row["trajectory_rms_km"], abs(row["range_error_km"]),
         abs(row["apex_error_km"]), abs(row["group_path_relative_error"]) * 100]
        for row in selected
    ])
    figure, axis = plt.subplots(figsize=(12, 6))
    image = axis.imshow(np.log10(np.maximum(values, 1e-12)).T,
                        aspect="auto", cmap="viridis")
    axis.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
    axis.set_yticks(range(4), (
        "Trajectory RMS km", "|Range error| km",
        "|Apex error| km", "|Group error| %",
    ))
    for row in range(values.shape[1]):
        for column in range(values.shape[0]):
            axis.text(column, row, f"{values[column, row]:.3g}",
                      ha="center", va="center", color="white", fontsize=8)
    figure.colorbar(image, ax=axis, label="log10(metric)")
    axis.set_title(
        "A4-C1G-1 | solver relative to spherical linear-grid integral")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def run_c1g1() -> str:
    previous = json.loads(
        (OUT / "c1g_0" / "c1g_0_freeze_manifest.json").read_text(
            encoding="utf-8"))
    if previous["status"] != "passed_c1g0_frozen_admission":
        raise RuntimeError("C1G-0 admission is not passed")
    directory = OUT / "c1g_1"
    directory.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    global context_global
    context_global = context
    profiles = _vertical_profiles(context)
    solver_cases = _solver_cases(context)
    exact_paths = {}
    rows = []
    series = {}

    for case_id, branch, anchor, nominal_angle in CASES:
        for profile_name, ne_function in profiles.items():
            nominal = _exact_solution(
                ne_function, context["radius"], context["frequency_MHz"],
                elevation_deg=nominal_angle)
            exact_paths[(case_id, profile_name, "nominal", "shared")] = nominal
        for solver in ("PHaRLAP", "AH3D"):
            case = solver_cases[(case_id, solver)]
            metrics = _solver_metrics(case, context, profiles["linear"])
            series[(case_id, solver)] = {
                "impact": metrics.pop("impact_series"),
                "plane": metrics.pop("plane_angle_series_deg"),
                "dispersion": metrics.pop("dispersion_series"),
            }
            for profile_name, ne_function in profiles.items():
                for launch_reference in ("nominal", "state"):
                    exact_key = (
                        case_id, profile_name, launch_reference,
                        "shared" if launch_reference == "nominal" else solver,
                    )
                    exact = exact_paths.get(exact_key)
                    if exact is None:
                        exact = _exact_solution(
                            ne_function, context["radius"],
                            context["frequency_MHz"],
                            impact_km=metrics["initial_impact_km"])
                        exact_paths[exact_key] = exact
                    rms, maximum = _trajectory_error(case, exact)
                    rows.append({
                        "case_id": case_id,
                        "branch": branch,
                        "anchor": anchor,
                        "solver": solver,
                        "nominal_elevation_deg": nominal_angle,
                        "state_elevation_deg":
                            metrics["initial_elevation_from_state_deg"],
                        "profile": profile_name,
                        "launch_reference": launch_reference,
                        "impact_km": exact["impact_km"],
                        "trajectory_rms_km": rms,
                        "trajectory_max_km": maximum,
                        "range_error_km":
                            metrics["ground_range_km"]
                            - exact["ground_range_km"],
                        "apex_error_km":
                            metrics["apex_km"] - exact["apex_km"],
                        "geom_path_relative_error":
                            metrics["geom_path_km"] / exact["geom_path_km"] - 1.0,
                        "group_path_relative_error":
                            metrics["group_path_km"] / exact["group_path_km"] - 1.0,
                        "impact_relative_drift_max":
                            metrics["impact_relative_drift_max"],
                        "plane_normal_drift_max_deg":
                            metrics["plane_normal_drift_max_deg"],
                        "dispersion_abs_max": metrics["dispersion_abs_max"],
                        "quadrature_range_uncertainty_km":
                            2.0 * context["radius"]
                            * exact["quadrature_uncertainty"]["half_angle_rad"],
                        "quadrature_geom_uncertainty_km":
                            exact["quadrature_uncertainty"]["geom_path_km"],
                        "quadrature_group_uncertainty_km":
                            exact["quadrature_uncertainty"]["group_path_km"],
                    })

    metrics_csv = directory / "c1g_1_spherical_benchmark_metrics.csv"
    _write_csv(metrics_csv, rows)
    arrays = {}
    for (case_id, solver), values in series.items():
        prefix = f"{case_id}_{solver.lower()}"
        arrays[f"{prefix}_impact_km"] = values["impact"]
        arrays[f"{prefix}_plane_angle_deg"] = values["plane"]
        arrays[f"{prefix}_dispersion"] = values["dispersion"]
    samples_npz = directory / "c1g_1_invariant_samples.npz"
    np.savez_compressed(samples_npz, **arrays)

    trajectory_png = directory / "c1g_1_spherical_trajectories.png"
    invariant_png = directory / "c1g_1_spherical_invariants.png"
    matrix_png = directory / "c1g_1_spherical_benchmark_matrix.png"
    _plot_trajectories(rows, exact_paths, solver_cases, context, trajectory_png)
    _plot_invariants(rows, series, invariant_png)
    _plot_matrix(rows, matrix_png)
    image_checks = {
        path.name: _image_check(path)
        for path in (trajectory_png, invariant_png, matrix_png)
    }

    nominal_linear = [
        row for row in rows
        if row["profile"] == "linear"
        and row["launch_reference"] == "nominal"
    ]
    high = [row for row in nominal_linear if row["branch"] == "high_ray"]
    low = [row for row in nominal_linear if row["branch"] == "low_ray"]

    def solver_gate(solver: str, selected: list[dict]) -> bool:
        subset = [row for row in selected if row["solver"] == solver]
        return bool(len(subset) == 2 and all(
            row["trajectory_rms_km"] <= 0.15
            and abs(row["range_error_km"]) <= 0.2
            for row in subset))

    ph_high = solver_gate("PHaRLAP", high)
    ah_high = solver_gate("AH3D", high)
    low_controls = all(
        row["trajectory_rms_km"] <= 0.15
        and abs(row["range_error_km"]) <= 0.2 for row in low)
    unique_solver = (
        "PHaRLAP" if ph_high and not ah_high else
        "AH3D" if ah_high and not ph_high else None
    )
    max_uncertainty = max(
        max(row["quadrature_range_uncertainty_km"],
            row["quadrature_geom_uncertainty_km"],
            row["quadrature_group_uncertainty_km"])
        for row in rows)
    checks = {
        "four_cases_two_solvers_three_profiles_two_launch_references":
            len(rows) == 48,
        "quadrature_uncertainty_le_0_01_km": max_uncertainty <= 0.01,
        "three_valid_images": all(
            value["passed"] for value in image_checks.values()),
        "both_low_controls_pass_linear_nominal": low_controls,
        "exactly_one_solver_closes_both_high": unique_solver is not None,
    }
    if unique_solver and low_controls:
        status = "identified_c1g1_spherical_forward_equation_candidate"
        next_stage = "stop_long_branches_and_confirm_with_state_reference"
    else:
        status = "completed_c1g1_no_unique_spherical_candidate"
        next_stage = "continue_c1g2_launch_geometry_audit"
    summary = {
        "status": status,
        "created_utc": _utc(),
        "scope": "B0_spherical_Chapman_same_angle_no_root_solving",
        "unique_solver_matching_linear_spherical_integral": unique_solver,
        "next_stage": next_stage,
        "checks": checks,
        "max_quadrature_uncertainty_km": max_uncertainty,
        "image_checks": image_checks,
        "best_rows": sorted(nominal_linear,
                            key=lambda row: row["trajectory_rms_km"]),
        "source_tag": _hash_dict(_sources(input_path))[:12],
    }
    summary_json = directory / "c1g_1_spherical_benchmark.json"
    _write_json(summary_json, summary)
    artifacts = [
        summary_json, metrics_csv, samples_npz,
        trajectory_png, invariant_png, matrix_png,
    ]
    _manifest(1, status, checks, artifacts, _sources(input_path),
              "B=0 spherical Chapman; exact integral and invariants; "
              "same nominal angles; no root solving; no production changes")
    return status


def _effective_ne_profile(case: dict, context: dict):
    height = np.asarray(case["height_km"], dtype=float)
    phase_index = np.asarray(case["phase_index_vector"], dtype=float)
    apex_index = int(np.argmax(height))
    ascending_height = height[:apex_index + 1]
    ascending_n2 = np.sum(phase_index[:apex_index + 1] ** 2, axis=1)
    keep = np.r_[True, np.diff(ascending_height) > 1e-8]
    ascending_height = ascending_height[keep]
    ascending_n2 = ascending_n2[keep]
    effective_ne = (
        (1.0 - ascending_n2)
        * (context["frequency_MHz"] * 1e6) ** 2 / cfg.K_FP**2
    )
    interpolator = PchipInterpolator(
        ascending_height, effective_ne, extrapolate=True)

    def query(height_km):
        return np.maximum(interpolator(height_km), 0.0)

    return query, ascending_height, ascending_n2


def _plot_effective_medium(rows: list[dict], samples: dict,
                           solver_cases: dict, exact_paths: dict,
                           context: dict, paths: tuple[Path, Path, Path]) -> None:
    medium_path, trajectory_path, launch_path = paths
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False)
    for axis, (case_id, _, _, _) in zip(axes.flat, CASES):
        values = samples[case_id]
        axis.plot(values["height_km"], values["expected_n2"], color="#1769aa",
                  label="Frozen linear-grid n2")
        axis.plot(values["height_km"], values["effective_n2"], color="#c7362f",
                  ls="--", label="PHaRLAP |p|^2")
        axis.set_title(case_id)
        axis.set(xlabel="Height (km)", ylabel="n^2")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle(
        "A4-C1G-2 | PHaRLAP effective no-field medium vs frozen input")
    figure.tight_layout()
    figure.savefig(medium_path, dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    for axis, (case_id, _, _, _) in zip(axes.flat, CASES):
        case = solver_cases[(case_id, "PHaRLAP")]
        exact = exact_paths[case_id]
        ph_xyz = c1q._xyz(case["position_ecef_km"], context)
        exact_xyz = c1q._xyz(exact["position_ecef_km"], context)
        axis.plot(ph_xyz[:, 0], ph_xyz[:, 2], color="#111111",
                  label="PHaRLAP")
        axis.plot(exact_xyz[:, 0], exact_xyz[:, 2], color="#c7362f",
                  ls="--", label="Integral of PHaRLAP effective n2")
        axis.set_title(case_id)
        axis.set(xlabel="Great-circle range (km)", ylabel="Height (km)")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle(
        "A4-C1G-2 | state-derived effective-medium integral closure")
    figure.tight_layout()
    figure.savefig(trajectory_path, dpi=180)
    plt.close(figure)

    selected = [row for row in rows if row["branch"] == "high_ray"]
    labels = [row["case_id"] for row in selected]
    values = np.array([
        [abs(row["nominal_minus_state_elevation_deg"]),
         row["effective_trajectory_rms_km"],
         abs(row["effective_range_error_km"]),
         row["plasma_term_ratio_median"]]
        for row in selected
    ])
    figure, axes = plt.subplots(1, 4, figsize=(15, 4.5))
    titles = (
        "|Nominal-state angle| deg", "Effective RMS km",
        "|Effective range error| km", "Median Xeff/Xinput",
    )
    for index, axis in enumerate(axes):
        axis.bar(labels, values[:, index], color=("#1769aa", "#c7362f"))
        axis.set_title(titles[index])
        axis.tick_params(axis="x", rotation=20)
        axis.grid(axis="y", alpha=0.25)
    figure.suptitle(
        "A4-C1G-2 | launch exclusion and effective-medium closure")
    figure.tight_layout()
    figure.savefig(launch_path, dpi=180)
    plt.close(figure)


def run_c1g2() -> str:
    previous = json.loads(
        (OUT / "c1g_1" / "c1g_1_freeze_manifest.json").read_text(
            encoding="utf-8"))
    if previous["status"] != "identified_c1g1_spherical_forward_equation_candidate":
        raise RuntimeError("C1G-1 did not identify a unique spherical candidate")
    directory = OUT / "c1g_2"
    directory.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    context = prep3._load_context(input_path)
    global context_global
    context_global = context
    profiles = _vertical_profiles(context)
    solver_cases = _solver_cases(context)
    rows = []
    samples = {}
    exact_paths = {}

    for case_id, branch, anchor, nominal_angle in CASES:
        case = solver_cases[(case_id, "PHaRLAP")]
        effective_ne, height, effective_n2 = _effective_ne_profile(case, context)
        metrics = _solver_metrics(case, context, profiles["linear"])
        exact = _exact_solution(
            effective_ne, context["radius"], context["frequency_MHz"],
            impact_km=metrics["initial_impact_km"])
        exact_paths[case_id] = exact
        rms, maximum = _trajectory_error(case, exact)
        expected_ne = profiles["linear"](height)
        expected_n2 = (
            1.0 - cfg.K_FP**2 * expected_ne
            / (context["frequency_MHz"] * 1e6) ** 2)
        expected_x = 1.0 - expected_n2
        effective_x = 1.0 - effective_n2
        significant = expected_x > 0.05
        ratio = effective_x[significant] / expected_x[significant]
        descending_height = np.asarray(case["height_km"], dtype=float)[
            int(np.argmax(case["height_km"])):]
        descending_n2 = np.sum(
            np.asarray(case["phase_index_vector"], dtype=float)[
                int(np.argmax(case["height_km"])):] ** 2, axis=1)
        order = np.argsort(descending_height)
        descending_height = descending_height[order]
        descending_n2 = descending_n2[order]
        unique = np.r_[True, np.diff(descending_height) > 1e-8]
        descending_interp = PchipInterpolator(
            descending_height[unique], descending_n2[unique],
            extrapolate=False)
        overlap = (
            (height >= descending_height[unique][0])
            & (height <= descending_height[unique][-1]))
        branch_difference = (
            effective_n2[overlap] - descending_interp(height[overlap]))
        row = {
            "case_id": case_id,
            "branch": branch,
            "anchor": anchor,
            "nominal_elevation_deg": nominal_angle,
            "state_elevation_deg":
                metrics["initial_elevation_from_state_deg"],
            "nominal_minus_state_elevation_deg":
                nominal_angle - metrics["initial_elevation_from_state_deg"],
            "effective_trajectory_rms_km": rms,
            "effective_trajectory_max_km": maximum,
            "effective_range_error_km":
                metrics["ground_range_km"] - exact["ground_range_km"],
            "effective_apex_error_km":
                metrics["apex_km"] - exact["apex_km"],
            "effective_geom_relative_error":
                metrics["geom_path_km"] / exact["geom_path_km"] - 1.0,
            "effective_group_relative_error":
                metrics["group_path_km"] / exact["group_path_km"] - 1.0,
            "plasma_term_ratio_median": float(np.median(ratio)),
            "plasma_term_ratio_p05": float(np.percentile(ratio, 5)),
            "plasma_term_ratio_p95": float(np.percentile(ratio, 95)),
            "ascending_descending_n2_rms": float(np.sqrt(np.mean(
                branch_difference**2))),
            "impact_relative_drift_max":
                metrics["impact_relative_drift_max"],
            "quadrature_range_uncertainty_km":
                2.0 * context["radius"]
                * exact["quadrature_uncertainty"]["half_angle_rad"],
        }
        rows.append(row)
        samples[case_id] = {
            "height_km": height,
            "effective_n2": effective_n2,
            "expected_n2": expected_n2,
        }

    metrics_csv = directory / "c1g_2_launch_effective_medium_metrics.csv"
    _write_csv(metrics_csv, rows)
    samples_npz = directory / "c1g_2_effective_medium_samples.npz"
    np.savez_compressed(samples_npz, **{
        f"{case_id}_{name}": value
        for case_id, values in samples.items()
        for name, value in values.items()
    })
    medium_png = directory / "c1g_2_effective_medium.png"
    trajectory_png = directory / "c1g_2_effective_integral_trajectories.png"
    launch_png = directory / "c1g_2_launch_exclusion.png"
    _plot_effective_medium(
        rows, samples, solver_cases, exact_paths, context,
        (medium_png, trajectory_png, launch_png))
    images = {
        path.name: _image_check(path)
        for path in (medium_png, trajectory_png, launch_png)
    }
    high = [row for row in rows if row["branch"] == "high_ray"]
    low = [row for row in rows if row["branch"] == "low_ray"]
    launch_excluded = all(
        abs(row["nominal_minus_state_elevation_deg"]) <= 1e-10
        for row in rows)
    high_effective_closed = all(
        row["effective_trajectory_rms_km"] <= 0.15
        and abs(row["effective_range_error_km"]) <= 0.2
        for row in high)
    low_effective_closed = all(
        row["effective_trajectory_rms_km"] <= 0.15
        and abs(row["effective_range_error_km"]) <= 0.2
        for row in low)
    checks = {
        "four_pharlap_cases": len(rows) == 4,
        "nominal_and_state_launch_identical": launch_excluded,
        "both_high_close_with_state_effective_medium": high_effective_closed,
        "both_low_controls_close_with_state_effective_medium":
            low_effective_closed,
        "three_valid_images": all(value["passed"] for value in images.values()),
    }
    if all(checks.values()):
        status = "identified_c1g2_pharlap_effective_no_field_medium_candidate"
    else:
        status = "completed_c1g2_effective_medium_not_closed"
    summary = {
        "status": status,
        "created_utc": _utc(),
        "scope": "PHaRLAP_output_state_effective_n2_spherical_integral",
        "checks": checks,
        "records": rows,
        "image_checks": images,
        "interpretation": (
            "Launch angle is excluded. PHaRLAP preserves the spherical "
            "Bouguer invariant but follows an effective no-field n2 that "
            "differs from the shared frozen linear-grid n2."
            if all(checks.values()) else
            "State-derived effective medium did not uniquely close all cases."
        ),
    }
    summary_json = directory / "c1g_2_launch_effective_medium.json"
    _write_json(summary_json, summary)
    artifacts = [
        summary_json, metrics_csv, samples_npz,
        medium_png, trajectory_png, launch_png,
    ]
    _manifest(2, status, checks, artifacts, _sources(input_path),
              "state-based launch audit and PHaRLAP effective no-field n2; "
              "no explicit-state rerun; no WGS84/IGRF; no production changes")
    return status


def run_c1g6() -> str:
    stage1 = json.loads(
        (OUT / "c1g_1" / "c1g_1_spherical_benchmark.json").read_text(
            encoding="utf-8"))
    stage2 = json.loads(
        (OUT / "c1g_2" / "c1g_2_launch_effective_medium.json").read_text(
            encoding="utf-8"))
    if (stage1["status"]
            != "identified_c1g1_spherical_forward_equation_candidate"
            or stage2["status"]
            != "identified_c1g2_pharlap_effective_no_field_medium_candidate"):
        raise RuntimeError("C1G-1/2 did not satisfy the early-stop attribution gate")
    directory = OUT / "c1g_6"
    directory.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    stage0_manifest = json.loads(
        (OUT / "c1g_0" / "c1g_0_freeze_manifest.json").read_text(
            encoding="utf-8"))
    current_sources = _sources(input_path)
    production_keys = ("ah_solver", "ah_dispersion", "ah_provider", "config")
    checks = {
        "ah3d_unique_match_to_shared_spherical_integral":
            stage1["unique_solver_matching_linear_spherical_integral"] == "AH3D",
        "pharlap_effective_medium_closes_all_four_cases":
            stage2["checks"]["both_high_close_with_state_effective_medium"]
            and stage2["checks"]["both_low_controls_close_with_state_effective_medium"],
        "launch_geometry_excluded":
            stage2["checks"]["nominal_and_state_launch_identical"],
        "quadrature_uncertainty_passed":
            stage1["checks"]["quadrature_uncertainty_le_0_01_km"],
        "production_sources_unchanged": all(
            stage0_manifest["source_sha256"][key] == current_sources[key]
            for key in production_keys),
        "formal_c1_failure_retained": json.loads(
            FINAL_C1.read_text(encoding="utf-8"))["status"]
            == "failed_c1_chapman_high_low_external_comparison",
        "c1g3_to_c1g5_early_stop_allowed":
            "若C1G-1已由一维积分和守恒量唯一定位主因"
            in PLAN.read_text(encoding="utf-8"),
        "report_exists": REPORT.is_file(),
    }
    if not all(checks.values()):
        raise RuntimeError(f"C1G-6 final gate failed: {checks}")
    high_linear = [
        row for row in stage1["best_rows"] if row["branch"] == "high_ray"
    ]
    high_effective = [
        row for row in stage2["records"] if row["branch"] == "high_ray"
    ]
    closure = []
    for baseline, effective in zip(
            sorted((row for row in high_linear if row["solver"] == "PHaRLAP"),
                   key=lambda row: row["case_id"]),
            sorted(high_effective, key=lambda row: row["case_id"])):
        closure.append({
            "case_id": baseline["case_id"],
            "shared_input_integral_rms_km": baseline["trajectory_rms_km"],
            "effective_medium_integral_rms_km":
                effective["effective_trajectory_rms_km"],
            "trajectory_closure_percent": 100.0 * (
                1.0 - effective["effective_trajectory_rms_km"]
                / baseline["trajectory_rms_km"]),
            "shared_input_range_error_km": baseline["range_error_km"],
            "effective_medium_range_error_km":
                effective["effective_range_error_km"],
            "range_closure_percent": 100.0 * (
                1.0 - abs(effective["effective_range_error_km"])
                / abs(baseline["range_error_km"])),
        })
    result = {
        "status": "identified_spherical_equation_or_invariant_defect",
        "created_utc": _utc(),
        "affected_side": (
            "PHaRLAP_raytrace_3d_sp_OX_mode_0_effective_medium_or_dispersion_"
            "semantics; core source unavailable"),
        "ah3d_production_repair": "not_permitted_not_needed_by_c1g_evidence",
        "formal_c1_status": "failed_c1_chapman_high_low_external_comparison",
        "c2_gate": "manual_decision_required_c1_still_failed",
        "checks": checks,
        "closure": closure,
        "skipped_by_early_stop": {
            "C1G-3": "sphere_WGS84_differential_not_needed_for_B0_root_cause",
            "C1G-4": "B_to_zero_branch_scan_not_needed_for_shared_B0_input_attribution",
            "C1G-5": "local_angle_scan_not_needed_after_two_high_anchor_closure",
        },
        "interpretation": (
            "AH3D follows the shared spherical Chapman n2 integral. PHaRLAP "
            "preserves spherical invariants but follows a different effective "
            "no-field n2 reconstructed from its output state. The hidden MEX "
            "prevents source-level separation of interpolation and dispersion."
        ),
        "risks": [
            "PHaRLAP remains the project external standard, so formal C1 stays failed.",
            "The conclusion is limited to B=0 spherical horizontally uniform Chapman.",
            "PHaRLAP core source is unavailable; exact internal formula is unresolved.",
            "No evidence supports changing AH3D to reproduce the external effective n2.",
        ],
        "next_step": (
            "User decision: retain C1 failure and run A4-C2 only as "
            "diagnostic_with_c1_failed, or seek PHaRLAP vendor clarification."
        ),
    }
    result_json = directory / "c1g_6_final_attribution.json"
    _write_json(result_json, result)

    labels = list(checks)
    values = [1.0 if checks[label] else 0.0 for label in labels]
    figure, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].barh(range(len(labels)), values, color="#187b3f")
    axes[0].set_yticks(range(len(labels)),
                       [label.replace("_", " ") for label in labels])
    axes[0].set_xlim(0.0, 1.05)
    axes[0].set_xticks((0.0, 1.0), ("FAIL", "PASS"))
    axes[0].set_title("Final attribution gates")
    case_labels = [item["case_id"] for item in closure]
    x = np.arange(len(case_labels))
    width = 0.35
    axes[1].bar(
        x - width / 2,
        [item["shared_input_integral_rms_km"] for item in closure],
        width, label="PHaRLAP vs shared input integral", color="#c7362f")
    axes[1].bar(
        x + width / 2,
        [item["effective_medium_integral_rms_km"] for item in closure],
        width, label="PHaRLAP vs effective n2 integral", color="#1769aa")
    axes[1].axhline(0.15, color="#111111", ls="--", lw=1.0)
    axes[1].set_xticks(x, case_labels, rotation=15)
    axes[1].set(ylabel="Trajectory RMS (km)",
                title="High-ray attribution closure")
    axes[1].legend(fontsize=8)
    axes[1].grid(axis="y", alpha=0.25)
    figure.suptitle(
        "A4-C1G final | B=0 spherical Chapman | formal C1 failure retained")
    figure.tight_layout()
    final_png = directory / "c1g_6_final_attribution.png"
    figure.savefig(final_png, dpi=180)
    plt.close(figure)
    if not _image_check(final_png)["passed"]:
        raise RuntimeError("invalid C1G final visualization")
    artifacts = [result_json, final_png, REPORT]
    _manifest(6, result["status"], checks, artifacts, current_sources,
              "final attribution; C1G-3/4/5 skipped by frozen early-stop rule; "
              "formal C1 failure retained; no production repair")
    return result["status"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage", choices=("0", "1", "2", "6", "all"), default="all")
    args = parser.parse_args()
    if args.stage in {"0", "all"}:
        print(f"C1G-0: {run_c1g0()}")
    if args.stage in {"1", "all"}:
        print(f"C1G-1: {run_c1g1()}")
    if args.stage in {"2", "all"}:
        print(f"C1G-2: {run_c1g2()}")
    if args.stage in {"6", "all"}:
        print(f"C1G-6: {run_c1g6()}")
    return 0


context_global = None


if __name__ == "__main__":
    raise SystemExit(main())
