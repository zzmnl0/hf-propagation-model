"""Independent P2 evidence and numerical-gate verification.

This module reads archived native arrays directly. It does not import the
experiment runner or analysis module, and a verified archive is not a physics
acceptance. Historical A2 and P1 records remain immutable.
"""
from __future__ import annotations

import ast
import datetime as dt
import hashlib
import importlib.util
import json
import shutil
import tempfile
import traceback
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.integrate import quad, solve_ivp
from scipy.io import loadmat


CODE_ROOT = Path("D:/code11/IRI01/IRI03/TID/code")
A2_ROOT = CODE_ROOT / "output2/a4_pharlap/c1/c1m/a2_20260928_v1/a2"
P1_ROOT = CODE_ROOT / "output2/a4_pharlap/c1/c1m/p1_20261007_v1"
P1_DELIVERY_SHA = "d5bcb9578e57081c44c2dc3dbce735a015b5b39e442851b6d828c49ceecf2452"
A2_DELIVERY_SHA = "e85508c0c83e21df67de93b387905621214f5fb53682aba289164a0459d87a1f"
INTERNAL_GRID_SHA = "92f046c429c7365a06dbdf0572134762524cc24eb238d0dd1f69410bb1531d4e"
S_KM = 100.0
R_KM = 6371.0
C_KM_S = 299792.458
LIMITS = {"rms_m": 10.0, "max_m": 50.0, "landing_m": 50.0, "group_rel": 1e-5}
PH_PROFILES = ((1e-8, .01, 10.0), (1e-10, .001, 5.0), (1e-12, .001, 1.0))
METRICS = tuple(LIMITS)
P0_FROZEN_RUNNER_SHA = "973ec41896fe5828322a9cde2660f055c519c0d0a6ff67af9eb184806e4f99d1"
_P0_DECISION = None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _check(ok: bool, label: str, errors: list[str]) -> None:
    if not ok:
        errors.append(label)


def _native_vector(obj: object, field: str) -> np.ndarray:
    value = np.asarray(getattr(obj, field), dtype=float).reshape(-1)
    if value.size == 0 or not np.all(np.isfinite(value)):
        raise ValueError("invalid native field: " + field)
    return value


def _native_xyz(obj: object, names: tuple[str, str, str]) -> np.ndarray:
    return np.column_stack([_native_vector(obj, field) for field in names])


def _validate_trace(trace: dict, label: str) -> None:
    s, xyz, p = trace["s"], trace["xyz"], trace["p"]
    if (s.ndim != 1 or s.size < 2 or xyz.shape != (s.size, 3)
            or p.shape != xyz.shape or np.any(np.diff(s) <= 0.)
            or abs(float(s[0])) > 1e-12 or float(s[-1]) < S_KM
            or np.any(np.linalg.norm(p, axis=1) <= 0.)
            or any(np.any(~np.isfinite(v)) for v in (s, xyz, p))):
        raise ValueError("invalid fixed-100km native trace: " + label)
    if "group" in trace:
        if (trace["group"].shape != s.shape
                or np.any(np.diff(trace["group"]) <= 0.)
                or not np.all(np.isfinite(trace["group"]))):
            raise ValueError("invalid native group coordinate: " + label)
    elif (trace.get("termination") != "fixed_geometric_path_completed"
          or float(s[-1]) != S_KM or trace.get("group_total_km", 0.) <= 0.):
        raise ValueError("AH native group total is not the fixed 100km section: " + label)


def _ph(path: Path) -> tuple[dict, dict]:
    raw = loadmat(path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
    state, out = raw["ray_state_vec"], raw["ray_path_data"]
    s0 = _native_vector(state, "geometrical_path") / 1000.
    g0 = _native_vector(state, "group_path") / 1000.
    trace = {"s": s0 - s0[0], "xyz": _native_xyz(state, ("pos_x", "pos_y", "pos_z")) / 1000.,
             "p": _native_xyz(state, ("dir_x", "dir_y", "dir_z")), "group": g0 - g0[0],
             "s_absolute": s0, "group_absolute": g0,
             "path_s": _native_vector(out, "geometric_distance"),
             "path_group": _native_vector(out, "group_range"),
             "ne_out": _native_vector(out, "electron_density") * 1e6,
             "n_out": _native_vector(out, "refractive_index"),
             "ng_out": _native_vector(out, "group_refractive_index")}
    count = s0.size
    if any(v.shape[0] != count for v in trace.values()):
        raise ValueError("PH native path/state counts differ: " + path.name)
    for field in ("geomag_x", "geomag_y", "geomag_z", "collision_frequency"):
        if np.any(_native_vector(out, field) != 0.):
            raise ValueError("PH B0/collisionless contract violated: " + field)
    _validate_trace(trace, path.name)
    return trace, raw


def _ah(path: Path, frequency_MHz: float) -> tuple[dict, dict]:
    meta = _json(path)
    raw_path = path.with_suffix(".npz")
    if _sha(raw_path) != meta["raw_npz_sha256"]:
        raise ValueError("AH archived raw NPZ identity failed: " + path.name)
    with np.load(raw_path, allow_pickle=False) as raw:
        trace = {"s": np.asarray(raw["geom_path_km"], dtype=float),
                 "xyz": np.asarray(raw["xyz_km"], dtype=float),
                 "p": np.asarray(raw["wavevector_km_inv"], dtype=float)
                      * C_KM_S / (2*np.pi*frequency_MHz*1e6),
                 "exported_chord_tangent": np.asarray(raw["group_tangent_ecef"], dtype=float),
                 "ne_out": np.asarray(raw["Ne_m3_native"], dtype=float),
                 "n2_out": np.asarray(raw["n2_native"], dtype=float),
                 "group_total_km": float(meta["group_path_km_native"]),
                 "termination": meta["termination_status"]}
    _validate_trace(trace, path.name)
    return trace, meta


def _reconstruct(trace: dict, samples: np.ndarray, method: str) -> tuple[np.ndarray, float]:
    """p is c*k/omega; position derivatives are p/|p|, not p or grad(Ne)."""
    s, xyz, p = trace["s"], trace["xyz"], trace["p"]
    q = np.asarray(samples, dtype=float)
    if q[0] < s[0] or q[-1] > s[-1] or q[-1] != S_KM:
        raise ValueError("reconstruction must cover exactly the full 100km interval")
    if method == "pchip":
        positions = PchipInterpolator(s, xyz, axis=0, extrapolate=False)(q)
        group = (PchipInterpolator(s, trace["group"], extrapolate=False)(S_KM)
                 if "group" in trace else trace["group_total_km"])
    elif method in ("direction_hermite", "derivative_hermite"):
        norm = np.linalg.norm(p, axis=1)
        positions = CubicHermiteSpline(s, xyz, p/norm[:, None], axis=0,
                                      extrapolate=False)(q)
        group = (CubicHermiteSpline(s, trace["group"], 1./norm,
                                   extrapolate=False)(S_KM)
                 if "group" in trace else trace["group_total_km"])
    else:
        raise ValueError("unknown prescribed reconstruction: " + method)
    if np.any(~np.isfinite(positions)) or not np.isfinite(group):
        raise ValueError("nonfinite reconstructed fixed-100km values")
    return positions, float(group)


def _pair(a: dict, b: dict, method: str, step_km: float) -> dict:
    count = int(round(S_KM/step_km))
    if count < 1 or count*step_km != S_KM:
        raise ValueError("sampling step must partition the fixed100km interval")
    q = np.linspace(0., S_KM, count+1)
    xa, ga = _reconstruct(a, q, method)
    xb, gb = _reconstruct(b, q, method)
    error_m = np.linalg.norm(xa-xb, axis=1)*1000.
    return {"rms_m": float(np.sqrt(np.mean(error_m**2))), "max_m": float(error_m.max()),
            "landing_m": float(error_m[-1]), "group_rel": abs(ga-gb)/abs(ga)}


def _polynomial_chord_upper(s: np.ndarray, values: np.ndarray, spline) -> tuple[float, float]:
    """Bound a prescribed spline's distance from each native endpoint chord.

    Component extrema include all real derivative roots in each interval.
    Combining component maxima with the Euclidean norm remains an upper bound.
    The physical curve-to-chord term is checked separately.
    """
    vectors = np.asarray(values, dtype=float)
    if vectors.ndim == 1:
        vectors = vectors[:, None]
    coefficients = spline.c
    if coefficients.ndim == 2:
        coefficients = coefficients[:, :, None]
    deviation, gap = 0., 0.
    for index in range(len(s)-1):
        if s[index] >= S_KM:
            break
        width = float(s[index+1]-s[index])
        end = min(width, S_KM-float(s[index]))
        slope = (vectors[index+1]-vectors[index])/width
        maxima = []
        for component in range(vectors.shape[1]):
            c = coefficients[:, index, component].copy()
            c[2] -= slope[component]
            c[3] -= vectors[index, component]
            extrema = [0., end]
            for root in np.roots(np.array([3*c[0], 2*c[1], c[2]])):
                if abs(float(np.imag(root))) < 1e-12 and 0. < float(np.real(root)) < end:
                    extrema.append(float(np.real(root)))
            maxima.append(max(abs(float(np.polyval(c, point))) for point in extrema))
        deviation = max(deviation, float(np.linalg.norm(maxima)))
        gap = max(gap, width)
    return float(np.nextafter(deviation, np.inf)), gap


def _reconstruction_upper(trace: dict, method: str, acceleration_upper: float,
                          group_second_derivative_upper: float) -> dict:
    """A deterministic native interpolation bound given verified derivative bounds.

    A C2 curve with ||x''||<=M differs from its endpoint chord by <=M*h^2/8.
    This result requires independent field/trajectory support for M and for any
    native-node numerical error. Analytic implementation controls alone do not
    establish those prerequisites.
    """
    if (not np.isfinite(acceleration_upper) or acceleration_upper < 0.
            or not np.isfinite(group_second_derivative_upper)
            or group_second_derivative_upper < 0.):
        raise ValueError("invalid reconstruction derivative bounds")
    s, xyz, norm = trace["s"], trace["xyz"], np.linalg.norm(trace["p"], axis=1)
    if method == "pchip":
        curve = PchipInterpolator(s, xyz, axis=0, extrapolate=False)
    elif method in ("direction_hermite", "derivative_hermite"):
        curve = CubicHermiteSpline(s, xyz, trace["p"]/norm[:, None], axis=0, extrapolate=False)
    else:
        raise ValueError("unknown reconstruction upper-bound method")
    polynomial_deviation_km, gap_km = _polynomial_chord_upper(s, xyz, curve)
    position_chord_km = acceleration_upper*gap_km**2/8.
    output = {"position_upper_m": 1000.*(position_chord_km+polynomial_deviation_km),
              "position_chord_upper_m": 1000.*position_chord_km,
              "position_polynomial_chord_upper_m": 1000.*polynomial_deviation_km,
              "max_native_gap_km": gap_km,
              "group_upper_km": 0., "group_polynomial_chord_upper_km": 0.,
              "native_node_error_included": False}
    if "group" in trace:
        group_curve = (PchipInterpolator(s, trace["group"], extrapolate=False) if method == "pchip"
                       else CubicHermiteSpline(s, trace["group"], 1./norm, extrapolate=False))
        group_deviation, group_gap = _polynomial_chord_upper(s, trace["group"], group_curve)
        output["group_polynomial_chord_upper_km"] = group_deviation
        output["group_upper_km"] = group_second_derivative_upper*group_gap**2/8.+group_deviation
    return output


def _independent_decision(actual: dict, d12: dict, d23: dict, reconstruction: dict,
                          eligible: bool) -> dict:
    """No roundoff waiver and no interpolation-choice spread as a budget."""
    checks = {}
    for metric, limit in LIMITS.items():
        values = [actual[metric], reconstruction[metric]]
        values += [ref[side][metric] for ref in (d12, d23) for side in ("ph", "ah")]
        if any(not isinstance(v, (int, float)) or not np.isfinite(v) or v < 0. for v in values):
            raise ValueError("invalid numerical-gate evidence: " + metric)
        sides = {side: {"monotonic": bool(d23[side][metric] <= d12[side][metric]),
                        "last_change_small": bool(d23[side][metric] <= .2*limit)}
                 for side in ("ph", "ah")}
        budget = 3.*(d23["ph"][metric]+d23["ah"][metric]+reconstruction[metric])
        checks[metric] = {"physical_pass": bool(actual[metric] <= limit),
                          "explained_by_numeric": bool(actual[metric] <= budget),
                          "three_times_budget": budget, "sides": sides}
    passed = bool(eligible) and all(row["physical_pass"] and row["explained_by_numeric"]
                                   and all(all(s.values()) for s in row["sides"].values())
                                   for row in checks.values())
    return {"passed": passed, "eligible": bool(eligible), "checks": checks,
            "roundoff_waiver_applied": False, "interpolation_method_difference_in_budget": False}


def _decision(actual: dict, d12: dict, d23: dict, reconstruction: dict,
              eligible: bool) -> dict:
    """Reuse the exact P0 gate and cross-check its result independently."""
    global _P0_DECISION
    if _P0_DECISION is None:
        source = P1_ROOT / "sources/runner.py"
        if _sha(source) != P0_FROZEN_RUNNER_SHA:
            raise ValueError("P0 frozen rule source identity failed")
        tree = ast.parse(source.read_text(encoding="utf-8-sig"))
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                        and node.name == "numeric_decision_v2")
        namespace = {"np": np, "LIMIT": {k: v for k, v in LIMITS.items() if k != "landing_m"}}
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
        _P0_DECISION = namespace["numeric_decision_v2"]
    expected = _independent_decision(actual, d12, d23, reconstruction, eligible)
    original = _P0_DECISION(actual, d12, d23, reconstruction, eligible=eligible)
    if original != expected:
        raise ValueError("independent numerical gate differs from frozen P0 rules")
    return original


def _initial_matches(xyz0: np.ndarray, p0: np.ndarray, trace: dict,
                     exact: bool = False) -> bool:
    if exact:
        return bool(np.array_equal(xyz0, trace["xyz"][0]) and np.array_equal(p0, trace["p"][0]))
    return bool(np.max(abs(trace["xyz"][0]-xyz0)) < 1e-10
                and np.max(abs(trace["p"][0]-p0)) < 1e-14)


def _verify_delivery(root: Path, filename: str, expected_sha: str,
                     expected_count: int, errors: list[str]) -> dict:
    manifest = root / filename
    _check(_sha(manifest) == expected_sha, filename+": manifest identity", errors)
    records = _json(manifest)["files" if filename == "delivery_sha256.json" else "records"]
    _check(len(records) == expected_count, filename+": record count", errors)
    n_ok, n_bytes = 0, 0
    for record in records:
        path = root / record["path"]
        try:
            path.resolve().relative_to(root.resolve())
        except ValueError:
            errors.append(filename+": out-of-scope record "+record["path"])
            continue
        good = path.is_file() and path.stat().st_size == record["bytes"] and _sha(path) == record["sha256"]
        _check(good, filename+": archived file identity "+record["path"], errors)
        n_ok += int(good)
        n_bytes += record["bytes"] if good else 0
    return {"manifest_sha256": _sha(manifest), "records": len(records), "verified": n_ok,
            "verified_bytes": n_bytes, "historical_external_current_tool_sha_compared": False}


def _verify_legacy(errors: list[str]) -> tuple[dict, dict]:
    archive_manifest = A2_ROOT / "sha256_manifest.json"
    _check(_sha(archive_manifest) == A2_DELIVERY_SHA, "A2 historical manifest identity", errors)
    historical = {str(Path(r["path"]).resolve()): r for r in _json(archive_manifest)["records"]}
    records, traces = [], {}

    def historical_file(path: Path) -> str:
        key = str(path.resolve())
        if key not in historical:
            raise ValueError("original A2 file missing from historical manifest: "+key)
        record = historical[key]
        digest = _sha(path)
        _check(digest == record["sha256"] and path.stat().st_size == record["bytes"],
               "A2 historical file identity: "+path.name, errors)
        return digest

    base = _json(A2_ROOT / "input_manifest.json")
    grid_path = Path(base["grid"]["internal_common"]["grid_mat"])
    grid_sha = historical_file(grid_path)
    _check(grid_sha == INTERNAL_GRID_SHA, "internal profile SHA", errors)
    grid = loadmat(grid_path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
    h, ne = np.asarray(grid["height_km"]), np.asarray(grid["ne_profile_m3"])
    _check(np.array_equal(h, np.arange(0., 462., 2.)), "internal grid height nodes", errors)
    _check(np.array_equal(grid["Ne_m3"], np.broadcast_to(ne, grid["Ne_m3"].shape)),
           "internal radial density array", errors)
    for field in ("Bx_T", "By_T", "Bz_T", "collision_freq_Hz"):
        _check(np.all(grid[field] == 0.), "internal grid B0/collisionless "+field, errors)
    _check(float(grid["earth_radius_m"]) == R_KM*1000., "internal sphere radius", errors)
    for variant in ("common", "identified"):
        directory = A2_ROOT if variant == "common" else A2_ROOT / "identified"
        manifest_path = directory / "input_manifest.json"
        historical_file(manifest_path)
        manifest = _json(manifest_path)
        input_sha = historical_file(directory / "input.mat")
        _check(input_sha == manifest["input_mat_sha256"], variant+": input identity", errors)
        jobs = [j for j in manifest["jobs"] if j["stage"] == "internal"]
        _check(len(jobs) == 6, variant+": six original internal jobs", errors)
        _check({j["elevation_deg"] for j in jobs} == {22., 39.5}, variant+": original angles", errors)
        for job in jobs:
            label = job["job_id"]
            j0 = {k: v for k, v in job.items() if k != "job_sha256"}
            digest = hashlib.sha256(json.dumps(j0, sort_keys=True, separators=(",", ":"),
                                               allow_nan=False).encode()).hexdigest()
            _check(digest == job["job_sha256"], label+": canonical job identity", errors)
            _check(tuple(job["tol"]) == PH_PROFILES[job["level"]-1], label+": legal PH profile", errors)
            _check(job["frequency_MHz"] == 10. and job["origin_height_km"] == 100.
                   and job["explicit_state"] is True, label+": original launch contract", errors)
            path = directory / "pharlap_raw" / (label+".mat")
            ph_sha = historical_file(path)
            ph, raw = _ph(path)
            _check(str(raw["job_sha256"]) == job["job_sha256"]
                   and str(raw["grid_sha256"]) == INTERNAL_GRID_SHA, label+": PH embedded identities", errors)
            initial = raw["initial"]
            p0 = np.array([float(getattr(initial, field)) for field in ("dir_x", "dir_y", "dir_z")])
            xyz0 = np.array([float(getattr(initial, field)) for field in ("pos_x", "pos_y", "pos_z")]) / 1000.
            _check(_initial_matches(xyz0, p0, ph, exact=True),
                   label+": PH initial/native state", errors)
            coefficient = 80.6404 if variant == "common" else 80.61638494312905
            _check(job["launch_k2"] == coefficient, label+": variant coefficient", errors)
            n2_initial = 1.-coefficient*float(ne[50])/(10e6)**2
            _check(abs(float(p0@p0)-n2_initial) < 2e-15, label+": initial dispersion", errors)
            ah_path = directory / "ah/formal" / variant / "hermite" / (label+".json")
            ah_json_sha = historical_file(ah_path)
            ah_raw_sha = historical_file(ah_path.with_suffix(".npz"))
            ah, meta = _ah(ah_path, job["frequency_MHz"])
            ah_coefficient = 80.64040000000001 if variant == "common" else coefficient
            _check(meta["job_sha256"] == job["job_sha256"] and meta["grid_sha256"] == INTERNAL_GRID_SHA
                   and meta["coefficient_K2"] == ah_coefficient and meta["interpolation"] == "hermite",
                   label+": AH native identity/medium", errors)
            _check(_initial_matches(xyz0, p0, ah),
                   label+": matched AH/PH initial state", errors)
            traces[(variant, float(job["elevation_deg"]), int(job["level"]))] = {"ph": ph, "ah": ah}
            records.append({"job_id": label, "variant": variant, "level": int(job["level"]),
                            "angle_deg": float(job["elevation_deg"]), "coefficient_K2": coefficient,
                            "ph_sha256": ph_sha, "ah_json_sha256": ah_json_sha, "ah_npz_sha256": ah_raw_sha,
                            "ph_native_end_km": float(ph["s"][-1]), "ah_native_end_km": float(ah["s"][-1]),
                            "ph_initial_p2": float(p0@p0), "n2_initial_from_input": n2_initial})
    baseline_report = A2_ROOT / "identified/analysis/a2_identified_formal_identified_hermite_internal.json"
    baseline_sha = historical_file(baseline_report)
    baseline_rows = []
    for row in _json(baseline_report)["rows"]:
        if row["level"] != 3:
            continue
        angle = float(row["elevation_deg"])
        actual = _pair(traces[("identified", angle, 3)]["ph"], traces[("identified", angle, 3)]["ah"], "pchip", .5)
        reported = {"rms_m": row["position_rms_m"], "max_m": row["position_max_m"],
                    "landing_m": row["position_endpoint_m"], "group_rel": row["group_relative_error"]}
        for metric in METRICS:
            _check(actual[metric] == reported[metric], "old A2 baseline exact reproduction "+metric, errors)
        baseline_rows.append({"angle_deg": angle, "method": "pchip", "query_step_km": .5,
                              "recomputed": actual, "old_report": reported, "max_difference": max(abs(actual[k]-reported[k]) for k in METRICS)})
    return {"archive_manifest_sha256": _sha(archive_manifest), "grid_sha256": grid_sha,
            "original_records": records, "frozen100km": S_KM,
            "old_A2_report_sha256": baseline_sha, "old_A2_baseline_reproductions": baseline_rows,
            "historical_external_current_tool_sha_compared": False}, traces


def _verify_p2_sources(run: Path, errors: list[str], analysis_root: Path) -> dict:
    """Each executed component is bound before its own execution, by snapshot."""
    prepare_path = run / "p2_manifest.json"
    analysis_path = analysis_root / "analysis_manifest.json"
    prepare, analysis = _json(prepare_path), _json(analysis_path)
    _check(prepare["schema"] == "a4-c1m-p2-prepare-v1", "P2 prepare schema", errors)
    _check(analysis["schema"] in ("a4-c1m-p2-analysis-source-v1", "a4-c1m-p2-analysis-source-v2"),
           "P2 analysis source schema", errors)
    _check(analysis["prepared_manifest_sha256"] == _sha(prepare_path), "analysis->prepare manifest link", errors)
    _check(prepare["arc_km"] == S_KM and prepare["limits"] == LIMITS, "P2 frozen interval/thresholds", errors)
    _check(tuple(tuple(v) for v in prepare["PH_profiles"]) == PH_PROFILES, "P2 frozen three legal PH profiles", errors)
    _check(prepare["original_fd_km"] == .5 and prepare["candidate_fd_km"] == .125,
           "P2 baseline/candidate FD", errors)
    _check(prepare["constants"] == {"C": 80.61638494312905,
                                   "coefficient_source": str(A2_ROOT / "identified/coefficient_freeze.json"),
                                   "R_km": R_KM, "C_km_s": C_KM_S}, "P2 frozen physical constants", errors)
    _check(prepare["auxiliary_common_is_formal"] is False and prepare["original_A2_provenance_complete"] is False,
           "P2 auxiliary/original provenance scope", errors)
    _check(prepare["analysis_implementation"]["prepared_at_this_stage"] is False,
           "P2 records later-bound analysis source honestly", errors)
    source_rows = []
    for stage, block, directory in (("prepare", prepare["sources"], run),
                                    ("analysis", analysis["sources"], analysis_path.parent)):
        for name, record in block.items():
            if record["snapshot"]:
                path = directory / record["snapshot"]
                path.resolve().relative_to(directory.resolve())
                good = path.is_file() and _sha(path) == record["sha256"]
                checked = "frozen_snapshot"
            else:
                path = Path(record["path"])
                good = path.is_file() and _sha(path) == record["sha256"]
                checked = "original_runtime_binary"
            _check(good, stage+": source identity "+name, errors)
            source_rows.append({"stage": stage, "source": name, "checked": checked,
                                "sha256": record["sha256"], "verified": good})
    _check({"runner", "p2_prepare_dispatcher", "ah_solver", "ah_dispersion", "ah_provider"}
           <= set(prepare["sources"]), "P2 prepare source coverage", errors)
    _check("core" in analysis["sources"], "P2 analysis core source coverage", errors)
    for label in ("production_modified", "fitted_parameters", "reconstruction_method_difference_in_budget"):
        _check(analysis[label] is False, "P2 analysis contract: "+label, errors)
    if analysis["schema"].endswith("-v2"):
        for name in ("previous_analysis_manifest", "previous_artifact_manifest"):
            record = analysis[name]
            _check(_sha(Path(record["path"])) == record["sha256"], "P2 v2 prior provenance link: "+name, errors)
        for record in analysis["historical_v1_records"]:
            _check(_sha(Path(record["path"])) == record["sha256"],
                   "P2 v2 historical v1 artifact retained: "+record["path"], errors)
        _check(len(analysis["reused_AH_candidates"]) == 6, "P2 v2 reuses exactly six candidate AH records", errors)
        for record in analysis["reused_AH_candidates"]:
            for name in ("meta", "raw"):
                _check(_sha(Path(record[name])) == record[name+"_sha256"], "P2 v2 candidate reuse identity", errors)
        _check(analysis["quad_extra_levels"] == [4, 5] and analysis["quad_extra_angle_deg"] == 39.5
               and analysis["maximum_additional_quad_profiles"] == 2 and analysis["uncertainty_multiplier"] == 1.
               and analysis["numerical_estimates_are_not_interval_certified"] is True
               and analysis["all_v1_failures_preserved"] is True and analysis["PH_recomputed"] is False
               and analysis["AH_recomputed"] is False, "P2 v2 supplement frozen scope", errors)
    input_count = 0
    for filename, record in prepare["historical_inputs"].items():
        path = Path(filename)
        good = path.is_file() and _sha(path) == record["sha256"]
        if "bytes" in record:
            good = good and path.stat().st_size == record["bytes"]
        _check(good, "P2 historical input identity: "+filename, errors)
        input_count += 1
    for filename, digest in prepare["prepared_inputs"].items():
        path = Path(filename)
        path.resolve().relative_to(run.resolve())
        _check(path.is_file() and _sha(path) == digest, "P2 prepared input identity: "+filename, errors)
    exit_ref = prepare["exit_component"]
    _check(_sha(Path(exit_ref["manifest"])) == exit_ref["sha256"], "P2 exit manifest reference", errors)
    _check(_sha(run / "inputs/internal_common.mat") == INTERNAL_GRID_SHA, "P2 internal grid copy", errors)
    selected = _json(run / "inputs/internal_jobs.json")
    _check(len(selected["primary"]) == 6 and len(selected["auxiliary"]) == 6,
           "P2 selected original primary/auxiliary count", errors)
    for formal, items in ((True, selected["primary"]), (False, selected["auxiliary"])):
        expected_variant = "identified" if formal else "common"
        for item in items:
            _check(item["formal"] is formal and item["common_arc_km"] == S_KM
                   and item["original_fd_km"] == .5, "P2 selected job contract", errors)
            job = item["original_job"]
            path = Path(item["PH_raw"])
            expected_folder = A2_ROOT / "identified" if formal else A2_ROOT
            _check(path == expected_folder / "pharlap_raw" / (job["job_id"]+".mat"),
                   "P2 original PH path binding", errors)
            _check(Path(item["AH_meta"]) == expected_folder / "ah/formal" / expected_variant
                   / "hermite" / (job["job_id"]+".json"), "P2 original matched AH path binding", errors)
            _check(item["input_Ne100_m3"] == 7911590468.012695,
                   "P2 selected density from input grid", errors)
    return {"prepare_manifest_sha256": _sha(prepare_path), "analysis_manifest_sha256": _sha(analysis_path),
            "source_records": source_rows, "historical_inputs_verified_count": input_count,
            "later_components_bound_at_execution": True,
            "live_earlier_tool_changes_treated_as_archive_tamper": False}


def _verify_production(run: Path, errors: list[str]) -> dict:
    original = _json(P1_ROOT / "p1_manifest.json")
    prepare = _json(run / "p2_manifest.json")
    records = []
    for key in ("ah_solver", "ah_dispersion", "ah_provider", "config", "utils", "pharlap_mex"):
        old, new = original["sources"][key], prepare["sources"][key]
        path = Path(new["path"])
        good = old["sha256"] == new["sha256"] and path.is_file() and _sha(path) == old["sha256"]
        _check(good, "production/MEX unchanged from P1: "+key, errors)
        records.append({"source": key, "path": str(path), "P1_sha256": old["sha256"],
                        "P2_sha256": new["sha256"], "current_sha256": _sha(path), "unchanged": good})
    return {"records": records, "production_and_mex_unchanged": all(r["unchanged"] for r in records),
            "A2_old_preparation_provenance_gap_resolved": False,
            "current_verification_does_not_establish_historical_preexecution_freeze": True}


def _verify_exit(run: Path, errors: list[str]) -> dict:
    """Verify every exit restart node, rather than extrapolating a first step."""
    root = run / "exit"
    retry = root / "retry_nomem_v1"
    summary_path = retry / "analysis_v3/summary.json"
    summary = _json(summary_path)
    _check(summary["diagnostic_execution_complete"] is True, "exit diagnostic completion", errors)
    _check(all(v is True for v in summary["checks"].values()), "exit v3 all checks", errors)
    for field in ("physical_acceptance_passed", "global_passed", "closed_core_location_identified",
                  "independent_ODE_modified_to_fit_PH"):
        _check(summary[field] is False, "exit limited claim: "+field, errors)
    _check(summary["earlier_outputs_preserved"] is True and bool(summary["supersedes"]),
           "exit v3 preserves/supersedes early interpretations", errors)
    for older in (root / "analysis_v1/summary.json", retry / "analysis_v2/summary.json"):
        _check(older.is_file(), "exit earlier analysis preserved "+str(older), errors)
    manifest = _json(root / "manifest.json")
    _check(manifest["source_freeze_precedes_input_generation"] is True,
           "exit prepare source freeze precedes input generation", errors)
    for record in manifest["sources"]:
        snapshot = root / record["snapshot"]
        snapshot.resolve().relative_to(root.resolve())
        _check(snapshot.is_file() and _sha(snapshot) == record["sha256"],
               "exit frozen preparation source "+record["path"], errors)
    for field in ("inputs_sha256", "external_sha256"):
        for name, digest in manifest[field].items():
            path = Path(name)
            _check(path.is_file() and _sha(path) == digest, "exit "+field+": "+name, errors)
    retry_manifest = _json(retry / "manifest.json")
    _check(retry_manifest["source_freeze_precedes_input_generation"] is True
           and retry_manifest["parameter_changed"] is False
           and retry_manifest["original_failed_execution_preserved"] is True
           and retry_manifest["retry_count"] == 1, "exit exact retry authorization/scope", errors)
    for field in ("source_sha256", "input_sha256"):
        for name, digest in retry_manifest[field].items():
            path = Path(name)
            if not path.is_absolute():
                path = retry / path
            _check(path.is_file() and _sha(path) == digest, "exit retry "+field+": "+name, errors)
    v2_source = _json(retry / "execution_analysis_source_freeze_v2.json")
    _check(_sha(Path(v2_source["source"])) == v2_source["sha256"]
           and v2_source["preparation_source_v1_preserved"] is True,
           "exit retry v2 executed analysis source", errors)
    sequence_manifest = _json(retry / "sequence_manifest_v3.json")
    _check(sequence_manifest["preparation_sources_frozen_before_input_hashes"] is True,
           "exit v3 sequence source freeze", errors)
    for field in ("sources_sha256", "inputs_sha256"):
        for name, digest in sequence_manifest[field].items():
            path = Path(name)
            _check(path.is_file() and _sha(path) == digest, "exit sequence v3 "+field+": "+name, errors)
    jobs = _json(root / "jobs.json")["jobs"]
    retry_jobs = _json(retry / "jobs.json")["jobs"]
    _check(len(jobs) == 12 and len(retry_jobs) == 1 and retry_jobs[0] == jobs[0],
           "exit retry preserves the exact original job", errors)
    raw_records = _json(retry / "analysis_v2/raw_sources.json")["records"]
    raw_map = {r["job_id"]: r for r in raw_records}
    _check(len(raw_records) == len(raw_map) == 12 and set(raw_map) == {j["job_id"] for j in jobs},
           "exit complete merged raw-source mapping", errors)
    original_execution = _json(root / "execution.json")
    retry_execution = _json(retry / "execution.json")
    _check(summary["original_batch_returncode"] == original_execution["returncode"]
           and summary["retry_returncode"] == retry_execution["returncode"],
           "exit preserved execution attempt returncodes", errors)
    rows, natural = [], {}
    for job in jobs:
        label = job["job_id"]
        data = {k: v for k, v in job.items() if k != "job_sha256"}
        canonical = hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                              allow_nan=False).encode()).hexdigest()
        _check(canonical == job["job_sha256"], label+": canonical job identity", errors)
        _check(tuple(job["tol"]) == PH_PROFILES[job["level"]-1]
               and job["frequency_MHz"] == 10. and job["diagnostic_only"] is True,
               label+": legal diagnostic profile", errors)
        source = Path(job["source_raw_mat"])
        _check(_sha(source) == job["source_raw_sha256"], label+": original natural PH identity", errors)
        raw_record = raw_map[label]
        path = Path(raw_record["raw_mat"])
        expected_root = retry if label == jobs[0]["job_id"] else root
        _check(path == expected_root / "pharlap_raw" / (label+".mat"), label+": merged attempt source", errors)
        _check(_sha(path) == raw_record["sha256"], label+": restart raw identity", errors)
        raw = loadmat(path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
        for field in ("job_sha256", "grid_sha256", "exporter_sha256", "mex_sha256"):
            _check(str(raw[field]) == job[field], label+": embedded "+field, errors)
        for field, value in job["initial"].items():
            native_field = "geometrical_path" if field == "geom_path" else field
            _check(float(getattr(raw["initial"], field)) == value
                   and float(_native_vector(raw["ray_state_vec"], native_field)[0]) == value,
                   label+": inherited initial field "+field, errors)
        state, out = raw["ray_state_vec"], raw["ray_path_data"]
        p = _native_xyz(state, ("dir_x", "dir_y", "dir_z"))
        p2 = np.sum(p*p, axis=1)
        s, group = _native_vector(state, "geometrical_path"), _native_vector(state, "group_path")
        n, ng = _native_vector(out, "refractive_index"), _native_vector(out, "group_refractive_index")
        ne = _native_vector(out, "electron_density")*1e6
        h = _native_vector(out, "height")
        _check(len(s) >= 5 and np.all(np.diff(s) > 0.) and np.all(np.diff(group) > 0.),
               label+": complete monotonically advancing native sequence", errors)
        _check(np.array_equal(group, _native_vector(state, "indep_var")),
               label+": all-node group/independent coordinate equality", errors)
        path_s, path_group = _native_vector(out, "geometric_distance"), _native_vector(out, "group_range")
        _check(np.max(abs(s-path_s*1000.)) < 5e-9 and np.max(abs(group-path_group*1000.)) < 5e-9,
               label+": path/state all-node coordinate units", errors)
        slope = np.diff(group)/np.diff(s)
        role = job["position_role"]
        if role in ("first_outside", "next_outside"):
            _check(np.all(ne == 0.) and np.all(h < job["h_boundary_km"]),
                   label+": full exterior vacuum sequence", errors)
            _check(abs(float(p2[0]-p2[1])) < 1e-14 and abs(float(p2[0]-1.)) > 1e-6
                   and np.max(abs(p2[2:]-1.)) < 1e-13,
                   label+": p magnitude changes at node2 and stays unit", errors)
            _check(abs(float(n[0]-np.sqrt(p2[0]))) < 1e-13
                   and abs(float(ng[0]-1./np.sqrt(p2[0]))) < 1e-13
                   and np.max(abs(n[1:]-1.)) < 1e-13 and np.max(abs(ng[1:]-1.)) < 1e-13,
                   label+": path n/ng changes at node1 and stays vacuum", errors)
            _check(abs(float(slope[0]-1./np.sqrt(p2[0]))) < 1e-11
                   and np.max(abs(slope[1:]-1.)) < 1e-11,
                   label+": first G slope inherited then all later slopes unit", errors)
        elif role == "first_outside_unit":
            _check(np.all(ne == 0.) and np.max(abs(p2-1.)) < 1e-13
                   and np.max(abs(n-1.)) < 1e-13 and np.max(abs(ng-1.)) < 1e-13
                   and np.max(abs(slope-1.)) < 1e-11,
                   label+": normalized restart stays unit/vacuum for all nodes", errors)
        else:
            _check(role == "last_inside", label+": expected exit role", errors)
        rows.append({"job_id": label, "role": role, "level": job["level"],
                     "native_nodes": len(s), "raw_sha256": raw_record["sha256"],
                     "attempt": raw_record["execution_attempt"], "p2_all_nodes": p2.tolist(),
                     "n_all_nodes": n.tolist(), "ng_all_nodes": ng.tolist(),
                     "group_slope_all_intervals": slope.tolist(),
                     "height_km_all_nodes": h.tolist(), "Ne_m3_all_nodes": ne.tolist()})
        if job["level"] not in natural:
            original = loadmat(source, squeeze_me=True, struct_as_record=False, mat_dtype=True)
            st, op = original["ray_state_vec"], original["ray_path_data"]
            ns = _native_vector(st, "geometrical_path")
            np2 = np.sum(_native_xyz(st, ("dir_x", "dir_y", "dir_z"))**2, axis=1)
            nh = _native_vector(op, "height")
            apex = int(np.argmax(nh))
            exterior = np.flatnonzero(nh[apex:] < job["h_boundary_km"])+apex
            crossing = int(exterior[0])
            _check(crossing > 0 and crossing+2 < len(nh), "natural complete exit sequence", errors)
            nn, nng = _native_vector(op, "refractive_index"), _native_vector(op, "group_refractive_index")
            nG = _native_vector(st, "group_path")
            _check(abs(float(np2[crossing]-1.)) > 1e-6 and np.max(abs(np2[crossing+1:]-1.)) < 1e-13,
                   "natural exterior p unit from subsequent native node", errors)
            _check(np.max(abs(nn[crossing:]-1.)) < 1e-13 and np.max(abs(nng[crossing:]-1.)) < 1e-13,
                   "natural first and subsequent exterior path n/ng unit", errors)
            _check(np.max(abs(np.diff(nG)[crossing:]/np.diff(ns)[crossing:]-1.)) < 1e-11,
                   "natural all exterior G slopes unit", errors)
            start = crossing-1
            natural[job["level"]] = {"source_raw_sha256": job["source_raw_sha256"],
                                      "first_exterior_index": crossing,
                                      "p2_last_inside_through_end": np2[start:].tolist(),
                                      "n_last_inside_through_end": nn[start:].tolist(),
                                      "ng_last_inside_through_end": nng[start:].tolist(),
                                      "Gprime_last_inside_through_end": (np.diff(nG)[start:]/np.diff(ns)[start:]).tolist(),
                                      "crossing_interval_G_semantics_resolved": False}
    return {"latest_summary_sha256": _sha(summary_path), "original_attempt_returncode": summary["original_batch_returncode"],
            "retry_returncode": summary["retry_returncode"], "restart_sequences": rows,
            "natural_sequences": natural, "all_nodes_verified": True,
            "physical_acceptance_passed": False, "closed_core_location_identified": False,
            "earlier_interpretations_superseded_and_preserved": True}


def _sample_native(trace: dict, stride: int) -> dict:
    """Preserve PH brackets and the next real support node for PCHIP slopes."""
    s = trace["s"]
    right = int(np.searchsorted(s, S_KM, side="left"))
    indices = set(range(0, right+1, stride)) | {0, right}
    if right > 0:
        indices.add(right-1)
    if right+1 < len(s):
        indices.add(right+1)
        indices.add(right-1)
    selected = np.array(sorted(indices), dtype=int)
    output = {}
    for key, value in trace.items():
        output[key] = value[selected] if isinstance(value, np.ndarray) and value.shape[0] == s.size else value
    _validate_trace(output, "native node stride"+str(stride))
    return output


def _near_metrics(given: dict, computed: dict, label: str, errors: list[str]) -> None:
    _check(set(given) == set(METRICS), label+": four fixed metrics", errors)
    for key in METRICS:
        # This checks serialized evidence against independently recomputed data.
        # It is not a relaxation of any scientific gate or convergence comparison.
        absolute = 1e-14 if key == "group_rel" else 1e-7
        value = given.get(key, float("nan"))
        _check(np.isfinite(value) and np.isclose(value, computed[key], rtol=2e-9, atol=absolute),
               label+": recomputed "+key, errors)


def _verify_AH_derivative_semantics(run: Path, original: dict, errors: list[str]) -> dict:
    rows = []
    for angle in (22., 39.5):
        for fd in (.5, .125):
            for level in (1, 2, 3):
                tag = "22" if angle == 22. else "39p5"
                name = f"internal_identified_f10_e{tag}_t{level}.json"
                path = (A2_ROOT / "identified/ah/formal/identified/hermite" / name if fd == .5
                        else run / "ah/fd0.125" / name)
                trace, meta = _ah(path, 10.)
                xyz, p, s = trace["xyz"], trace["p"], trace["s"]
                chord = np.empty_like(xyz)
                chord[0] = xyz[1]-xyz[0]
                chord[-1] = xyz[-1]-xyz[-2]
                chord[1:-1] = xyz[2:]-xyz[:-2]
                chord /= np.linalg.norm(chord, axis=1)[:, None]
                exported = trace["exported_chord_tangent"]
                _check(np.max(abs(chord-exported)) < 1e-13,
                       name+": exported group tangent is a normalized native chord", errors)
                physical = p/np.linalg.norm(p, axis=1)[:, None]
                angles = np.arctan2(np.linalg.norm(np.cross(exported, physical), axis=1),
                                    np.sum(exported*physical, axis=1))
                norms = np.linalg.norm(exported-physical, axis=1)
                index = int(np.argmax(angles))
                rows.append({"angle_deg": angle, "fd_km": fd, "level": level,
                             "raw_npz_sha256": meta["raw_npz_sha256"], "native_nodes": len(s),
                             "maximum_direction_difference_norm": float(np.max(norms)),
                             "maximum_angle_rad": float(angles[index]), "maximum_at_s_km": float(s[index]),
                             "first_angle_rad": float(angles[0]), "last_angle_rad": float(angles[-1]),
                             "interior_maximum_angle_rad": float(np.max(angles[1:-1])),
                             "first_step_km": float(s[1]-s[0]), "last_step_km": float(s[-1]-s[-2])})
    return {"records": rows,
            "raw_group_tangent_is_normalized_chord": True,
            "raw_group_tangent_is_direct_Hamiltonian_drds": False,
            "physical_B0_isotropic_drds_is_p_over_norm_p": True,
            "physical_Hermite_uses_p_over_norm_p": True,
            "export_field_semantics_need_correction": True,
            "AH_Hamiltonian_error_established_by_this_finding": False,
            "production_patch_implemented": False,
            "proposal": "Export true evaluated drds separately and label the existing normalized chord field explicitly."}


def _refinement_passes(d12: dict, d23: dict) -> bool:
    return all(d23[metric] <= d12[metric] and d23[metric] <= .2*LIMITS[metric] for metric in METRICS)


def _load_reference(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as raw:
        output = {"s": np.asarray(raw["s_km"], dtype=float), "xyz": np.asarray(raw["xyz_km"], dtype=float),
                  "p": np.asarray(raw["p"], dtype=float), "group": np.asarray(raw["group_km"], dtype=float),
                  "tangent": (np.asarray(raw["tangent"], dtype=float) if "tangent" in raw
                              else np.asarray(raw["p"], dtype=float)/np.linalg.norm(raw["p"], axis=1)[:, None])}
    _validate_trace(output, path.name)
    return output


def _reference_at(reference: dict, query: np.ndarray) -> dict:
    indices = np.searchsorted(reference["s"], query)
    if np.any(indices == len(reference["s"])) or np.max(abs(reference["s"][indices]-query)) > 1e-11:
        raise ValueError("native/fine query missing from independently frozen radial-reference arrays")
    return {key: value[indices] for key, value in reference.items()}


def _density_spline():
    path = A2_ROOT / "grids/internal_common.mat"
    raw = loadmat(path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
    height, density = raw["height_km"], raw["ne_profile_m3"]
    return CubicHermiteSpline(height, density, np.gradient(density, height), extrapolate=False)


def _field_extrema(spline, lower: float, upper: float) -> tuple[float, float]:
    density_max, gradient_max = 0., 0.
    for index in range(len(spline.x)-1):
        left, right = float(spline.x[index]), float(spline.x[index+1])
        if right < lower or left > upper:
            continue
        lo, hi = max(lower, left)-left, min(upper, right)-left
        c = spline.c[:, index]
        candidates = [lo, hi]
        for root in np.roots([3*c[0], 2*c[1], c[2]]):
            if abs(float(np.imag(root))) < 1e-12 and lo < float(np.real(root)) < hi:
                candidates.append(float(np.real(root)))
        density_max = max(density_max, *(float(np.polyval(c, point)) for point in candidates))
        candidates = [lo, hi]
        if c[0] != 0.:
            root = -c[1]/(3*c[0])
            if lo < root < hi:
                candidates.append(root)
        gradient_max = max(gradient_max, *(abs(float(np.polyval([3*c[0], 2*c[1], c[2]], point)))
                                           for point in candidates))
    return density_max, gradient_max


def _second_derivative_max(spline) -> float:
    maximum = 0.
    for index in range(len(spline.x)-1):
        if spline.x[index] >= S_KM:
            break
        width = min(S_KM, spline.x[index+1])-spline.x[index]
        for coordinate in (0., width):
            value = 6*spline.c[0, index]*coordinate+2*spline.c[1, index]
            maximum = max(maximum, float(np.linalg.norm(np.atleast_1d(value))))
    return maximum


def _verify_reference(run: Path, original: dict, errors: list[str], root: Path) -> dict:
    report = _json(root / "radial_ODE_controls.json")
    statuses = {float(r["angle_deg"]): r for r in report["reference_status"]}
    uncertainty_report = _json(root / "reference_numerical_uncertainty.json") if (root / "reference_numerical_uncertainty.json").is_file() else None
    uncertainty_cases = {float(c["angle_deg"]): c for c in uncertainty_report["cases"]} if uncertainty_report else {}
    if uncertainty_report:
        _check(uncertainty_report["numerical_estimates_not_certified_intervals"] is True
               and uncertainty_report["formal_multiplier"] == 1.
               and uncertainty_report["historical_v1_multiplier"] == 8.
               and uncertainty_report["arbitrary_allowance_in_v2_budget"] is False,
               "P2 v2 explicit numerical-estimate provenance", errors)
    spline = _density_spline()
    alpha = 80.61638494312905/(10e6)**2
    output = {}
    for angle in (22., 39.5):
        status = statuses[angle]
        active = status["quadrature_control"].get("active_window_levels", [1, 2, 3])
        _check(active in ([1, 2, 3], [3, 4, 5]), "prespecified three-level quadrature window", errors)
        reference = [_load_reference(root / f"reference/e{angle:g}_level{level}.npz") for level in active]
        q = np.linspace(0., 100., 401)
        for trace in reference:
            _check(np.array_equal(trace["s"], q), "radial reference fixed query100 grid", errors)
            _check(np.max(abs(trace["tangent"]-trace["p"]/np.linalg.norm(trace["p"], axis=1)[:, None])) < 1e-12,
                   "radial reference normalized dr/ds", errors)
        d12 = _pair(reference[0], reference[1], "derivative_hermite", .25)
        d23 = _pair(reference[1], reference[2], "derivative_hermite", .25)
        _near_metrics(status["quadrature_control"]["delta12"], d12, f"reference e{angle:g} d12", errors)
        _near_metrics(status["quadrature_control"]["delta23"], d23, f"reference e{angle:g} d23", errors)
        quad_ok = _refinement_passes(d12, d23)
        _check(status["quadrature_control"]["passed"] == quad_ok, "strict quadrature convergence", errors)
        if active == [3, 4, 5]:
            control = status["quadrature_control"]
            _check(control["original_three_profiles_and_outputs_preserved"] is True
                   and control["extra_profile_count"] == 2 and control["no_additional_scan_permitted"] is True,
                   "two fixed additional quadrature controls preserve the original failure", errors)
            old_refs = [_load_reference(root / f"reference/e{angle:g}_level{i}.npz") for i in (1, 2, 3)]
            for pair, key in (((0, 1), "delta12"), ((1, 2), "delta23")):
                old_delta = _pair(old_refs[pair[0]], old_refs[pair[1]], "derivative_hermite", .25)
                _near_metrics(control["original_window_control"][key], old_delta, "preserved original quad "+key, errors)
            all_refs = [_load_reference(root / f"reference/e{angle:g}_level{i}.npz") for i in (1, 2, 3, 4, 5)]
            for increment in control["all_increment_controls"]:
                computed = _pair(all_refs[increment["from_level"]-1], all_refs[increment["to_level"]-1],
                                 "derivative_hermite", .25)
                _near_metrics(increment["metrics"], computed, "all fixed quad increments", errors)
        ph = original[("identified", angle, 3)]["ph"]
        initial = np.concatenate((ph["xyz"][0], ph["p"][0], [0.]))

        def derivative(_s, state):
            xyz, p = state[:3], state[3:6]
            radius, norm = np.linalg.norm(xyz), np.linalg.norm(p)
            gradient = float(spline.derivative()(radius-R_KM))*xyz/radius
            return np.concatenate((p/norm, -alpha*gradient/(2*norm), [1./norm]))

        # The checker uses Cartesian Hamiltonian DOP853. The analysis reference
        # uses radial invariants and quadrature; no analysis code import.
        solved = solve_ivp(derivative, (0., 100.), initial, method="DOP853", rtol=2e-12,
                           atol=2e-13, max_step=.25, t_eval=q)
        if not solved.success or solved.y.shape != (7, 401):
            raise ValueError("independent Cartesian physical-reference control failed")
        cross_max_m = float(np.max(np.linalg.norm(reference[2]["xyz"]-solved.y[:3].T, axis=1)))*1000.
        cross_group_rel = abs(reference[2]["group"][-1]-solved.y[6, -1])/abs(reference[2]["group"][-1])
        _check(cross_max_m < .01 and cross_group_rel < 1e-8,
               "independent Cartesian-vs-radial reference control", errors)
        fine = _load_reference(root / f"reference/native_queries_e{angle:g}.npz")
        n2 = 1.-alpha*spline(np.linalg.norm(fine["xyz"], axis=1)-R_KM)
        energy_drift = float(np.max(abs(np.sum(fine["p"]**2, axis=1)-n2)))
        angular = np.cross(fine["xyz"], fine["p"])
        angular_drift = float(np.max(np.linalg.norm(angular-angular[0], axis=1))/np.linalg.norm(angular[0]))
        invariants_ok = energy_drift <= 1e-10 and angular_drift <= 1e-10
        reported_cross = status["cross_check"]
        dop_paths = [root / f"reference/dop_e{angle:g}_level{i}.npz" for i in (1, 2, 3)]
        if all(path.is_file() for path in dop_paths):
            dop = [_load_reference(path) for path in dop_paths]
            for i, key in ((0, "delta12"), (1, "delta23")):
                _near_metrics(status["DOP_control"][key], _pair(dop[i], dop[i+1], "derivative_hermite", .25),
                              "independently archived DOP control "+key, errors)
            _near_metrics(reported_cross, _pair(reference[2], dop[2], "derivative_hermite", .25),
                          "quad/DOP archived cross check", errors)
        cross_ok = reported_cross["max_m"] <= status["cross_check_limits"]["max_m"] and reported_cross["group_rel"] <= status["cross_check_limits"]["group_rel"]
        established = bool(quad_ok and cross_ok and invariants_ok and not status["quad_warnings"])
        _check(status["established"] == established, "independent reference establishment", errors)
        output[angle] = {"established": established, "fine": fine, "delta12": d12, "delta23": d23,
                         "active_quad_window": active,
                         "quad_uncertainty_support": uncertainty_cases.get(angle),
                         "Cartesian_cross_max_m": cross_max_m, "Cartesian_cross_group_rel": cross_group_rel,
                         "energy_drift": energy_drift, "angular_momentum_vector_relative_drift": angular_drift}
    return output


_QUAD_SUPPORT_CACHE = {}


def _verify_quad_estimate(source: dict, model: dict, alpha: float, errors: list[str], label: str) -> dict:
    """Re-evaluate backend estimates from the recorded physical query and tolerances."""
    if source["s_km"] == 0.:
        return {"position_uncertainty_m": 0., "group_uncertainty_km": 0.}
    radius = source["radius_km"]
    start = float(np.linalg.norm(model["xyz"][0]))
    density = _density_spline()
    energy = float(model["p"][0]@model["p"][0])+alpha*float(density(start-R_KM))
    angular = float(np.linalg.norm(np.cross(model["xyz"][0], model["p"][0])))
    q0 = float(np.sqrt(energy-alpha*float(density(radius-R_KM))-(angular/radius)**2))
    _check(np.isclose(source["L_km"], angular, rtol=1e-12, atol=1e-10)
           and np.isclose(source["q_sqrt"], q0, rtol=1e-11, atol=1e-12),
           label+": quad estimate physical invariant", errors)
    _check(source["root_variable"] == "ECEF_radius_km" and source["root_rtol"] == 4*np.finfo(float).eps,
           label+": recorded root tolerance semantics", errors)
    cache_key = (start, radius, energy, angular, source["epsabs"], source["epsrel"])
    if cache_key not in _QUAD_SUPPORT_CACHE:
        knots = density.x+R_KM
        points = knots[(knots > start) & (knots < radius)]

        def integrand(r, kind):
            n2 = energy-alpha*float(density(r-R_KM))
            radial = np.sqrt(n2-(angular/r)**2)
            return (np.sqrt(n2)/radial if kind == 0 else angular/(r*r*radial) if kind == 1 else 1./radial)

        _QUAD_SUPPORT_CACHE[cache_key] = [quad(integrand, start, radius, args=(i,), epsabs=source["epsabs"],
                                             epsrel=source["epsrel"], points=points, limit=300) for i in (0, 1, 2)]
    values = _QUAD_SUPPORT_CACHE[cache_key]
    _check(abs(values[0][0]-source["s_km"]) <= 4*(source["root_uncertainty_km"]+values[0][1]+1e-12),
           label+": recorded arc-length root query", errors)
    for field, i in (("arc_quad_estimate_km", 0), ("phi_quad_estimate_rad", 1), ("group_quad_estimate_km", 2)):
        _check(np.isclose(source[field], values[i][1], rtol=1e-5,
                          atol=1e-15 if i == 1 else 1e-12), label+": independently repeated backend "+field, errors)
    root_error = source["root_xtol_km"]+source["root_rtol"]*abs(radius)
    dr_error = root_error+source["arc_quad_estimate_km"]
    position_error = 1000.*(dr_error*(1.+angular/(radius*q0))+radius*source["phi_quad_estimate_rad"])
    group_error = source["group_quad_estimate_km"]+dr_error/q0
    _check(np.isclose(source["root_uncertainty_km"], root_error, rtol=1e-12, atol=1e-15)
           and np.isclose(source["radial_uncertainty_estimate_km"], dr_error, rtol=1e-12, atol=1e-15)
           and np.isclose(source["position_uncertainty_m"], position_error, rtol=1e-10, atol=1e-10)
           and np.isclose(source["group_uncertainty_km"], group_error, rtol=1e-10, atol=1e-13),
           label+": backend/root error estimate composition", errors)
    return {"position_uncertainty_m": position_error, "group_uncertainty_km": group_error}


def _verify_reconstruction_proof(run: Path, angle: float, raw: dict, method: str, side: str,
                                  result: dict, label: str, errors: list[str], reference: dict) -> dict:
    proof = result["proof"]
    trace = _sample_native(raw, 1)
    exact_nodes = _reference_at(reference["fine"], trace["s"])
    model = {"s": trace["s"], "xyz": exact_nodes["xyz"], "p": exact_nodes["p"], "group": exact_nodes["group"]}
    q = np.linspace(0., 100., 1601)
    truth = _reference_at(reference["fine"], q)
    if method == "pchip":
        pos = PchipInterpolator(model["s"], model["xyz"], axis=0, extrapolate=False)
        group = PchipInterpolator(model["s"], model["group"], extrapolate=False)
    else:
        norm = np.linalg.norm(model["p"], axis=1)
        pos = CubicHermiteSpline(model["s"], model["xyz"], model["p"]/norm[:, None], axis=0, extrapolate=False)
        group = CubicHermiteSpline(model["s"], model["group"], 1./norm, extrapolate=False)
    actual = np.linalg.norm(pos(q)-truth["xyz"], axis=1)*1000.
    actual_G = np.abs(group(q)-truth["group"])
    position_derivative = np.linalg.norm(pos.derivative()(q)-truth["tangent"], axis=1)
    group_derivative = np.abs(group.derivative()(q)-1./np.linalg.norm(truth["p"], axis=1))
    alpha = 80.61638494312905/(10e6)**2
    _check(np.isclose(proof["alpha"], alpha, rtol=1e-14, atol=0.)
           and np.isclose(proof["C"], 80.61638494312905, rtol=1e-14, atol=0.)
           and proof["R_km"] == R_KM, label+": reconstruction physical constants", errors)
    lower, upper = proof["height_min_km"], proof["height_max_km"]
    _check(lower <= 100. and upper >= 100.+float(trace["s"][-1]), label+": field bound radial reach", errors)
    density_max, gradient_max = _field_extrema(_density_spline(), lower, upper)
    energy = float(model["p"][0]@model["p"][0])+alpha*float(_density_spline()(100.))
    n2_min = energy-alpha*density_max
    _check(n2_min > 0. and np.isclose(proof["n2_min"], n2_min, rtol=1e-11, atol=1e-13)
           and np.isclose(proof["Neprime_max_m3_per_km"], gradient_max, rtol=1e-11, atol=1e-3),
           label+": independent Hermite density/gradient extrema", errors)
    M2, GM2 = alpha*gradient_max/(2*n2_min), alpha*gradient_max/(2*n2_min**1.5)
    _check(np.isclose(proof["position_M2_per_km"], M2, rtol=1e-11, atol=1e-13)
           and np.isclose(proof["group_M2_per_km"], GM2, rtol=1e-11, atol=1e-13)
           and proof["h4_remainder_used"] is False, label+": C1-media second-derivative bounds", errors)
    second, gsecond = _second_derivative_max(pos), _second_derivative_max(group)
    fill = .0625/2*float(np.max(position_derivative))+.0625**2/8*(second+M2)
    group_fill = .0625/2*float(np.max(group_derivative))+.0625**2/8*(gsecond+GM2)
    assertions = {"node_maxgap_km": float(np.max(np.diff(trace["s"]))), "validation_step_km": .0625,
                  "poly_second_derivative_max": second, "group_poly_second_derivative_max": gsecond,
                  "fill_bound_m": fill*1000., "group_fill_bound_km": group_fill,
                  "actual_reference_interpolation_rms_m": float(np.sqrt(np.mean(actual**2))),
                  "actual_reference_interpolation_max_m": float(np.max(actual)),
                  "actual_reference_interpolation_endpoint_m": float(actual[-1]),
                  "actual_reference_group_interpolation_max_km": float(np.max(actual_G)),
                  "actual_reference_group_interpolation_endpoint_km": float(actual_G[-1])}
    for key, value in assertions.items():
        _check(np.isclose(proof[key], value, rtol=1e-8, atol=1e-10 if key.endswith("km") else 1e-7),
               label+": reference-only reconstruction proof "+key, errors)
    _check(proof["node_error_in_reconstruction_budget"] is False
           and proof["reference_samples_are_exact_model_queries"] is True
           and proof["quadrature_errors_are_numerical_estimates"] is True,
           label+": node propagation/reference reconstruction separation", errors)
    _check(result["bound_established"] == reference["established"], label+": reference-control support", errors)
    _, denominator = _reconstruct(trace, np.array([100.]), method)
    uncertainty, Guncertainty = proof["quad_uncertainty_m"], proof["group_quad_uncertainty_km"]
    _check(np.isfinite(uncertainty) and uncertainty >= 0. and np.isfinite(Guncertainty) and Guncertainty >= 0.,
           label+": explicit numerical quadrature estimates", errors)
    if "uncertainty_multiplier" in proof:
        support = reference["quad_uncertainty_support"]
        _check(support is not None and support["profile_level"] == reference["active_quad_window"][-1]
               and support["numerical_estimate_is_interval_certified"] is False
               and support["uncertainty_multiplier"] == 1., label+": quad support active profile", errors)
        if support is None:
            raise ValueError("quadrature support log missing")
        position_witness, group_witness = proof["position_uncertainty_source"], proof["group_uncertainty_source"]
        _check(position_witness in support["records"] and group_witness in support["records"],
               label+": uncertainty maxima link to immutable per-query logs", errors)
        _check(np.isclose(support["max_pos_uncertainty_m"],
                          max(row["position_uncertainty_m"] for row in support["records"]), rtol=1e-12, atol=1e-12)
               and np.isclose(support["max_group_uncertainty_km"],
                              max(row["group_uncertainty_km"] for row in support["records"]), rtol=1e-12, atol=1e-15),
               label+": global backend estimate maxima", errors)
        position_support = _verify_quad_estimate(proof["position_uncertainty_source"], model, alpha, errors, label)
        group_support = _verify_quad_estimate(proof["group_uncertainty_source"], model, alpha, errors, label)
        _check(proof["uncertainty_multiplier"] == 1.
               and np.isclose(uncertainty, position_support["position_uncertainty_m"], rtol=1e-9, atol=1e-10)
               and np.isclose(Guncertainty, group_support["group_uncertainty_km"], rtol=1e-9, atol=1e-13)
               and uncertainty == proof["max_pos_uncertainty_m_unmultiplied"]
               and Guncertainty == proof["max_group_uncertainty_km_unmultiplied"],
               label+": single supported quadrature estimates, no arbitrary multiplier", errors)
    else:
        errors.append(label+": raw quadrature estimate support missing")
    computed = {"rms_m": float(np.max(actual))+fill*1000.+uncertainty,
                "max_m": float(np.max(actual))+fill*1000.+uncertainty,
                "landing_m": float(actual[-1])+uncertainty,
                "group_rel": (float(actual_G[-1])+Guncertainty)/abs(denominator) if side == "ph" else 0.}
    _near_metrics(result["bounds"], computed, label+": independently reconstructed error bound", errors)
    raw_positions, _ = _reconstruct(trace, q, method)
    node_propagation = np.linalg.norm(raw_positions-pos(q), axis=1)*1000.
    _check(np.isclose(proof["propagated_node_error_query_rms_m"], np.sqrt(np.mean(node_propagation**2)), rtol=1e-8, atol=1e-7)
           and np.isclose(proof["propagated_node_error_query_max_m"], np.max(node_propagation), rtol=1e-8, atol=1e-7),
           label+": propagated-node differences reported separately", errors)
    return computed


def _verify_numerical(run: Path, original: dict, errors: list[str], root: Path) -> dict:
    gates = _json(root / "strict_gates.json")
    _check(gates["limits"] == LIMITS, "P2 strict thresholds unchanged", errors)
    _check(gates["sample_steps_km"] == [1., .5, .25]
           and gates["node_sampling_strides"] == [4, 2, 1], "P2 prespecified sampling controls", errors)
    groups = gates["groups"]
    references = _verify_reference(run, original, errors, root)
    expected_groups = {(angle, fd, method) for angle in (22., 39.5) for fd in (.5, .125)
                       for method in ("pchip", "derivative_hermite")}
    _check(len(groups) == 8 and {(g["angle_deg"], g["fd_km"], g["method"]) for g in groups} == expected_groups,
           "P2 all original/candidate angle/method groups", errors)
    rows = []
    for group in groups:
        label, angle, fd, method = group["group_id"], float(group["angle_deg"]), group["fd_km"], group["method"]
        _check(group["arc_km"] == S_KM and group["role"] == ("original_fd0.5" if fd == .5 else "candidate_fd0.125"),
               label+": original angle/FD/full100km contract", errors)
        levels = sorted(group["levels"], key=lambda row: row["level"])
        _check([row["level"] for row in levels] == [1, 2, 3], label+": three integration levels", errors)
        traces = {}
        for row in levels:
            level = int(row["level"])
            source = original[("identified", angle, level)]
            expected_job = f"internal_identified_f10_e{'22' if angle == 22. else '39p5'}_t{level}"
            _check(Path(row["PH_raw"]) == A2_ROOT / "identified/pharlap_raw" / (expected_job+".mat"),
                   label+": original PH native path", errors)
            expected_ah = (A2_ROOT / "identified/ah/formal/identified/hermite" / (expected_job+".json") if fd == .5
                           else run / "ah/fd0.125" / (expected_job+".json"))
            _check(Path(row["AH_meta"]) == expected_ah and Path(row["AH_raw"]) == expected_ah.with_suffix(".npz"),
                   label+": AH fixed-FD native path", errors)
            if fd == .5:
                ah, meta = _ah(expected_ah, 10.)
            else:
                ah, meta = _ah(expected_ah, 10.)
                _check(meta["coefficient_K2"] == 80.61638494312905
                       and meta["grid_sha256"] == INTERNAL_GRID_SHA and meta["interpolation"] == "hermite",
                       label+": candidate medium identity", errors)
                _check(_initial_matches(source["ph"]["xyz"][0], source["ph"]["p"][0], ah),
                       label+": candidate matched PH initial state", errors)
            _check(row["status"] == "completed", label+": recorded native completion", errors)
            traces[level] = {"ph": source["ph"], "ah": ah}
        actual = _pair(traces[3]["ph"], traces[3]["ah"], method, .25)
        _near_metrics(group["metrics"], actual, label+": cross residual", errors)
        d12, d23, sampling12, sampling23 = {}, {}, {}, {}
        for side in ("ph", "ah"):
            d12[side] = _pair(traces[1][side], traces[2][side], method, .25)
            d23[side] = _pair(traces[2][side], traces[3][side], method, .25)
            _near_metrics(group["integration"][side]["delta12"], d12[side], label+": integration "+side+" 12", errors)
            _near_metrics(group["integration"][side]["delta23"], d23[side], label+": integration "+side+" 23", errors)
            sampled = [_sample_native(traces[3][side], stride) for stride in (4, 2, 1)]
            sampling12[side] = _pair(sampled[0], sampled[1], method, .25)
            sampling23[side] = _pair(sampled[1], sampled[2], method, .25)
            _near_metrics(group["sampling"][side]["delta12"], sampling12[side], label+": native sampling "+side+" 12", errors)
            _near_metrics(group["sampling"][side]["delta23"], sampling23[side], label+": native sampling "+side+" 23", errors)
        _check([row["step_km"] for row in group["query_sampling"]] == [1., .5, .25],
               label+": three query-grid controls", errors)
        query_controls = []
        for row in group["query_sampling"]:
            computed = _pair(traces[3]["ph"], traces[3]["ah"], method, row["step_km"])
            _near_metrics(row["metrics"], computed, label+": query grid "+str(row["step_km"]), errors)
            query_controls.append(computed)
        reconstructed = {}
        for side in ("ph", "ah"):
            bound = group["reconstruction"][side]
            _check(type(bound["bound_established"]) is bool, label+": reconstruction establishment status", errors)
            reconstructed[side] = _verify_reconstruction_proof(run, angle, traces[3][side], method,
                                                               side, bound, label, errors, references[angle])
        resampling = {metric: reconstructed["ph"][metric]+reconstructed["ah"][metric] for metric in METRICS}
        budget = {metric: d23["ph"][metric]+d23["ah"][metric]+resampling[metric] for metric in METRICS}
        _near_metrics(group["budget"], budget, label+": supported single numerical budget", errors)
        integration_ok = all(_refinement_passes(d12[side], d23[side]) for side in ("ph", "ah"))
        sampling_ok = all(_refinement_passes(sampling12[side], sampling23[side]) for side in ("ph", "ah"))
        bounds_ok = all(group["reconstruction"][side]["bound_established"] for side in ("ph", "ah"))
        reference_ok = references[angle]["established"]
        eligible = bool(integration_ok and sampling_ok and bounds_ok and reference_ok)
        _check(group["method_difference_in_budget"] is False and group["FD_difference_in_budget"] is False
               and group["node_error_in_reconstruction_budget"] is False
               and group["AH_group_kind"] == "native_total_at_exact_100km"
               and group["PH_group_kind"] == "native_G_nodes_reconstructed_at_100_per_method",
               label+": native group and supported-budget scope", errors)
        _check(group["eligible"] == eligible, label+": independent eligibility", errors)
        decision = _decision(actual, d12, d23, resampling, eligible)
        _check(group["strict_passed"] == decision["passed"], label+": strict decision", errors)
        _check(group["decision"]["passed"] == decision["passed"] and group["decision"]["eligible"] == eligible,
               label+": exact frozen P0 gate decision", errors)
        for metric in METRICS:
            recorded = group["decision"]["checks"][metric]
            computed = decision["checks"][metric]
            _check(all(recorded[field] == computed[field] for field in ("physical_pass", "explained_by_numeric", "sides")),
                   label+": unrelaxed frozen P0 check "+metric, errors)
            _check(np.isclose(recorded["three_times_budget"], computed["three_times_budget"], rtol=2e-9,
                              atol=1e-14 if metric == "group_rel" else 1e-7),
                   label+": 3x supported budget "+metric, errors)
        rows.append({"group_id": label, "angle_deg": angle, "fd_km": fd, "method": method,
                     "computed_actual": actual, "integration_d12": d12, "integration_d23": d23,
                     "sampling_d12": sampling12, "sampling_d23": sampling23,
                     "query_controls": query_controls, "reconstruction": reconstructed,
                     "integration_convergent": integration_ok, "native_sampling_convergent": sampling_ok,
                     "reference_controls_passed": reference_ok, "supported_bounds_established": bounds_ok,
                     "single_budget": budget, "independent_P0_decision": decision,
                     "failure_checks": [metric for metric, check in decision["checks"].items()
                                        if not check["physical_pass"] or not check["explained_by_numeric"]
                                        or not all(all(v.values()) for v in check["sides"].values())]})
    return {"groups": rows, "strict_passed_count": sum(row["independent_P0_decision"]["passed"] for row in rows),
            "strict_passed": bool(rows) and all(row["independent_P0_decision"]["passed"] for row in rows),
            "physical_thresholds": LIMITS, "roundoff_waiver_applied": False,
            "interpolation_method_difference_in_budget": False,
            "node_propagation_error_added_to_resampling_budget": False,
            "reference_controls": {str(angle): {key: value for key, value in ref.items()
                                                if key not in ("fine", "quad_uncertainty_support")}
                                   for angle, ref in references.items()}}


def self_check() -> list[str]:
    errors = []
    zero = dict.fromkeys(METRICS, 0.)
    good = {key: value*1e-3 for key, value in LIMITS.items()}
    first = {side: {key: value*2 for key, value in good.items()} for side in ("ph", "ah")}
    last = {side: dict(good) for side in ("ph", "ah")}
    _check(_decision(good, first, last, zero, True)["passed"], "selfcheck supported convergent case", errors)
    bad = {side: dict(good) for side in ("ph", "ah")}
    bad["ph"]["group_rel"] *= 3
    _check(not _decision(good, first, bad, zero, True)["passed"], "selfcheck PH nonconvergence", errors)
    _check(not _decision(good, first, last, zero, False)["passed"], "selfcheck unsupported bound", errors)
    s = np.array([0., 50., 100.])
    unit = np.array([.6, .8, 0.])
    trace = {"s": s, "xyz": s[:, None]*unit, "p": np.tile(.9*unit, (3, 1)), "group": s/.9}
    positions, group = _reconstruct(trace, np.linspace(0., 100., 201), "direction_hermite")
    _check(np.max(abs(positions-np.linspace(0., 100., 201)[:, None]*unit)) < 3e-14,
           "selfcheck derivative is unit tangent", errors)
    _check(abs(group-100./.9) < 1e-13, "selfcheck native group reconstruction", errors)
    bound = _reconstruction_upper(trace, "direction_hermite", 0., 0.)
    _check(bound["position_upper_m"] < 1e-10 and bound["group_upper_km"] < 1e-12,
           "selfcheck linear deterministic reconstruction bound", errors)
    radius = 250.
    angle = s/radius
    circle = {"s": s, "xyz": np.column_stack((radius*np.sin(angle), radius*(1.-np.cos(angle)), s*0.)),
              "p": .9*np.column_stack((np.cos(angle), np.sin(angle), s*0.)), "group": s/.9}
    query = np.linspace(0., S_KM, 10001)
    truth = np.column_stack((radius*np.sin(query/radius), radius*(1.-np.cos(query/radius)), query*0.))
    for method in ("pchip", "direction_hermite"):
        pred, _ = _reconstruct(circle, query, method)
        observed = float(np.max(np.linalg.norm(pred-truth, axis=1)))*1000.
        upper = _reconstruction_upper(circle, method, 1./radius, 0.)["position_upper_m"]
        _check(observed <= upper, "selfcheck analytic circle reconstruction bound "+method, errors)
    try:
        _reconstruct(trace, np.linspace(0., 99., 201), "direction_hermite")
    except ValueError:
        pass
    else:
        errors.append("selfcheck shortened comparison allowed")
    short = {**trace, "s": np.array([0., 50., 99.])}
    try:
        _validate_trace(short, "selfcheck shortened native path")
    except ValueError:
        pass
    else:
        errors.append("selfcheck shortened native path allowed")
    _check(not _initial_matches(trace["xyz"][0], trace["p"][0]+np.array([0., 0., 1e-5]),
                                trace), "selfcheck mismatched initial state", errors)
    scratch = CODE_ROOT / "output"
    with tempfile.TemporaryDirectory(prefix="c1m_p2_verifier_", dir=scratch) as folder:
        root = Path(folder)
        root.resolve().relative_to(scratch.resolve())
        raw = root / "raw.bin"
        raw.write_bytes(b"native-state-fixture")
        manifest = root / "delivery_sha256.json"
        manifest.write_text(json.dumps({"files": [{"path": raw.name, "sha256": _sha(raw),
                                                  "bytes": raw.stat().st_size}]}), encoding="utf-8")
        frozen = _sha(manifest)
        clean = []
        _verify_delivery(root, manifest.name, frozen, 1, clean)
        _check(not clean, "selfcheck untampered delivery", errors)
        raw.write_bytes(b"tampered-state")
        tampered = []
        _verify_delivery(root, manifest.name, frozen, 1, tampered)
        _check(bool(tampered), "selfcheck archived data tamper detection", errors)
    return errors


def _verify_impl(run: Path, output: Path, analysis_root: Path) -> list[str]:
    report = {"schema": "a4-c1m-p2-independent-verification-v1", "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source_verified": False, "contract_verified": False, "strict_passed": False,
              "global_physics_accepted": False, "old_A2_preparation_provenance_gap_resolved": False,
              "scope_limits": [
                  "Reference exact-node interpolation bounds exclude propagated native-node errors.",
                  "Quadrature backend error estimates are numerical estimates, not interval certificates.",
                  "Source integrity verified now does not establish old A2 pre-execution source freezing.",
                  "The exit core normalization location and first crossing interval remain unresolved."]}
    errors = []
    source_end = 0
    for name, function in (
        ("self_check", lambda: {"errors": self_check()}),
        ("P1_delivery", lambda: _verify_delivery(P1_ROOT, "delivery_sha256.json", P1_DELIVERY_SHA, 1285, errors)),
        ("legacy_A2", lambda: _verify_legacy(errors)),
        ("P2_component_sources", lambda: _verify_p2_sources(run, errors, analysis_root)),
        ("production_and_MEX", lambda: _verify_production(run, errors)),
        ("exit_complete_sequences", lambda: _verify_exit(run, errors)),
    ):
        try:
            result = function()
            if name == "legacy_A2":
                report[name], original = result
            else:
                report[name] = result
            if name == "self_check":
                errors.extend(result["errors"])
        except Exception as exception:
            errors.append(name+": "+type(exception).__name__+": "+str(exception))
            report[name] = {"exception": traceback.format_exc()}
        source_end = len(errors)
    report["source_verified"] = source_end == 0
    try:
        if "original" not in locals():
            raise ValueError("Original A2 native pairs could not be loaded")
        report["AH_derivative_export_semantics"] = _verify_AH_derivative_semantics(run, original, errors)
        report["numerical_verification"] = _verify_numerical(run, original, errors, analysis_root)
        report["strict_passed"] = bool(report["numerical_verification"]["strict_passed"])
    except Exception as exception:
        errors.append("numerical_verification: "+type(exception).__name__+": "+str(exception))
        report["numerical_verification"] = {"exception": traceback.format_exc()}
    report["contract_verified"] = len(errors) == source_end
    report["errors"] = errors
    report["error_count"] = len(errors)
    path = output / "verification.json"
    path.write_text(json.dumps(report, ensure_ascii=True, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return errors


def verify(run: Path, analysis_directory: str = "analysis/p2_independent") -> list[str]:
    """Freeze a complete verifier and execute that snapshot; retain every version."""
    run = Path(run).resolve()
    root = run / analysis_directory
    root.resolve().relative_to(run.resolve())
    version = 1
    while (root / f"verification_v{version}").exists():
        version += 1
    output = root / f"verification_v{version}"
    sources = output / "sources"
    sources.mkdir(parents=True)
    original = Path(__file__).resolve()
    snapshot = sources / original.name
    shutil.copy2(original, snapshot)
    manifest = {"schema": "a4-c1m-p2-verification-source-v1",
                "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "version": version,
                "source_path": str(original), "source_sha256": _sha(original),
                "snapshot": str(snapshot.relative_to(output)), "snapshot_sha256": _sha(snapshot),
                "prepare_manifest_sha256": _sha(run / "p2_manifest.json"),
                "analysis_manifest_sha256": _sha(root / "analysis_manifest.json"),
                "execution_uses_frozen_snapshot": True, "earlier_verification_versions_preserved": True}
    if manifest["source_sha256"] != manifest["snapshot_sha256"]:
        raise ValueError("verifier source changed during freeze")
    manifest_path = output / "verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    spec = importlib.util.spec_from_file_location(f"c1m_p2_frozen_verifier_v{version}", snapshot)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    errors = module._verify_impl(run, output, root)
    index = {"schema": "a4-c1m-p2-verification-index-v1", "version": version,
             "manifest": str(manifest_path), "manifest_sha256": _sha(manifest_path),
             "result": str(output / "verification.json"), "result_sha256": _sha(output / "verification.json"),
             "errors": len(errors), "source_execution_from_snapshot": True}
    (root / f"verification_index_v{version}.json").write_text(json.dumps(index, indent=2)+"\n", encoding="utf-8")
    # The first index is mirrored for the original main-entry contract. Later
    # versions retain their own results and do not overwrite an earlier failure.
    for name in ("verification_manifest.json", "verification.json"):
        if not (root / name).exists():
            shutil.copy2(output / name, root / name)
    return errors
