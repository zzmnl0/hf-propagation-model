"""Independent read-only P3 evidence and numerical acceptance verification."""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.io import loadmat

LIMITS = {"rms_m": 10., "max_m": 50., "landing_m": 50., "group_rel": 1e-5}
PH_PROFILES = ((1e-8, .01, 10.), (1e-10, .001, 5.), (1e-12, .001, 1.))


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def local_import_closure(entries, root):
    """Follow local static imports without importing or executing experiment code."""
    todo, visited = [Path(p).resolve() for p in entries], set()
    while todo:
        path = todo.pop()
        if path in visited:
            continue
        visited.add(path)
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        for node in ast.walk(tree):
            names = []
            bases = [path.parent, root, root/"tests"]
            if isinstance(node, ast.Import):
                names = [x.name for x in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module] if node.module else []
                if node.level:
                    base = path.parent
                    for _ in range(node.level-1):
                        base = base.parent
                    bases = [base]
            for name in names:
                for base in bases:
                    target = base.joinpath(*name.split("."))
                    candidates = (target.with_suffix(".py"), target/"__init__.py")
                    found = next((p for p in candidates if p.is_file() and p.is_relative_to(root)), None)
                    if found:
                        todo.append(found.resolve())
                        break
    return sorted(visited)


def vector(obj, name):
    return np.asarray(getattr(obj, name), dtype=float).reshape(-1)


def ph_trace(path):
    raw = loadmat(path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
    state, output = raw["ray_state_vec"], raw["ray_path_data"]
    trace = {
        "s": vector(state, "geometrical_path") / 1000.,
        "xyz": np.column_stack([vector(state, n) for n in ("pos_x", "pos_y", "pos_z")]) / 1000.,
        "p": np.column_stack([vector(state, n) for n in ("dir_x", "dir_y", "dir_z")]),
        "group": vector(state, "group_path") / 1000.,
        "termination": int(np.asarray(raw["ray_data"].ray_label).reshape(-1)[-1]),
    }
    trace["s"] -= trace["s"][0]
    trace["group"] -= trace["group"][0]
    for name in ("geomag_x", "geomag_y", "geomag_z", "collision_frequency"):
        if np.any(vector(output, name) != 0.):
            raise ValueError("PH violates B0 collisionless scope: " + name)
    return trace, raw


def ah_trace(path):
    with np.load(path, allow_pickle=False) as raw:
        return {k: np.array(raw[k], dtype=float) for k in ("s", "xyz", "p", "group")}


def validate_trace(t, scope):
    s = t["s"]
    if (s.ndim != 1 or len(s) < 2 or np.any(np.diff(s) <= 0.)
            or abs(s[0]) > 1e-10 or t["xyz"].shape != (len(s), 3)
            or t["p"].shape != t["xyz"].shape or t["group"].shape != s.shape
            or any(not np.all(np.isfinite(t[k])) for k in ("s", "xyz", "p", "group"))
            or np.any(np.linalg.norm(t["p"], axis=1) <= 0.)
            or np.any(np.diff(t["group"]) <= 0.)):
        raise ValueError("Invalid native trace")
    if scope == "first_100_km" and s[-1] < 100. - 1e-9:
        raise ValueError("Native trace does not cover 100 km")


def curves(t):
    norm = np.linalg.norm(t["p"], axis=1)
    return (CubicHermiteSpline(t["s"], t["xyz"], t["p"] / norm[:, None], axis=0, extrapolate=False),
            CubicHermiteSpline(t["s"], t["group"], 1. / norm, extrapolate=False))


def metrics(a, b, scope, spacing=.25):
    """Full-path endpoints use each natural end; local endpoints both use 100 km."""
    end = 100. if scope == "first_100_km" else min(a["s"][-1], b["s"][-1])
    q = np.linspace(0., end, int(np.ceil(end / spacing)) + 1)
    ap, ag = curves(a)
    bp, bg = curves(b)
    distance = np.linalg.norm(ap(q) - bp(q), axis=1) * 1000.
    if scope == "first_100_km":
        endpoint = np.linalg.norm(ap(100.) - bp(100.)) * 1000.
        ga, gb = float(ag(100.)), float(bg(100.))
    else:
        endpoint = np.linalg.norm(a["xyz"][-1] - b["xyz"][-1]) * 1000.
        ga, gb = float(a["group"][-1]), float(b["group"][-1])
    return {"rms_m": float(np.sqrt(np.mean(distance**2))), "max_m": float(np.max(distance)),
            "landing_m": float(endpoint), "group_rel": abs(ga-gb)/max(abs(ga), 1e-15)}


def decision(actual, controls, resampling, eligible):
    """No propagation bias or method difference enters the numerical budget."""
    by_metric = {}
    for key, limit in LIMITS.items():
        finite = all(np.isfinite(v) and v >= 0. for v in (
            actual[key], resampling[key], *[c[d][key] for c in controls for d in ("delta12", "delta23")]))
        convergence = all(c["delta23"][key] <= c["delta12"][key]
                          and c["delta23"][key] <= .2*limit for c in controls)
        budget = 3.*(sum(c["delta23"][key] for c in controls) + resampling[key])
        by_metric[key] = {"budget": budget, "passed": bool(finite and convergence
            and resampling[key] <= .2*limit and actual[key] <= limit and actual[key] <= budget)}
    return {"metrics": by_metric, "passed": bool(eligible and all(x["passed"] for x in by_metric.values()))}


def self_check():
    s = np.arange(0., 102., 1.)
    straight = {"s": s, "xyz": np.column_stack((s, s*0, s*0)),
                "p": np.tile([.8, 0., 0.], (len(s), 1)), "group": s/.8}
    validate_trace(straight, "first_100_km")
    assert all(v == 0. for v in metrics(straight, straight, "first_100_km").values())
    short = {k: v[:-2] for k, v in straight.items()}
    try:
        validate_trace(short, "first_100_km")
    except ValueError:
        pass
    else:
        raise AssertionError("premature termination accepted")
    full = metrics(straight, short, "natural_full_path")
    assert full["rms_m"] == 0. and full["landing_m"] == 2000.
    assert full["group_rel"] > .01
    theta = s/500.
    circle = {"s": s, "xyz": np.column_stack((500.*np.sin(theta), 500.*(1-np.cos(theta)), s*0)),
              "p": .8*np.column_stack((np.cos(theta), np.sin(theta), s*0)), "group": s/.8}
    q = np.arange(.5, 100., 1.)
    exact = np.column_stack((500.*np.sin(q/500.), 500.*(1-np.cos(q/500.)), q*0))
    assert np.max(np.linalg.norm(curves(circle)[0](q)-exact, axis=1))*1000. < 1e-6
    small = {k: v*.01 for k, v in LIMITS.items()}
    zero = {k: 0. for k in LIMITS}
    controls = [{"delta12": small, "delta23": zero}]*2
    assert decision(zero, controls, zero, True)["passed"] is True
    assert decision(zero, controls, zero, False)["passed"] is False
    assert decision(small, controls, zero, True)["passed"] is False
    assert decision(zero, [{"delta12": zero, "delta23": small}]*2, zero, True)["passed"] is False
    return {"straight": True, "curved": True, "premature_termination": True,
            "natural_endpoints": True, "budget_failure": True, "nonmonotonic_failure": True}


def verify_r5(run):
    directory = run/"experiments"
    manifest = read(directory/"manifest.json")
    jobs = manifest["jobs"]
    rows = []

    def sequence(items, side, label, kind, settings):
        result = {"case_id": label, "kind": kind, "side": side, "settings": settings,
                  "jobs": [j["job_id"] for j in items], "passed": False}
        try:
            traces = []
            for job in items:
                if side == "ah":
                    meta = read(directory/"ah"/(job["job_id"]+".json"))
                    if meta["status"] != "completed":
                        raise ValueError(meta.get("error", "AH execution failed"))
                    path = directory/"ah"/(job["job_id"]+".npz")
                    if sha(path) != meta["raw_sha256"]:
                        raise ValueError("AH hash mismatch")
                    trace = ah_trace(path)
                else:
                    trace, _ = ph_trace(directory/"pharlap_raw"/(job["job_id"]+".mat"))
                validate_trace(trace, job["scope"])
                traces.append(trace)
            changes = [metrics(a, b, items[0]["scope"]) for a, b in zip(traces[:-1], traces[1:])]
            monotonic = {k: all(b[k] <= a[k] for a, b in zip(changes[:-1], changes[1:])) for k in LIMITS}
            small = {k: changes[-1][k] <= .2*LIMITS[k] for k in LIMITS}
            result.update(changes=changes, monotonic=monotonic, last_change_small=small,
                          passed=bool(len(items) >= 3 and all(monotonic.values()) and all(small.values())))
        except (ValueError, OSError, KeyError) as exc:
            result["error"] = str(exc)
        rows.append(result)

    for base in (j for j in jobs if j["role"] in ("R1", "R2") and j["level"] == 3):
        others = [j for j in jobs if j.get("ph_job_id") == base["job_id"]]
        if not others:
            continue
        items = sorted(others+[base], key=lambda j: -j["spatial_fd_base_km"])
        sequence(items, "ah", base["job_id"], "spatial_fd", [j["spatial_fd_base_km"] for j in items])
    groups = {}
    for job in jobs:
        if job["role"] == "R5_grid" and job["level"] == 3:
            groups.setdefault((job["case_id"].split("_d")[0], job["elevation_deg"]), []).append(job)
    for key, group in groups.items():
        step = {j["job_id"]: float(np.diff(loadmat(j["grid_mat"], squeeze_me=True)["height_km"])[0]) for j in group}
        items = sorted(group, key=lambda j: -step[j["job_id"]])
        for side in ("ah", "ph"):
            sequence(items, side, str(key), "height_grid", [step[j["job_id"]] for j in items])
    return {"schema": "c1m-p3-r5-independent-gates-v1", "cases": rows,
            "passed_count": sum(r["passed"] for r in rows), "case_count": len(rows),
            "r5_passed": bool(rows and all(r["passed"] for r in rows)),
            "used_in_numerical_uncertainty_budget": False, "roundoff_waiver_applied": False}


def verify(run, analysis):
    errors, identities, rows = [], [], []

    def check(ok, message):
        if not bool(ok):
            errors.append(message)

    def identity(path, expected, label):
        path = Path(path)
        actual = sha(path) if path.is_file() else None
        identities.append({"path": str(path.resolve()), "expected_sha256": expected,
                           "actual_sha256": actual, "passed": actual == expected})
        check(actual == expected, label + ": hash identity")

    def close_map(actual, reported, label):
        for key in LIMITS:
            check(key in reported and np.isclose(actual[key], reported[key], rtol=2e-8, atol=1e-13),
                  label + ": " + key)

    d = run / "experiments"
    manifest = read(d / "manifest.json")
    summary = read(analysis / "summary.json")
    parent_analysis = Path(summary["parent_summary"]).parent if "parent_summary" in summary else analysis
    prior = read(run / "p0_prior_evidence_verification.json")
    check(prior["passed"] is True and len(prior["records"]) == 219, "P0 prior evidence identities")
    for item in prior["records"]:
        identity(item["path"], item["expected_sha256"], "preserved historical evidence")
    diagnostic_path = run / "production_diagnostics/ph_native_smooth_entry/independent_check.json"
    diagnostic = read(diagnostic_path)
    check(diagnostic["passed"] is True and diagnostic["scientific_scope_passed"] is False,
          "smooth output-density diagnostic does not claim propagation equivalence")
    identity(diagnostic_path.parent / "summary.json", diagnostic["summary_sha256"], "smooth diagnostic report")
    r0 = read(run/"r0_validation.json")
    identity(r0["suite"], r0["suite_sha256"], "R0 production suite")
    suite = read(r0["suite"])
    check(suite["full_suite_passed"] is True and r0["suite_passed"] is True, "complete R0 suite")
    with np.load(run/"before/physics.npz") as previous, np.load(run/"after/physics.npz") as current:
        check(set(previous.files) == set(current.files), "default baseline array coverage")
        for key in previous.files:
            check(np.array_equal(previous[key], current[key]), "default baseline exact equality: "+key)
    check(read(run/"before/physics.json") == read(run/"after/physics.json"), "default baseline metadata")
    check(r0["default_regression"]["passed"] is True
          and r0["power_helper_regression"]["passed"] is True, "default regression boolean gates")
    jobs = {j["job_id"]: j for j in manifest["jobs"]}
    extra_path = d/"refinements_v1/manifest.json"
    extra = read(extra_path) if extra_path.is_file() else {"jobs": []}
    all_jobs = {**jobs, **{j["job_id"]: j for j in extra["jobs"]}}
    for job in extra["jobs"]:
        clone = dict(job)
        claimed = clone.pop("job_sha256")
        check(canonical(clone) == claimed and job["level"] in (4, 5), job["job_id"] + ": allowed additional identity")
    check(len(jobs) == len(manifest["jobs"]), "unique job IDs")
    check(manifest["coefficient"] == 80.61638494312905, "frozen coefficient")
    identity(d / "input.mat", manifest["input_sha256"], "MAT input")
    for path, expected in manifest["source_sha256"].items():
        identity(path, expected, "frozen preparation source")
    identity(analysis / "analysis_source.py", summary["analysis_sha256"], "analysis source snapshot")
    ah_sources = read(d / "ah_source_sha256.json")
    for path, expected in ah_sources.items():
        identity(d / "ah_source" / Path(path).name, expected, "AH execution source snapshot")
    root = Path(__file__).resolve().parents[1]
    before = read(run/"source_before_manifest.json")
    expected_sources = {str((root/item["path"]).resolve()): item["sha256"] for item in before["records"]}
    expected_sources.update({str(Path(p).resolve()): value for p, value in ah_sources.items()})
    for path, value in manifest["source_sha256"].items():
        expected_sources.setdefault(str((root/"tests"/Path(path).name).resolve()), value)
    expected_sources[str((root/"tests/_a4_c1m_p3_compare.py").resolve())] = read(parent_analysis/"summary.json")["analysis_sha256"]
    entrypoint = root/"tests/_a4_c1m_p3_compare.py"
    if parent_analysis != analysis:
        entrypoint = root/"tests/_a4_c1m_p3_refinement_acceptance.py"
        expected_sources[str(entrypoint.resolve())] = summary["analysis_sha256"]
    p2_manifest = read(run.parent/"p2_20261008_v1/p2_manifest.json")
    historical_ode = p2_manifest["sources"]["historical_P1_ODE"]
    expected_sources[str(Path(historical_ode["path"]).resolve())] = historical_ode["sha256"]
    entries = [entrypoint, Path(historical_ode["path"])]
    for source in local_import_closure(entries, root):
        expected = expected_sources.get(str(source))
        check(expected is not None, "analysis dependency frozen identity: "+str(source))
        if expected:
            identity(source, expected, "analysis transitive dependency")
    mex_paths = set()
    for path in d.glob("*_pharlap_summary.json"):
        execution = read(path)
        if execution.get("mex_path"):
            mex_paths.add(execution["mex_path"])
    check(len(mex_paths) == 1, "one actual PH MEX selected across runs")
    actual_mex_sha = sha(next(iter(mex_paths))) if len(mex_paths) == 1 else None
    for job in jobs.values():
        independent = dict(job)
        claimed = independent.pop("job_sha256")
        check(canonical(independent) == claimed, job["job_id"] + ": canonical identity")
        identity(job["grid_mat"], job["grid_sha256"], job["job_id"] + ": grid")
        check(tuple(job["tol"]) == PH_PROFILES[job["level"]-1], job["job_id"] + ": legal PH profile")
        if job["role"] != "R5_fd":
            path = d / "pharlap_raw" / (job["job_id"]+".mat")
            try:
                raw = loadmat(path, squeeze_me=True, struct_as_record=False, mat_dtype=True)
                check(str(raw["job_sha256"]) == job["job_sha256"], job["job_id"] + ": PH immutable job")
                check(str(raw["grid_sha256"]) == job["grid_sha256"], job["job_id"] + ": PH immutable grid")
                check(str(raw["mex_sha256"]) == actual_mex_sha, job["job_id"] + ": actual MEX identity")
                check(str(raw["exporter_sha256"]) == job["exporter_sha256"], job["job_id"] + ": actual exporter identity")
                identities.append({"path": str(path.resolve()), "actual_sha256": sha(path),
                                   "embedded_job_sha256": str(raw["job_sha256"]),
                                   "passed": str(raw["job_sha256"]) == job["job_sha256"]
                                   and str(raw["grid_sha256"]) == job["grid_sha256"]
                                   and str(raw["mex_sha256"]) == actual_mex_sha
                                   and str(raw["exporter_sha256"]) == job["exporter_sha256"]})
            except (OSError, KeyError) as exc:
                errors.append(job["job_id"] + ": missing PH raw evidence: " + str(exc))
    contract_path = run/"analysis/p3_contract_v3/report.json"
    if contract_path.is_file():
        contract = read(contract_path)
        grid_by_hash = {j["grid_sha256"]: j["grid_mat"] for j in jobs.values()}
        check(set(grid_by_hash) == {x["grid_sha256"] for x in contract["grid_readbacks"]},
              "readback covers every unique input grid")
        for item in contract["grid_readbacks"]:
            identity(item["path"], item["sha256"], "PH input readback")
            original_grid = loadmat(grid_by_hash[item["grid_sha256"]], squeeze_me=True)
            readback = loadmat(item["path"], squeeze_me=True)
            for field, claimed in item["fields_equal"].items():
                actual = np.array_equal(original_grid[field], readback[field])
                check(claimed is bool(actual) and actual, "PH readback exact array: "+field)
        check(contract["grid_readbacks_passed"] is True, "grid readback gate")
    check(summary.get("global_passed") is False, "global acceptance remains false before P4-P6")
    formal = []
    for row in summary["cases"]:
        label = row["case_id"]
        start_errors = len(errors)
        check(type(row["strict_passed"]) is bool, label + ": real JSON boolean")
        case_jobs = [jobs[k] for k in row["jobs"]]
        check(sorted(j["level"] for j in case_jobs) == [1, 2, 3], label + ": three integration levels")
        check(all(j["scope"] == row["scope"] and j["role"] == row["role"] for j in case_jobs), label + ": scope/role identity")
        case_jobs.sort(key=lambda j: j["level"])
        traces = {"ah": [], "ph": []}
        load_failures = []
        for index, job in enumerate(case_jobs):
            name = job["job_id"]
            try:
                ah_id = row.get("levels", {}).get("ah", [{}]*3)[index].get("job_id", name)
                ah_job = all_jobs[ah_id]
                ah_directory = d/"refinements_v1" if ah_job["level"] > 3 else d
                meta_path = ah_directory / "ah" / (ah_id + ".json")
                meta = read(meta_path)
                check(meta["job_sha256"] == ah_job["job_sha256"], ah_id + ": AH job identity")
                check(meta["source_sha256"] == ah_sources, name + ": AH execution sources")
                if meta["status"] != "completed":
                    raise ValueError("AH failed: " + meta.get("error", "unknown"))
                raw_path = meta_path.with_suffix(".npz")
                identity(raw_path, meta["raw_sha256"], name + ": AH raw")
                ah = ah_trace(raw_path)
                ah["termination"] = meta["termination"]
                ah["meta"] = meta
                validate_trace(ah, row["scope"])
                traces["ah"].append(ah)
                ph_path = d / "pharlap_raw" / (job.get("ph_job_id", name) + ".mat")
                ph, raw = ph_trace(ph_path)
                check(str(raw["job_sha256"]) == job["job_sha256"], name + ": PH job identity")
                check(str(raw["grid_sha256"]) == job["grid_sha256"], name + ": PH grid identity")
                check(str(raw["exporter_sha256"]) == job["exporter_sha256"], name + ": PH exporter identity")
                validate_trace(ph, row["scope"])
                traces["ph"].append(ph)
            except (ValueError, OSError, KeyError, AttributeError) as exc:
                load_failures.append(name + ": " + str(exc))
        if load_failures:
            check(row["strict_passed"] is False and "error" in row, label + ": failed raw cannot pass")
            rows.append({"case_id": label, "scientific_passed": False,
                         "raw_failures": load_failures, "verification_passed": len(errors) == start_errors})
        elif row["scope"] != "first_100_km" and "actual" not in row:
            check(row["strict_passed"] is False, label + ": unresolved full path cannot pass")
            rows.append({"case_id": label, "scientific_passed": False,
                         "unresolved": row.get("error"), "verification_passed": len(errors) == start_errors})
        else:
            try:
                for side in ("ah", "ph"):
                    actual_levels = [x["level"] for x in row["levels"][side]]
                    check(actual_levels == [1, 2, 3] or (side == "ah" and actual_levels == [3, 4, 5]),
                          label + ": reported levels " + side)
                    for job, trace, reported in zip(case_jobs, traces[side], row["levels"][side]):
                        identifier = reported.get("job_id", job["job_id"])
                        reported_job = all_jobs[identifier]
                        basedir = d/"refinements_v1" if side == "ah" and reported_job["level"] > 3 else d
                        expected_path = basedir / ("ah" if side == "ah" else "pharlap_raw") / (
                            identifier + (".npz" if side == "ah" else ".mat"))
                        check(Path(reported["raw_path"]).resolve() == expected_path.resolve(), label + ": reported raw path " + side)
                        identity(expected_path, reported["raw_sha256"], label + ": analysis raw binding " + side)
                        check(reported["termination"] == trace["termination"], label + ": native termination " + side)
                        check(reported["native_s_end_km"] == float(trace["s"][-1])
                              and reported["native_group_end_km"] == float(trace["group"][-1]), label + ": native endpoints " + side)
                formal_step = row.get("formal_query_step_km", .25)
                actual = metrics(traces["ah"][-1], traces["ph"][-1], row["scope"], formal_step)
                close_map(actual, row["actual"], label + ": actual recomputed")
                ctrl = {}
                for side in ("ah", "ph"):
                    t = traces[side]
                    ctrl[side] = {"delta12": metrics(t[0], t[1], row["scope"], formal_step),
                                  "delta23": metrics(t[1], t[2], row["scope"], formal_step)}
                    for key in ("delta12", "delta23"):
                        close_map(ctrl[side][key], row["integration"][side][key], label + ": " + side + " " + key)
                steps = row["resampling"]["steps_km"]
                check(steps in ([1., .5, .25], [.25, .125, .0625]), label + ": legal sampling levels")
                sampled = [metrics(traces["ah"][-1], traces["ph"][-1], row["scope"], h) for h in steps]
                for x, y in zip(sampled, row["resampling"]["metrics"]):
                    close_map(x, y, label + ": sampling metrics")
                sampling = {k: abs(sampled[2][k]-sampled[1][k]) for k in LIMITS}
                sampling12 = {k: abs(sampled[1][k]-sampled[0][k]) for k in LIMITS}
                close_map(sampling, row["resampling"]["final_bound"], label + ": sampling final bound")
                close_map(sampling12, row["resampling"]["delta12"], label + ": sampling first change")
                reconstruction = row["reconstruction"]["bound"]
                for side in ("ah", "ph"):
                    item = row["reconstruction"]["sides"][side]
                    proof = item["proof"]
                    check(item["bound_established"] is True and proof["node_error_in_reconstruction_budget"] is False
                          and proof["uncertainty_multiplier"] == 1., label + ": uninflated interpolation-only bound " + side)
                    position_bound = (proof["actual_reference_interpolation_max_m"]
                                      + proof["fill_bound_m"] + proof["quad_uncertainty_m"])
                    endpoint_bound = proof["actual_reference_interpolation_endpoint_m"] + proof["quad_uncertainty_m"]
                    group_bound = ((proof["actual_reference_group_interpolation_endpoint_km"]
                                    + proof["group_quad_uncertainty_km"])
                                   / abs(float(curves(traces[side][-1])[1](100.))) if side == "ph" else 0.)
                    close_map({"rms_m": position_bound, "max_m": position_bound,
                               "landing_m": endpoint_bound, "group_rel": group_bound},
                              item["bounds"], label + ": reconstruction equation " + side)
                summed_bound = {k: sum(row["reconstruction"]["sides"][side]["bounds"][k]
                                      for side in ("ah", "ph")) for k in LIMITS}
                close_map(summed_bound, reconstruction, label + ": reconstruction side sum")
                quad_path = analysis / (label+"_quadrature.json")
                quad = read(quad_path if quad_path.is_file() else parent_analysis/quad_path.name)
                if "field_reduction" in row["independent_reference"]:
                    grid = loadmat(case_jobs[0]["grid_mat"], squeeze_me=True)
                    heights = np.asarray(grid["height_km"]).reshape(-1)
                    density = np.asarray(grid["ne_profile_m3"]).reshape(-1)
                    spline = CubicHermiteSpline(heights, density, np.gradient(density, heights, edge_order=2))
                    coeff = spline.c
                    exact_affine = bool(np.all(coeff[:2] == 0.) and np.ptp(coeff[2]) == 0.
                        and np.array_equal(coeff[3], coeff[3, 0]+coeff[2, 0]*(heights[:-1]-heights[0])))
                    reduction = row["independent_reference"]["field_reduction"]
                    check(reduction["exact_global_affine"] is exact_affine,
                          label + ": exact affine reduction from original grid")
                    check(reduction["only_reference_knots_changed"] is True,
                          label + ": reduction preserves production grid")
                reference_control = row["independent_reference"]["convergence"]
                lower_warnings = []
                for level in (3, 4):
                    lower_path = analysis/(label+f"_quadrature_l{level}.json")
                    if not lower_path.is_file():
                        lower_path = parent_analysis/lower_path.name
                    if lower_path.is_file():
                        lower_warnings.extend(read(lower_path)["warnings"])
                reference_valid = all(reference_control["delta23"][k] <= reference_control["delta12"][k]
                    and reference_control["delta23"][k] <= .2*LIMITS[k] for k in LIMITS) and not quad["warnings"] and not lower_warnings
                check(row["reconstruction"]["validated"] is bool(reference_valid), label + ": reference validation gate")
                effective = {k: sampling[k]+reconstruction[k] for k in LIMITS}
                budget = {k: 3.*(ctrl["ah"]["delta23"][k]+ctrl["ph"]["delta23"][k]+effective[k]) for k in LIMITS}
                close_map(budget, row["budget"], label + ": independent budget sum")
                initial = all(np.linalg.norm(a["xyz"][0]-p["xyz"][0])*1000. <= 1e-3
                              and np.linalg.norm(a["p"][0]-p["p"][0]) <= 1e-10
                              for a, p in zip(traces["ah"], traces["ph"]))
                if row["scope"] == "first_100_km":
                    scope = all(t["termination"] == "fixed_geometric_path_completed"
                                and abs(t["s"][-1]-100.) <= 1e-8 for t in traces["ah"])
                else:
                    scope = all(a["termination"] == "ground" and p["termination"] == 1
                                for a, p in zip(traces["ah"], traces["ph"]))
                criteria = {
                    "thresholds": all(actual[k] <= LIMITS[k] for k in LIMITS),
                    "budget": all(actual[k] <= budget[k] for k in LIMITS),
                    "integration": all(c["delta23"][k] <= c["delta12"][k] and c["delta23"][k] <= .2*LIMITS[k]
                                       for c in ctrl.values() for k in LIMITS),
                    "sampling": all(sampling[k] <= .2*LIMITS[k] and sampling[k] <= sampling12[k] for k in LIMITS),
                    "reconstruction": row["reconstruction"]["validated"] is True
                        and all(np.isfinite(reconstruction[k]) and 0 <= reconstruction[k] <= .2*LIMITS[k] for k in LIMITS),
                    "initial_state": initial, "scope": scope,
                }
                recomputed_pass = bool(all(criteria.values()))
                check(row["strict_passed"] is recomputed_pass, label + ": strict gate recomputed")
                for key, value in criteria.items():
                    check(row["checks"].get(key) is bool(value), label + ": gate " + key)
                rows.append({"case_id": label, "scientific_passed": recomputed_pass,
                             "criteria": {k: bool(v) for k, v in criteria.items()},
                             "actual": actual, "budget": budget,
                             "verification_passed": len(errors) == start_errors})
            except (ValueError, KeyError, IndexError) as exc:
                errors.append(label + ": numerical verification failed: " + str(exc))
        if row["role"] in ("R1", "R2", "R3", "R4"):
            formal.append(row)
    check(len(formal) == 24 and summary["case_count"] == 24, "24 required formal cases")
    check(summary["formal_passed_count"] == sum(r["strict_passed"] for r in formal), "formal pass count")
    check(summary["p3_scope_passed"] is (len(formal) == 24 and all(r["strict_passed"] for r in formal)), "P3 gate")
    r5 = verify_r5(run)
    check(not summary["p3_scope_passed"] or r5["r5_passed"], "P3 completion includes R5 independent controls")
    return {"schema": "c1m-p3-independent-verification-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
            "verifier_sha256": sha(__file__), "summary_sha256": sha(analysis / "summary.json"),
            "self_checks": self_check(), "identities": identities, "cases": rows, "r5": r5,
            "errors": errors, "error_count": len(errors), "verification_passed": not errors,
            "scientific_scope_passed": summary["p3_scope_passed"], "global_passed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--analysis", type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    analysis = args.analysis.resolve() if args.analysis else run / "analysis/p3_v1"
    result = verify(run, analysis)
    version = 1
    while (analysis / f"verification_v{version}").exists():
        version += 1
    out = analysis / f"verification_v{version}"
    out.mkdir()
    (out / "verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (out / "verifier_source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"output": str(out), "error_count": result["error_count"], "verification_passed": result["verification_passed"]}))
    return 0 if result["verification_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
