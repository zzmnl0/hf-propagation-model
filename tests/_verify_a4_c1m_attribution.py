"""Read-only C1M freeze and uniform-medium evidence checks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests import _run_a4_c1m_attribution as experiment


def self_check() -> None:
    with tempfile.TemporaryDirectory(prefix='c1m_json_') as folder:
        output = Path(folder)/'values.json'
        experiment.write_json(output, {'false':np.bool_(False),'true':np.bool_(True)})
        values = json.loads(output.read_text(encoding='utf-8'))
        assert values['false'] is False and values['true'] is True
    s = np.array([0., 25., 50., 100.])
    points = np.column_stack((2*s+5, -3*s, s/2))
    target = np.linspace(0, 100, 401)
    expected = np.column_stack((2*target+5, -3*target, target/2))
    assert np.allclose(experiment.resample(s, points, target), expected, atol=1e-12)
    for bad in (np.array([0., 25., 50., 99.]), np.array([0., 40., 30., 100.])):
        try:
            experiment.resample(bad, points, target)
        except ValueError:
            pass
        else:
            raise AssertionError("truncated/nonmonotone path was accepted")
    # A constant uniform medium has an arbitrary conserved initial |p|.
    # Computing K from that |p| recovers the supplied initialization, not an
    # independently identified coefficient of the propagation implementation.
    for k2 in (0., 80.6404, 80.6163849431291):
        n = np.sqrt(1-k2*1e11/1e14)
        assert abs((1-n*n)*1e14/1e11-k2) < 1e-10
    print("C1M analysis self-check: PASS")


def verify(run: Path) -> list[str]:
    failures = []

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    baseline = json.loads((run / 'a0_manifest.json').read_text(encoding='utf-8'))
    # Preserve exact A0 tool snapshots after checker-path/JSON-type fixes.
    # All production physics sources must still match their live files.
    sources = experiment.source_files()
    for name, expected in baseline['source_sha256'].items():
        source = sources[name]
        if name in {'verifier', 'runner', 'matlab_source'} and experiment.sha(source) != expected:
            source = run / (f'{name}_frozen_a0.m' if name == 'matlab_source'
                            else f'{name}_frozen_a0.py')
        check(source.is_file() and experiment.sha(source) == expected,
              f'frozen source mismatch: {name}')
    check(experiment.sha(run/'input.mat') == baseline['input_mat_sha256'], 'input MAT drift')
    summary = json.loads((run / "pharlap_summary.json").read_text(encoding="utf-8"))
    check(summary["status"] == "completed_uniform_export", "MATLAB export incomplete")
    check(summary.get("coefficient_identifiable") is False, "unsupported coefficient identification")
    check(summary.get("calibration_k2") is None, "uniform calibration must not claim intrinsic K")
    check(abs(summary["launch_k2"]-experiment.K_FP**2) < 1e-12, "common launch coefficient mismatch")
    data = loadmat(run / "input.mat", simplify_cells=True)
    back = loadmat(run / "pharlap_raw" / "input_readback.mat", simplify_cells=True)
    for key, value in data.items():
        if not key.startswith("__"):
            check(np.array_equal(value, back[key]), f"MATLAB array readback mismatch: {key}")

    records = summary["records"]
    check(len(records) == 72, "expected 24 cases x 3 PH levels")
    expected_ids = {f"case_n{ni}_f{fi}_e{ei}_l{li}" for ni in range(1,5)
                    for fi in range(1,4) for ei in range(1,3) for li in range(1,4)}
    check({r["case_id"] for r in records} == expected_ids, "PH case identities mismatch")
    for record in records:
        raw = loadmat(record["raw_mat"], simplify_cells=True)
        state, path = raw["ray_state_vec"], raw["ray_path_data"]
        arc = np.asarray(state["geometrical_path"], float)/1000
        check(arc[0] == 0 and arc[-1] >= 100 and np.all(np.diff(arc)>0),
              f"invalid common arc: {record['case_id']}")
        n = np.asarray(path["refractive_index"],float)
        ng = np.asarray(path["group_refractive_index"],float)
        density = np.asarray(path["electron_density"],float)*1e6
        mask = arc <= 100
        check(np.max(np.abs(density[mask]-record["density_m3"])) <= max(1e-5,record["density_m3"]*1e-12),
              f"density mismatch: {record['case_id']}")
        check(np.max(np.abs(n[mask]-record["initial_index"])) < 1e-10,
              f"uniform phase index mismatch: {record['case_id']}")
        check(np.max(np.abs(ng[mask]*n[mask]-1)) < 1e-10,
              f"uniform group index mismatch: {record['case_id']}")
        p = np.column_stack([np.asarray(state[k],float) for k in ('dir_x','dir_y','dir_z')])
        ei = int(record["case_id"].split('_')[3][1:])-1
        target = data["launch_unit_ecef"][ei]*record["initial_index"]
        check(np.linalg.norm(p[0]-target) < 1e-12, f"initial wave normal mismatch: {record['case_id']}")

    controls = summary["initialization_records"]
    check(len(controls) == 9, "expected nine initialization controls")
    check({r["OX_mode"] for r in controls} == {-1,0,1}, "missing B=0 mode control")
    for record in controls:
        raw = loadmat(record["raw_mat"], simplify_cells=True)
        n = np.asarray(raw["ray_path_data"]["refractive_index"],float)
        check(np.max(np.abs(n-record["initial_index"])) < 1e-10,
              f"initialization control changed: {record['case_id']}")

    # Audit repeats from measured arrays, not legacy JSON string booleans.
    for record in summary['cache_records']:
        ni = 2 if record['density_m3']==1e10 else 4
        reference = run/'pharlap_raw'/f'case_n{ni}_f2_e1_l3.mat'
        first = loadmat(reference,simplify_cells=True)['ray_state_vec']
        repeat = loadmat(record['raw_mat'],simplify_cells=True)['ray_state_vec']
        for field in ('pos_x','pos_y','pos_z','group_path','geometrical_path'):
            check(np.array_equal(first[field],repeat[field]),f'repeat drift: {record["case_id"]} {field}')

    analysis = run / "a1"
    report = json.loads((analysis / "a1_report.json").read_text(encoding="utf-8"))
    check(len(report["cases"]) == 24, "missing A1 case reports")
    for case in report['cases']:
        levels = [np.load(analysis/f'case_{case["case_id"]}_l{li}_samples.npz')
                  for li in range(case['levels']-2,case['levels']+1)]
        for side in ('ah','ph'):
            increments = []
            for low,high in zip(levels,levels[1:]):
                delta = np.linalg.norm(high[side+'_xyz_km']-low[side+'_xyz_km'],axis=1)*1000
                key = 'ah_group_native_100km_km' if side=='ah' else 'ph_group_km'
                first,last = (float(np.asarray(z[key]).reshape(-1)[-1]) for z in (low,high))
                increments.append({'rms_m':float(np.sqrt(np.mean(delta**2))),
                                   'max_m':float(delta.max()),'group_rel':abs(last-first)/abs(last)})
            for key,floor in experiment.ROUND_FLOOR.items():
                first,last = [row[key] for row in increments]
                check(last <= first or max(first,last) <= floor,
                      f'neither decreasing nor roundoff-limited: {case["case_id"]} {side} {key}')
        for level in levels:
            level.close()
    check(report["a0_manifest_sha256"] == experiment.sha(run / "a0_manifest.json"), "A0 manifest drift")
    check(report["pharlap_summary_sha256"] == experiment.sha(run / "pharlap_summary.json"), "PH summary drift")
    manifest_path = analysis / "a1_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    # Schema deliberately mirrors the small experiment's named output map.
    artifacts = manifest.get("artifacts", manifest.get("sha256", {}))
    if not artifacts:
        failures.append("A1 manifest has no artifact hashes")
    else:
        for name, expected in artifacts.items():
            check(experiment.sha(analysis/name) == expected, f"artifact drift: {name}")
    for name, expected in report["pharlap_raw_sha256"].items():
        raw_path = run / "pharlap_raw" / (name if name.endswith('.mat') else name+'.mat')
        check(experiment.sha(raw_path) == expected, f"raw drift: {name}")
    check(baseline["status"] == "frozen_a0", "A0 freeze not complete")
    return failures


def _a2_jobs(path: Path) -> list[dict]:
    jobs = loadmat(path, simplify_cells=True)["a2_jobs"]
    return [jobs] if isinstance(jobs, dict) else list(jobs)


def _verify_a2_identified_holdout(
    root: Path, base_jobs: list[dict], freeze: dict
) -> tuple[list[dict], list[str]]:
    """Check the 8/15 MHz frozen-K2 control without running either tracer."""
    failures: list[str] = []

    def check(ok: bool, message: str) -> None:
        if not ok:
            failures.append(message)

    directory = root / "identified_holdout"
    if not directory.is_dir():
        return [], ["A2 identified 8/15 MHz holdout directory missing"]
    manifest = json.loads((directory / "input_manifest.json").read_text(
        encoding="utf-8"))
    jobs = manifest.get("jobs", [])
    selected = {
        str(job["job_id"]): job for job in base_jobs
        if job["stage"] == "affine"
        and float(job["frequency_MHz"]) in (8.0, 15.0)
        and str(job["case_id"]).endswith("d0p25")
    }
    ids = [str(job["job_id"]) for job in jobs]
    check(manifest.get("schema") == "a4-c1m-a2-identified-holdout-input-v1",
          "A2 identified holdout manifest schema mismatch")
    check(len(jobs) == len(set(ids)) == len(selected) == 48
          and set(ids) == set(selected),
          "A2 identified holdout must contain the 48 frozen 8/15 MHz jobs")
    check(experiment.sha(directory / "input.mat") ==
          manifest.get("input_mat_sha256"),
          "A2 identified holdout input MAT SHA mismatch")
    check(experiment.sha(root / "input_manifest.json") ==
          manifest.get("base_input_manifest_sha256"),
          "A2 identified holdout base-manifest SHA mismatch")
    check(experiment.sha(root / "identified" / "coefficient_freeze.json") ==
          manifest.get("coefficient_freeze_sha256"),
          "A2 identified holdout coefficient-freeze SHA mismatch")
    coefficient = float(freeze["identified_K2"])
    check(float(manifest.get("identified_K2", float("nan"))) == coefficient,
          "A2 identified holdout coefficient differs from training freeze")
    mat_jobs = _a2_jobs(directory / "input.mat")
    mat_by_id = {str(job["job_id"]): job for job in mat_jobs}
    check(len(mat_jobs) == len(mat_by_id) == 48 and
          set(mat_by_id) == set(ids),
          "A2 identified holdout MAT job identities mismatch")

    summary = json.loads((directory / "pharlap_summary.json").read_text(
        encoding="utf-8"))
    records = summary.get("records", [])
    record_by_id = {str(row["job_id"]): row for row in records}
    check(summary.get("schema_version") == "c1m-a2-export-v1" and
          summary.get("status") == "completed_a2_export" and
          summary.get("failures") == 0 and
          len(records) == len(record_by_id) == 48 and
          set(record_by_id) == set(ids),
          "A2 identified holdout PH export incomplete")
    execution_path = directory / "execution_record.json"
    if execution_path.is_file():
        execution = json.loads(execution_path.read_text(encoding="utf-8"))
        exporter = (directory / execution.get("exporter_source", "")).resolve()
        check(execution.get("matlab_exit_code") == 0 and
              execution.get("record_count") == 48 and
              execution.get("failed_records") == 0 and
              exporter == (root / "sources" /
                           "a4_c1m_uniform_export.m").resolve(),
              "A2 identified holdout MATLAB execution/source mismatch")
    else:
        failures.append("A2 identified holdout MATLAB execution record missing")

    grids_checked: set[str] = set()
    for job in jobs:
        name = str(job["job_id"])
        try:
            base = selected[name]
            mat_job = mat_by_id[name]
            record = record_by_id[name]
            check(all(job[key] == base[key] for key in (
                "job_id", "case_id", "stage", "grid_mat", "grid_sha256",
                "frequency_MHz", "elevation_deg", "origin_height_km",
                "explicit_state", "level")) and
                  np.array_equal(np.asarray(job["tol"], float),
                                 np.asarray(base["tol"], float)),
                f"A2 identified holdout base-job identity drift: {name}")
            check(bool(job["explicit_state"]) and
                  float(job["launch_k2"]) == coefficient and
                  float(job["frequency_MHz"]) in (8.0, 15.0) and
                  int(job["level"]) in (1, 2, 3),
                  f"A2 identified holdout initial-state settings wrong: {name}")
            unsigned = {key: value for key, value in job.items()
                        if key != "job_sha256"}
            digest = hashlib.sha256(json.dumps(
                unsigned, sort_keys=True, separators=(",", ":"),
                allow_nan=False).encode("utf-8")).hexdigest()
            check(job["job_sha256"] == digest and
                  mat_job["job_sha256"] == digest and
                  mat_job["grid_sha256"] == job["grid_sha256"] and
                  float(mat_job["launch_k2"]) == coefficient and
                  float(mat_job["frequency_MHz"]) ==
                  float(job["frequency_MHz"]) and
                  float(mat_job["elevation_deg"]) ==
                  float(job["elevation_deg"]) and
                  float(mat_job["origin_height_km"]) == 100.0 and
                  bool(mat_job["explicit_state"]) and
                  np.array_equal(np.asarray(mat_job["tol"], float).reshape(-1),
                                 np.asarray(job["tol"], float)),
                  f"A2 identified holdout job/MAT SHA or tolerance drift: {name}")
            grid_path = Path(job["grid_mat"])
            grid_hash = str(job["grid_sha256"])
            if grid_hash not in grids_checked:
                check(experiment.sha(grid_path) == grid_hash,
                      f"A2 identified holdout grid SHA drift: {name}")
                readback_path = directory / "pharlap_raw" / (
                    f"grid_readback_{grid_hash[:16]}.mat")
                if readback_path.is_file():
                    original = loadmat(grid_path, simplify_cells=True)
                    readback = loadmat(readback_path, simplify_cells=True)
                    check(all(key in readback and
                              np.array_equal(value, readback[key])
                              for key, value in original.items()
                              if not key.startswith("__")),
                          f"A2 identified holdout grid readback drift: {name}")
                else:
                    failures.append(
                        f"A2 identified holdout grid readback missing: {name}")
                grids_checked.add(grid_hash)
            raw_path = directory / "pharlap_raw" / f"{name}.mat"
            check(Path(record["raw_mat"]).resolve() == raw_path.resolve() and
                  record["job_sha256"] == digest and
                  record["grid_sha256"] == grid_hash and
                  not record.get("error_identifier") and
                  not record.get("error_message"),
                  f"A2 identified holdout PH record mismatch: {name}")
            raw = loadmat(raw_path, simplify_cells=True)
            check(str(raw["job_sha256"]) == digest and
                  str(raw["grid_sha256"]) == grid_hash and
                  str(raw["job"]["job_id"]) == name and
                  float(raw["job"]["launch_k2"]) == coefficient,
                  f"A2 identified holdout raw identity drift: {name}")
            state, path, initial = (raw[key] for key in (
                "ray_state_vec", "ray_path_data", "initial"))
            p2 = sum(np.asarray(state[key], float).reshape(-1) ** 2
                     for key in ("dir_x", "dir_y", "dir_z"))
            ne = np.asarray(path["electron_density"], float).reshape(-1) * 1e6
            arc = np.asarray(state["geometrical_path"], float).reshape(-1)
            check(len(p2) == len(ne) == len(arc) ==
                  int(record["point_count"]) >= 2 and
                  np.all(np.diff(arc) > 0) and
                  arc[-1] - arc[0] >= 100000.0,
                  f"A2 identified holdout PH path incomplete: {name}")
            check(all(np.all(np.asarray(path[key], float) == 0)
                      for key in ("geomag_x", "geomag_y", "geomag_z",
                                  "collision_frequency")),
                  f"A2 identified holdout nonzero B/collision: {name}")
            grid = loadmat(grid_path, simplify_cells=True)
            h = np.asarray(grid["height_km"], float).reshape(-1)
            profile = np.asarray(grid["ne_profile_m3"], float).reshape(-1)
            ne0 = float(np.interp(float(job["origin_height_km"]), h, profile))
            n2 = 1.0 - coefficient * ne0 / (float(job["frequency_MHz"])*1e6)**2
            initial_wave = np.array([float(initial[key]) for key in
                                     ("dir_x", "dir_y", "dir_z")])
            first_wave = np.array([
                np.asarray(state[key], float).reshape(-1)[0] for key in
                ("dir_x", "dir_y", "dir_z")])
            check(n2 > 0 and abs(initial_wave @ initial_wave - n2) < 1e-10 and
                  np.linalg.norm(first_wave - initial_wave) < 1e-10,
                  f"A2 identified holdout frozen-K2 launch mismatch: {name}")
            delta_ne = ne - ne[0]
            mask = np.abs(delta_ne) >= 1e8
            check(np.count_nonzero(mask) >= 2,
                  f"A2 identified holdout coefficient samples absent: {name}")
            if np.count_nonzero(mask) >= 2:
                effective = -(float(job["frequency_MHz"])*1e6)**2 * (
                    p2[mask] - p2[0]) / delta_ne[mask]
                check(abs(float(np.median(effective)) - coefficient) < 1e-8,
                      f"A2 identified holdout K2 did not reproduce: {name}")
        except (KeyError, ValueError, OSError, TypeError) as exc:
            failures.append(f"A2 identified holdout unreadable: {name}: {exc}")
    return jobs, failures


def _a2_sample(s: np.ndarray, values: np.ndarray, target: np.ndarray,
               *, linear: bool = False) -> np.ndarray:
    s = np.asarray(s, float).reshape(-1)
    values = np.asarray(values, float)
    target = np.asarray(target, float).reshape(-1)
    if (len(s) < 2 or np.any(np.diff(s) <= 0) or
            target[0] < s[0] - 1e-8 or target[-1] > s[-1] + 1e-8):
        raise ValueError("raw ray does not span a monotone comparison arc")
    target = np.clip(target, s[0], s[-1])
    if not linear:
        return np.asarray(PchipInterpolator(s, values, axis=0)(target))
    if values.ndim == 1:
        return np.interp(target, s, values)
    return np.column_stack([
        np.interp(target, s, values[:, col])
        for col in range(values.shape[1])])


def _a2_position_metrics(first: np.ndarray, second: np.ndarray) -> dict:
    distance_m = np.linalg.norm(first - second, axis=1) * 1000.0
    return {"rms_m": float(np.sqrt(np.mean(distance_m ** 2))),
            "max_m": float(np.max(distance_m))}


def _a2_close(observed: object, expected: float, *, atol: float = 1e-8) -> bool:
    return (type(observed) in (float, int) and
            np.isfinite(observed) and
            abs(float(observed) - expected) <= atol +
            1e-10 * abs(expected))


def _a2_native_trace(dataset_dir: Path, ah_dir: Path, job: dict) -> dict:
    name = str(job["job_id"])
    raw = loadmat(dataset_dir / "pharlap_raw" / f"{name}.mat",
                  simplify_cells=True)
    state = raw["ray_state_vec"]
    ph_xyz = np.column_stack([
        np.asarray(state[key], float).reshape(-1)
        for key in ("pos_x", "pos_y", "pos_z")]) / 1000.0
    ph_s = np.asarray(state["geometrical_path"], float).reshape(-1) / 1000.0
    ph_group = np.asarray(state["group_path"], float).reshape(-1) / 1000.0
    ph_s = ph_s - ph_s[0]
    ph_group = ph_group - ph_group[0]
    ph_p2 = sum(np.asarray(state[key], float).reshape(-1) ** 2
                for key in ("dir_x", "dir_y", "dir_z"))
    ah_path = ah_dir / f"{name}.json"
    ah = json.loads(ah_path.read_text(encoding="utf-8"))
    trace = {"ph_s": ph_s, "ph_xyz": ph_xyz, "ph_group": ph_group,
             "ph_p2": ph_p2, "ah_record": ah}
    if ah["status"] == "completed":
        with np.load(ah_dir / f"{name}.npz") as data:
            trace["ah_s"] = np.asarray(data["geom_path_km"], float)
            trace["ah_xyz"] = np.asarray(data["xyz_km"], float)
    return trace


def _a2_gate_from_traces(levels: list[dict], traces: list[dict],
                         rows: list[dict]) -> dict:
    fixed = str(levels[0]["stage"]) in ("internal", "affine")
    limits = {"rms_m": 10.0, "max_m": 50.0, "landing_m": 50.0,
              "group_rel": 1e-5}
    floor = {"rms_m": 1e-7, "max_m": 1e-7, "landing_m": 1e-7,
             "group_rel": 1e-13}

    def ph_endpoint(trace: dict) -> tuple[np.ndarray, float]:
        if fixed:
            xyz = _a2_sample(trace["ph_s"], trace["ph_xyz"],
                             np.array([100.0]))[0]
            group = float(_a2_sample(trace["ph_s"], trace["ph_group"],
                                     np.array([100.0]))[0])
            return xyz, group
        return trace["ph_xyz"][-1], float(trace["ph_group"][-1])

    def increment(first: dict, second: dict) -> dict:
        result = {}
        for side in ("ph", "ah"):
            s1, xyz1 = first[f"{side}_s"], first[f"{side}_xyz"]
            s2, xyz2 = second[f"{side}_s"], second[f"{side}_xyz"]
            end = 100.0 if fixed else min(float(s1[-1]), float(s2[-1]))
            target = np.linspace(0.0, end, 201 if fixed else 501)
            result[side] = _a2_position_metrics(
                _a2_sample(s1, xyz1, target),
                _a2_sample(s2, xyz2, target))
        ph1_xyz, ph1_group = ph_endpoint(first)
        ph2_xyz, ph2_group = ph_endpoint(second)
        result["ph"]["landing_m"] = float(
            np.linalg.norm(ph1_xyz - ph2_xyz) * 1000.0)
        result["ph"]["group_rel"] = abs(ph1_group - ph2_group) / ph2_group
        result["ah"]["landing_m"] = float(np.linalg.norm(
            first["ah_xyz"][-1] - second["ah_xyz"][-1]) * 1000.0)
        result["ah"]["group_rel"] = abs(
            float(first["ah_record"]["group_path_km_native"]) -
            float(second["ah_record"]["group_path_km_native"])) / ph2_group
        return result

    d12 = increment(traces[0], traces[1])
    d23 = increment(traces[1], traces[2])
    fine = traces[2]
    end = (100.0 if fixed else
           min(float(fine["ph_s"][-1]), float(fine["ah_s"][-1])))

    def positions(step: float, linear: bool = False) -> tuple[np.ndarray, np.ndarray]:
        target = np.linspace(0.0, end, int(np.ceil(end / step)) + 1)
        return (
            _a2_sample(fine["ph_s"], fine["ph_xyz"], target, linear=linear),
            _a2_sample(fine["ah_s"], fine["ah_xyz"], target, linear=linear))

    ph05, ah05 = positions(0.5)
    ph025, ah025 = positions(0.25)
    phlin, ahlin = positions(0.5, linear=True)
    coarse = _a2_position_metrics(ph05, ah05)
    refined = _a2_position_metrics(ph025, ah025)
    method = _a2_position_metrics(ph05 - ah05, phlin - ahlin)
    resampling = {
        "rms_m": max(abs(refined["rms_m"] - coarse["rms_m"]),
                     method["rms_m"]),
        "max_m": max(abs(refined["max_m"] - coarse["max_m"]),
                     method["max_m"]),
        "landing_m": 0.0, "group_rel": 0.0}
    if fixed:
        ph_end, ph_group = ph_endpoint(fine)
        ph_linear = _a2_sample(fine["ph_s"], fine["ph_xyz"],
                               np.array([100.0]), linear=True)[0]
        linear_group = float(_a2_sample(
            fine["ph_s"], fine["ph_group"], np.array([100.0]),
            linear=True)[0])
        resampling["landing_m"] = float(np.linalg.norm(
            ph_end - ph_linear) * 1000.0)
        resampling["group_rel"] = abs(ph_group - linear_group) / ph_group

    actual = {
        "rms_m": float(rows[2]["position_rms_m"]),
        "max_m": float(rows[2]["position_max_m"]),
        "landing_m": float(rows[2]["position_endpoint_m"]),
        "group_rel": float(rows[2]["group_relative_error"])}
    checks = {}
    for key in limits:
        last = d23["ph"][key] + d23["ah"][key]
        previous = d12["ph"][key] + d12["ah"][key]
        budget = max(floor[key], 3.0 * (last + resampling[key]))
        checks[key] = {
            "residual": actual[key], "physical_limit": limits[key],
            "delta12_ph": d12["ph"][key], "delta12_ah": d12["ah"][key],
            "delta23_ph": d23["ph"][key], "delta23_ah": d23["ah"][key],
            "resampling": resampling[key], "three_sigma_budget": budget,
            "physical_pass": actual[key] <= limits[key],
            "explained_by_numeric": actual[key] <= budget,
            "self_convergence_pass": last <= 0.2 * limits[key],
            "monotonic_or_roundoff": (last <= previous + floor[key] or
                                      last <= floor[key])}
    passed = all(all(values[field] for field in (
        "physical_pass", "explained_by_numeric", "self_convergence_pass",
        "monotonic_or_roundoff")) for values in checks.values())
    return {"status": "pass" if passed else "needs_extra_levels",
            "checks": checks}


def _verify_a2_formal_summaries(root: Path) -> list[str]:
    """Recalculate formal comparisons and gates from immutable PH/AH rays."""
    failures: list[str] = []
    required = {
        "base": ("identified_hermite_all", "common_hermite_smooth",
                 "common_linear_all",
                 "identified_hermite_smooth_fd1p0",
                 "identified_hermite_smooth_fd0p25",
                 "identified_hermite_smooth_fd0p125"),
        "identified": ("identified_hermite_internal",),
        "identified_holdout": ("identified_hermite_affine",)}
    for dataset, suffixes in required.items():
        directory = root if dataset == "base" else root / dataset
        analysis_dir = directory / "analysis"
        for suffix in suffixes:
            path = analysis_dir / f"a2_{dataset}_formal_{suffix}.json"
            if not path.is_file():
                failures.append(f"A2 formal analysis missing: {path}")
        for path in sorted(analysis_dir.glob(
                f"a2_{dataset}_formal_*.json")):
            try:
                report = json.loads(path.read_text(encoding="utf-8"))
                manifest = json.loads((directory / "input_manifest.json").read_text(
                    encoding="utf-8"))
                subset = str(report["subset"])
                jobs = [job for job in manifest["jobs"] if
                        subset == "all" or job["stage"] == subset]
                jobs_by_id = {str(job["job_id"]): job for job in jobs}
                rows = report["rows"]
                rows_by_id = {str(row["job_id"]): row for row in rows}
                if (report.get("schema") != "a4-c1m-a2-analysis-v1" or
                        report.get("dataset") != dataset or
                        report.get("tag") != "formal" or
                        report.get("input_manifest_sha256") !=
                        experiment.sha(directory / "input_manifest.json") or
                        not (len(rows) == len(rows_by_id) == len(jobs)) or
                        set(rows_by_id) != set(jobs_by_id)):
                    failures.append(f"A2 analysis identity/coverage mismatch: {path}")
                    continue
                coefficient = float(report["coefficient_K2"])
                if (report["interpolation"] == "hermite" and
                        abs(coefficient - 80.6404) > 1e-10):
                    freeze_path = root / "identified" / "coefficient_freeze.json"
                    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
                    if (abs(coefficient - float(freeze["identified_K2"])) >
                            1e-12 or
                            report.get("coefficient_freeze_sha256") !=
                            experiment.sha(freeze_path)):
                        failures.append(f"A2 analysis frozen K2 mismatch: {path}")
                label = ("common" if abs(coefficient - 80.6404) < 1e-10
                         else "identified")
                fd = report.get("spatial_fd_base_km")
                base_name = (f"a2_{dataset}_formal_{label}_"
                             f"{report['interpolation']}_{subset}")
                tagged_name = (
                    f"{base_name}_fd{str(float(fd)).replace('.', 'p')}.json"
                    if fd is not None else "")
                tagged_fd = path.name == tagged_name
                if (path.name != f"{base_name}.json" and not tagged_fd):
                    failures.append(f"A2 formal analysis path/parameter mismatch: {path}")
                ah_dir = (directory / "ah" / "formal" / label /
                          str(report["interpolation"]))
                if tagged_fd:
                    ah_dir = ah_dir / f"fd{float(fd):g}"
                traces: dict[str, dict] = {}
                group_rows: dict[str, list[dict]] = {}
                group_jobs: dict[str, list[dict]] = {}
                for row in rows:
                    name = str(row["job_id"])
                    job = jobs_by_id[name]
                    key = (f"{job['case_id']}|f{float(job['frequency_MHz']):g}"
                           f"|e{float(job['elevation_deg']):.14g}")
                    group_rows.setdefault(key, []).append(row)
                    group_jobs.setdefault(key, []).append(job)
                    raw_path = directory / "pharlap_raw" / f"{name}.mat"
                    ah_path = ah_dir / f"{name}.json"
                    status = ("missing" if not ah_path.is_file() else
                              json.loads(ah_path.read_text(
                                  encoding="utf-8"))["status"])
                    if (row.get("ph_raw_sha256") != experiment.sha(raw_path) or
                            row.get("ah_status") != status or
                            row.get("stage") != job["stage"] or
                            int(row.get("level", -1)) != int(job["level"])):
                        failures.append(f"A2 analysis raw row mismatch: {name}")
                    if status == "missing":
                        continue
                    trace = _a2_native_trace(directory, ah_dir, job)
                    traces[name] = trace
                    ah = trace["ah_record"]
                    p2 = float(trace["ph_p2"][0])
                    ah_n2 = ah.get("initial_n2_prescribed")
                    same_initial = (ah_n2 is not None and
                                    abs(p2 - float(ah_n2)) <= 1e-10)
                    if (row.get("same_initial_state") is not same_initial or
                            not _a2_close(row.get("ph_initial_p2"), p2,
                                          atol=1e-12)):
                        failures.append(f"A2 analysis initial-state gate drift: {name}")
                    if status != "completed":
                        if any(field in row for field in (
                                "position_rms_m", "position_endpoint_m",
                                "group_relative_error")):
                            failures.append(
                                f"A2 rejected AH ray has comparison metrics: {name}")
                        continue
                    fixed = job["stage"] in ("internal", "affine")
                    end = (100.0 if fixed else min(
                        float(trace["ph_s"][-1]), float(trace["ah_s"][-1])))
                    arc = np.linspace(0.0, end, 201 if fixed else 501)
                    ph_xyz = _a2_sample(trace["ph_s"], trace["ph_xyz"], arc)
                    ah_xyz = _a2_sample(trace["ah_s"], trace["ah_xyz"], arc)
                    metric = _a2_position_metrics(ph_xyz, ah_xyz)
                    ph_endpoint = (ph_xyz[-1] if fixed else
                                   trace["ph_xyz"][-1])
                    landing = float(np.linalg.norm(
                        ph_endpoint - trace["ah_xyz"][-1]) * 1000.0)
                    ph_group = (float(_a2_sample(
                        trace["ph_s"], trace["ph_group"],
                        np.array([100.0]))[0]) if fixed else
                        float(trace["ph_group"][-1]))
                    group_rel = abs(float(ah["group_path_km_native"]) -
                                    ph_group) / ph_group
                    numeric_fields = {
                        "comparison_arc_km": end,
                        "position_rms_m": metric["rms_m"],
                        "position_max_m": metric["max_m"],
                        "position_endpoint_m": landing,
                        "ph_group_for_gate_km": ph_group,
                        "group_relative_error": group_rel}
                    if (row.get("comparison_fixed_100km") is not fixed or
                            row.get("group_comparison_kind") !=
                            ("fixed_arc_100km" if fixed else
                             "natural_terminal") or
                            any(not _a2_close(row.get(field), value)
                                for field, value in numeric_fields.items())):
                        failures.append(f"A2 analysis numeric row drift: {name}")
                gates = report["gates"]
                gates_by_case = {str(item["case"]): item for item in gates}
                if (not (len(gates) == len(gates_by_case) ==
                         len(group_rows)) or
                        set(gates_by_case) != set(group_rows)):
                    failures.append(f"A2 formal gate coverage mismatch: {path}")
                    continue
                for case, level_rows in group_rows.items():
                    gate = gates_by_case[case]
                    level_rows.sort(key=lambda row: int(row["level"]))
                    level_jobs = sorted(
                        group_jobs[case], key=lambda job: int(job["level"]))
                    if ([int(row["level"]) for row in level_rows] != [1, 2, 3]
                            or gate.get("levels") != [1, 2, 3]):
                        failures.append(f"A2 formal gate levels invalid: {case}")
                        continue
                    if any(row["ah_status"] != "completed"
                           for row in level_rows):
                        expected_status = "native_incomplete"
                        numeric = None
                    elif any(row.get("same_initial_state") is not True
                             for row in level_rows):
                        expected_status = "initial_state_mismatch"
                        numeric = None
                    else:
                        level_traces = [traces[str(job["job_id"])]
                                        for job in level_jobs]
                        numeric = _a2_gate_from_traces(
                            level_jobs, level_traces, level_rows)
                        expected_status = numeric["status"]
                    if (gate.get("gate") != expected_status or
                            gate.get("numeric", {}).get("status") !=
                            expected_status):
                        failures.append(f"A2 formal gate status drift: {case}")
                    if numeric is None:
                        continue
                    actual_checks = gate["numeric"].get("checks", {})
                    if set(actual_checks) != set(numeric["checks"]):
                        failures.append(f"A2 formal gate metrics missing: {case}")
                        continue
                    for field, expected in numeric["checks"].items():
                        actual = actual_checks[field]
                        for key, value in expected.items():
                            if isinstance(value, bool):
                                valid = type(actual.get(key)) is bool and (
                                    actual[key] is value)
                            else:
                                valid = _a2_close(actual.get(key), value)
                            if not valid:
                                failures.append(
                                    f"A2 formal gate {field}/{key} drift: {case}")
                if report.get("completed_ah_count") != sum(
                        row["ah_status"] == "completed" for row in rows):
                    failures.append(f"A2 formal AH completed count drift: {path}")
                if report.get("missing_ah_count") != sum(
                        row["ah_status"] == "missing" for row in rows):
                    failures.append(f"A2 formal AH missing count drift: {path}")
            except (KeyError, ValueError, OSError, TypeError, IndexError) as exc:
                failures.append(f"A2 formal analysis unreadable: {path}: {exc}")
    return failures


def _verify_a2_spatial_fd(root: Path) -> tuple[list[str], list[dict]]:
    """Compare the three spatial-difference scales at fixed finest ODE tolerance."""
    failures: list[str] = []
    results: list[dict] = []
    jobs = [job for job in json.loads(
        (root / "input_manifest.json").read_text(encoding="utf-8"))["jobs"]
        if job["stage"] == "smooth" and int(job["level"]) == 3]
    if len(jobs) != 4:
        return ["A2 spatial FD smooth four-angle inventory mismatch"], results
    limit = {"rms_m": 10.0, "max_m": 50.0, "landing_m": 50.0,
             "group_rel": 1e-5}
    floor = {"rms_m": 1e-7, "max_m": 1e-7, "landing_m": 1e-7,
             "group_rel": 1e-13}
    for job in jobs:
        name = str(job["job_id"])
        try:
            traces = []
            for fd in ("1", "0.25", "0.125"):
                directory = (root / "ah" / "formal" / "identified" /
                             "hermite" / f"fd{fd}")
                record = json.loads((directory / f"{name}.json").read_text(
                    encoding="utf-8"))
                with np.load(directory / f"{name}.npz") as data:
                    s = np.asarray(data["geom_path_km"], float)
                    xyz = np.asarray(data["xyz_km"], float)
                if (record["status"] != "completed" or
                        record["termination_status"] != "ground" or
                        abs(float(record["spatial_fd_base_km"]) - float(fd)) >
                        1e-12):
                    failures.append(f"A2 spatial FD native ray invalid: {name}/{fd}")
                traces.append((s, xyz,
                               float(record["group_path_km_native"])))
            raw = loadmat(root / "pharlap_raw" / f"{name}.mat",
                          simplify_cells=True)
            ph_group = float(np.asarray(
                raw["ray_state_vec"]["group_path"], float).reshape(-1)[-1] /
                1000.0)
            deltas = []
            for first, second in ((traces[0], traces[1]),
                                  (traces[1], traces[2])):
                s1, xyz1, group1 = first
                s2, xyz2, group2 = second
                arc = np.linspace(0.0, min(float(s1[-1]), float(s2[-1])),
                                  501)
                metric = _a2_position_metrics(
                    _a2_sample(s1, xyz1, arc),
                    _a2_sample(s2, xyz2, arc))
                metric["landing_m"] = float(
                    np.linalg.norm(xyz1[-1] - xyz2[-1]) * 1000.0)
                metric["group_rel"] = abs(group1 - group2) / ph_group
                deltas.append(metric)
            checks = {
                key: {"monotonic_or_roundoff":
                      deltas[1][key] <= deltas[0][key] + floor[key],
                      "self_convergence_pass":
                      deltas[1][key] <= 0.2 * limit[key]}
                for key in limit}
            if any(not all(flags.values()) for flags in checks.values()):
                failures.append(f"A2 spatial FD self-convergence failed: {name}")
            results.append({"job_id": name, "fd1_to_fd0p25": deltas[0],
                            "fd0p25_to_fd0p125": deltas[1],
                            "checks": checks})
        except (KeyError, ValueError, OSError, TypeError, IndexError) as exc:
            failures.append(f"A2 spatial FD unreadable: {name}: {exc}")
    return failures, results


def _verify_a2_cause_matrix(root: Path) -> list[str]:
    failures: list[str] = []
    analysis = root / "analysis"
    matrix_path = analysis / "a2_cause_matrix.json"
    if not matrix_path.is_file():
        return ["A2 cause matrix missing"]
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    if matrix.get("schema") != "a4-c1m-a2-cause-matrix-v1":
        failures.append("A2 cause matrix schema mismatch")
    snapshot = Path(matrix.get("analysis_runner_snapshot", ""))
    if (not snapshot.is_file() or
            snapshot.parent.resolve() !=
            (root / "sources" / "analysis").resolve() or
            experiment.sha(snapshot) != matrix.get("analysis_runner_sha256")):
        failures.append("A2 cause matrix analysis-source snapshot drift")
    replay_path = analysis / "a2_prepare_replay.json"
    if (not replay_path.is_file() or
            experiment.sha(replay_path) !=
            matrix.get("prepare_replay_sha256")):
        failures.append("A2 prepare replay SHA drift")
    else:
        replay = json.loads(replay_path.read_text(encoding="utf-8"))
        original = json.loads((root / "input_manifest.json").read_text(
            encoding="utf-8"))
        if (replay.get("schema") != "a4-c1m-a2-input-replay-v1" or
                replay.get("original_input_manifest_sha256") !=
                experiment.sha(root / "input_manifest.json") or
                replay.get("original_runner_sha256_at_prepare") !=
                original["runner_sha256_at_prepare"] or
                replay.get("physical_array_and_job_parameters_match") is
                not True or replay.get("grid_count") != 22 or
                replay.get("job_count") != 210 or replay.get("mismatches")):
            failures.append("A2 prepare replay evidence mismatch")
    freeze_path = root / "identified" / "coefficient_freeze.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if (abs(float(matrix["identified_K2"]) -
            float(freeze["identified_K2"])) > 1e-12 or
            matrix.get("coefficient_freeze_sha256") !=
            experiment.sha(freeze_path)):
        failures.append("A2 cause matrix coefficient-freeze drift")
    summary_paths = {
        "common_linear": analysis / "a2_base_formal_common_linear_all.json",
        "common_hermite": analysis / "a2_base_formal_common_hermite_smooth.json",
        "internal_common_hermite": analysis /
        "a2_base_formal_common_hermite_internal.json",
        "identified_linear": analysis /
        "a2_base_formal_identified_linear_smooth.json",
        "identified_hermite": analysis /
        "a2_base_formal_identified_hermite_all.json",
        "identified_hermite_fd025": analysis /
        "a2_base_formal_identified_hermite_smooth_fd0p25.json",
        "identified_hermite_fd0125": analysis /
        "a2_base_formal_identified_hermite_smooth_fd0p125.json",
        "identified_hermite_fd1": analysis /
        "a2_base_formal_identified_hermite_smooth_fd1p0.json",
        "internal_identified": root / "identified" / "analysis" /
        "a2_identified_formal_identified_hermite_internal.json",
        "holdout_identified": root / "identified_holdout" / "analysis" /
        "a2_identified_holdout_formal_identified_hermite_affine.json"}
    if set(matrix.get("summary_sha256", {})) != set(summary_paths):
        failures.append("A2 cause matrix summary inventory mismatch")
    summaries = {}
    for name, path in summary_paths.items():
        if (not path.is_file() or
                matrix["summary_sha256"].get(name) != experiment.sha(path)):
            failures.append(f"A2 cause matrix summary SHA drift: {name}")
            continue
        summaries[name] = json.loads(path.read_text(encoding="utf-8"))
    artifacts = {
        "per_case_csv_sha256": analysis / "a2_per_case.csv",
        "per_point_csv_sha256": analysis / "a2_per_point.csv",
        "difference_png_sha256": analysis / "a2_smooth_difference.png"}
    if set(matrix.get("artifacts", {})) != set(artifacts):
        failures.append("A2 cause matrix artifact inventory mismatch")
    for key, path in artifacts.items():
        if not path.is_file() or matrix["artifacts"].get(key) != experiment.sha(path):
            failures.append(f"A2 cause matrix artifact SHA drift: {key}")
    case_csv = artifacts["per_case_csv_sha256"]
    point_csv = artifacts["per_point_csv_sha256"]
    if not case_csv.is_file() or not point_csv.is_file():
        return failures
    with case_csv.open(newline="", encoding="utf-8") as stream:
        case_rows = list(csv.DictReader(stream))
    if len(case_rows) != 540:
        failures.append(f"A2 per-case CSV count: {len(case_rows)}/540")
    by_candidate = {}
    for name, summary in summaries.items():
        by_candidate[name] = {
            row["job_id"]: row for row in summary["rows"]}
    for row in case_rows:
        name, job_id = row["candidate"], row["job_id"]
        source = by_candidate.get(name, {}).get(job_id)
        if source is None or row["ah_status"] != source["ah_status"]:
            failures.append(f"A2 per-case CSV row identity drift: {name}/{job_id}")
            continue
        for csv_key, source_key in (
                ("position_rms_m", "position_rms_m"),
                ("position_max_m", "position_max_m"),
                ("position_endpoint_m", "position_endpoint_m"),
                ("group_relative_error", "group_relative_error")):
            if source.get(source_key) is None:
                if row[csv_key]:
                    failures.append(
                        f"A2 per-case CSV rejected-ray metric: {name}/{job_id}")
            elif not _a2_close(float(row[csv_key]), float(source[source_key])):
                failures.append(
                    f"A2 per-case CSV metric drift: {name}/{job_id}/{csv_key}")
    with point_csv.open(newline="", encoding="utf-8") as stream:
        point_count = sum(1 for _ in csv.DictReader(stream))
    if point_count != 122460:
        failures.append(f"A2 per-point CSV count: {point_count}/122460")
    if "identified_hermite_fd0125" in summaries:
        fine = summaries["identified_hermite_fd0125"]
        if (len(fine["gates"]) != 4 or
                any(gate["gate"] != "pass" for gate in fine["gates"])):
            failures.append("A2 finest FD smooth four-angle gates not all pass")
    smooth = matrix.get("smooth_level3", {})
    for name, summary in summaries.items():
        for row in summary["rows"]:
            if row["stage"] != "smooth" or int(row["level"]) != 3:
                continue
            angle = f"{float(row['elevation_deg']):.14g}"
            values = smooth.get(angle, {}).get(name)
            if values is None:
                failures.append(f"A2 cause matrix smooth cell missing: {angle}/{name}")
                continue
            for key, source in (
                    ("rms_m", "position_rms_m"),
                    ("max_m", "position_max_m"),
                    ("natural_landing_m", "position_endpoint_m"),
                    ("natural_group_rel", "group_relative_error")):
                if not _a2_close(values.get(key), float(row[source])):
                    failures.append(
                        f"A2 cause matrix smooth value drift: {angle}/{name}/{key}")
    if "identified_hermite" in summaries:
        rows = summaries["identified_hermite"]["rows"]
        for hmin in (20, 40, 60, 80):
            selected = [row for row in rows
                        if row["case_id"] == f"boundary_h{hmin:03d}"]
            cell = matrix.get("boundary", {}).get(str(hmin), {})
            expected_first = [
                float(row["ph_Ne_first_interval_rms_m3"])
                for row in selected]
            if (cell.get("job_count") != len(selected) or
                    cell.get("native_completed") != sum(
                        row["ah_status"] == "completed" for row in selected) or
                    cell.get("native_rejected") != sum(
                        row["ah_status"] == "physics_domain_rejected"
                        for row in selected) or
                    not _a2_close(cell.get("max_PH_Ne_below_grid_m3"),
                                  max(float(row["ph_Ne_below_grid_max_m3"])
                                      for row in selected)) or
                    len(cell.get("first_interval_PH_vs_Hermite_rms_m3", []))
                    != len(expected_first) or
                    any(not _a2_close(observed, expected, atol=1e-6)
                        for observed, expected in zip(
                            cell["first_interval_PH_vs_Hermite_rms_m3"],
                            expected_first))):
                failures.append(f"A2 cause matrix boundary drift: hmin={hmin}")
    if "holdout_identified" in summaries:
        held = summaries["holdout_identified"]
        if (len(held["gates"]) != 16 or
                any(gate["gate"] != "pass" for gate in held["gates"])):
            failures.append("A2 matched 8/15MHz holdout gates not all pass")
        held_matrix = matrix.get("matched_holdout", {})
        finest_rows = [row for row in held["rows"] if int(row["level"]) == 3]
        if (held_matrix.get("jobs") != len(held["rows"]) or
                held_matrix.get("gate_counts") != {"pass": 16} or
                not _a2_close(held_matrix.get("max_level3_rms_m"),
                              max(float(row["position_rms_m"])
                                  for row in finest_rows)) or
                not _a2_close(held_matrix.get("max_level3_group_rel"),
                              max(float(row["group_relative_error"])
                                  for row in finest_rows))):
            failures.append("A2 cause matrix holdout drift")
    return failures


def _verify_a2_production_identity(root: Path) -> list[str]:
    path = root / "production_and_mex_identity.json"
    if not path.is_file():
        return ["A2 production/MEX identity record missing"]
    record = json.loads(path.read_text(encoding="utf-8"))
    failures = []
    if record.get("all_unchanged_since_a0") is not True:
        failures.append("A2 production/MEX unchanged aggregate is not true")
    for name, item in record.get("files", {}).items():
        source = Path(item["path"])
        if (not source.is_file() or
                experiment.sha(source) != item.get("sha256") or
                item.get("sha256") != item.get("a0_sha256") or
                item.get("unchanged_since_a0") is not True):
            failures.append(f"A2 production/MEX hash changed: {name}")
    if not record.get("files"):
        failures.append("A2 production/MEX identity inventory empty")
    return failures


def _verify_a2_ah_files(root: Path, base_jobs: list[dict],
                        identified_jobs: list[dict],
                        holdout_jobs: list[dict]) -> list[str]:
    failures = []
    verified_snapshots: set[str] = set()
    for dataset, directory, jobs in (
        ("base", root / "ah" / "formal", base_jobs),
        ("identified", root / "identified" / "ah" / "formal", identified_jobs),
        ("identified_holdout", root / "identified_holdout" / "ah" / "formal",
         holdout_jobs),
    ):
        by_id = {job["job_id"]: job for job in jobs}
        if not directory.is_dir():
            failures.append(f"A2 {dataset} formal AH directory missing")
            continue
        required = directory / "identified" / "hermite"
        actual_ids = [path.stem for path in required.glob("*.json")]
        if (len(actual_ids) != len(set(actual_ids)) or
                set(actual_ids) != set(by_id) or
                len(actual_ids) != len(jobs)):
            failures.append(
                f"A2 {dataset} formal identified/Hermite AH coverage incomplete: "
                f"{len(actual_ids)}/{len(jobs)}")
        for path in sorted(directory.rglob("*.json")):
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
                job = by_id[record["job_id"]]
                if record.get("schema") != "a4-c1m-a2-ah-raw-v1":
                    failures.append(f"A2 AH schema mismatch: {path}")
                if (record.get("job_sha256") != job["job_sha256"] or
                    record.get("grid_sha256") != job["grid_sha256"]):
                    failures.append(f"A2 AH job/grid identity mismatch: {path}")
                parts = path.relative_to(directory).parts
                standard = (
                    len(parts) == 3 and parts[0] in ("common", "identified")
                    and parts[1] in ("linear", "hermite"))
                fd_variant = (
                    len(parts) == 4 and parts[0] in ("common", "identified")
                    and parts[1] in ("linear", "hermite")
                    and parts[2].startswith("fd")
                    and parts[2][2:].replace(".", "", 1).isdigit()
                    and float(parts[2][2:]) ==
                    float(record.get("spatial_fd_base_km", float("nan"))))
                if path.stem != job["job_id"] or not (standard or fd_variant):
                    failures.append(f"A2 AH artifact path mismatch: {path}")
                source = record.get("source_sha256", {})
                runner_sha = source.get("runner", "")
                snapshot = root / "sources" / "ah" / runner_sha[:16] / "source_manifest.json"
                if not runner_sha or not snapshot.is_file():
                    failures.append(f"A2 AH source snapshot missing: {path}")
                elif runner_sha not in verified_snapshots:
                    snapshot_data = json.loads(snapshot.read_text(encoding="utf-8"))
                    files = snapshot_data.get("files", {})
                    if (snapshot_data.get("schema") != "a4-c1m-a2-ah-source-snapshot-v1"
                        or set(files) != set(source)):
                        failures.append(f"A2 AH source snapshot schema mismatch: {snapshot}")
                    for key, value in files.items():
                        frozen_path = Path(value["snapshot"])
                        if (not frozen_path.is_file() or
                            experiment.sha(frozen_path) != value["sha256"] or
                            source.get(key) != value["sha256"]):
                            failures.append(f"A2 AH frozen source mismatch: {snapshot}/{key}")
                    verified_snapshots.add(runner_sha)
                status = record.get("status")
                if status not in {"completed", "physics_domain_rejected",
                                  "incomplete_max_path", "failed"}:
                    failures.append(f"A2 AH invalid native status: {path}")
                if status == "failed":
                    if not record.get("error"):
                        failures.append(f"A2 AH failure lacks reason: {path}")
                    continue
                if (status == "physics_domain_rejected" and
                    (record.get("termination_status") != "physics_domain_rejected" or
                     not record.get("failure_reason"))):
                    failures.append(f"A2 AH physics rejection not preserved: {path}")
                raw_path = path.with_suffix(".npz")
                if (not raw_path.is_file() or
                    experiment.sha(raw_path) != record.get("raw_npz_sha256")):
                    failures.append(f"A2 AH NPZ missing/hash mismatch: {path}")
                    continue
                with np.load(raw_path, allow_pickle=False) as raw:
                    keys = {"xyz_km", "geom_path_km", "wavevector_km_inv",
                            "group_tangent_ecef", "lat_lon_height", "Ne_m3_native",
                            "n2_native"}
                    if not keys.issubset(raw.files):
                        failures.append(f"A2 AH native state arrays missing: {raw_path}")
                        continue
                    n = len(raw["geom_path_km"])
                    s = np.asarray(raw["geom_path_km"], float)
                    if (n != record["point_count"] or raw["xyz_km"].shape != (n, 3)
                        or raw["wavevector_km_inv"].shape != (n, 3)
                        or not np.all(np.isfinite(s)) or
                        not np.all(np.diff(s) > 0)):
                        failures.append(f"A2 AH native state shape/arc invalid: {raw_path}")
                    if (status == "completed" and job["stage"] in ("affine", "internal")
                        and (record.get("termination_status") !=
                             "fixed_geometric_path_completed" or abs(s[-1] - 100) > 1e-8)):
                        failures.append(f"A2 AH fixed 100km path did not complete: {path}")
                    if (status == "completed" and job["stage"] in ("boundary", "smooth")
                        and record.get("termination_status") != "ground"):
                        failures.append(f"A2 AH full-ray did not land: {path}")
            except (KeyError, ValueError, OSError, TypeError) as exc:
                failures.append(f"A2 AH artifact unreadable: {path}: {exc}")
    return failures


def verify_a2(run: Path) -> list[str]:
    """Read frozen A2 inputs and PH raw without starting either solver."""
    failures: list[str] = []

    def check(ok: bool, message: str) -> None:
        if not ok:
            failures.append(message)

    root = run / "a2"
    manifest = json.loads((root / "input_manifest.json").read_text(encoding="utf-8"))
    check(manifest.get("schema") == "a4-c1m-a2-input-v1", "A2 input schema mismatch")
    check(experiment.sha(root / "input.mat") == manifest["input_mat_sha256"],
          "A2 input MAT SHA mismatch")
    check(experiment.sha(Path(manifest["source_c1_mat"])) ==
          manifest["source_c1_sha256"], "frozen C1 Chapman input drift")
    prepared_runner = root / "sources" / "runner_at_prepare.py"
    runner_source = (prepared_runner if prepared_runner.is_file()
                     else ROOT / "tests" / "_run_a4_c1m_attribution.py")
    check(experiment.sha(runner_source) == manifest["runner_sha256_at_prepare"],
          "A2 prepare runner SHA cannot be verified against frozen source")
    matlab_snapshot = root / "sources" / "a4_c1m_uniform_export.m"
    check(matlab_snapshot.is_file(), "A2 MATLAB source snapshot missing")
    jobs = _a2_jobs(root / "input.mat")
    identities = [str(job["job_id"]) for job in jobs]
    check(len(jobs) == 210 and len(set(identities)) == 210,
          "expected 210 distinct A2 jobs")
    stages = {stage: [job for job in jobs if job["stage"] == stage]
              for stage in ("boundary", "smooth", "internal", "affine")}
    check({key: len(value) for key, value in stages.items()} ==
          {"boundary": 48, "smooth": 12, "internal": 6, "affine": 144},
          "A2 stage job counts mismatch")
    train = [job for job in stages["affine"] if float(job["frequency_MHz"]) == 10.0]
    holdout = [job for job in stages["affine"] if float(job["frequency_MHz"]) in (8.0, 15.0)]
    check(len(train) == 96 and len(holdout) == 48 and
          set(identities) == {job["job_id"] for job in
                             stages["boundary"] + stages["smooth"] +
                             stages["internal"] + train + holdout},
          "A2 training/holdout split mismatch")
    check(all(int(job["explicit_state"]) == (job["stage"] in ("internal", "affine"))
              for job in jobs), "A2 launch mode mismatch")
    check(all(np.array_equal(np.asarray(job["tol"], float).reshape(-1),
                             experiment.TOL[int(job["level"]) - 1]) for job in jobs),
          "A2 PH tolerance differs from public legal levels")
    grid_by_hash: dict[str, dict] = {}
    for case, record in manifest["grid"].items():
        grid_path = Path(record["grid_mat"])
        digest = record["grid_sha256"]
        if not grid_path.is_file() or experiment.sha(grid_path) != digest:
            failures.append(f"A2 grid missing/drifted: {case}")
            continue
        grid = loadmat(grid_path, simplify_cells=True)
        readback_path = root / "pharlap_raw" / f"grid_readback_{digest[:16]}.mat"
        if not readback_path.is_file():
            failures.append(f"A2 MATLAB grid readback missing: {case}")
            continue
        readback = loadmat(readback_path, simplify_cells=True)
        for key, original in grid.items():
            if key.startswith("__"):
                continue
            check(key in readback and np.array_equal(original, readback[key]),
                  f"A2 MATLAB grid array mismatch: {case}/{key}")
        h = np.asarray(grid["height_km"], float).reshape(-1)
        bh = np.asarray(grid["geomag_height_km"], float).reshape(-1)
        ne = np.asarray(grid["Ne_m3"], float)
        profile = np.asarray(grid["ne_profile_m3"], float).reshape(-1)
        check(2 <= len(h) <= 401 and len(bh) <= 201 and
              np.all(np.diff(h) > 0), f"A2 grid dimensions invalid: {case}")
        check(ne.shape == (75, 9, len(h)) and
              np.array_equal(ne[0, 0], profile) and
              np.array_equal(ne[-1, -1], profile),
              f"A2 grid axis/order/radial values invalid: {case}")
        check(all(np.all(np.asarray(grid[key]) == 0) for key in
                  ("Bx_T", "By_T", "Bz_T", "collision_freq_Hz")),
              f"A2 nonzero B or collisions: {case}")
        if case.startswith("affine_") and case.endswith("d0p25"):
            check(len(h) == 401 and h[0] == 80 and h[-1] == 180,
                  f"A2 finest affine grid is not canonical 401 nodes: {case}")
        grid_by_hash[digest] = grid
    summary_path = root / "pharlap_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    check(summary.get("schema_version") == "c1m-a2-export-v1" and
          summary.get("status") == "completed_a2_export" and
          summary.get("failures") == 0, "A2 PH export incomplete")
    records = summary["records"]
    check(len(records) == len(jobs) and
          {record["job_id"] for record in records} == set(identities),
          "A2 PH raw identities mismatch")
    record_by_id = {record["job_id"]: record for record in records}
    check(all(type(record["reused"]) is bool for record in records),
          "A2 PH reused flags are not JSON booleans")
    state_fields = {"pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z",
                    "group_path", "geometrical_path", "phase_path",
                    "absorption", "indep_var", "ODE_step_size"}
    path_fields = {"lat", "lon", "height", "electron_density",
                   "refractive_index", "group_refractive_index",
                   "geometric_distance", "group_range", "geomag_x",
                   "geomag_y", "geomag_z", "collision_frequency"}
    for job in jobs:
        name = job["job_id"]
        record = record_by_id.get(name)
        if record is None:
            continue
        raw_path = root / "pharlap_raw" / f"{name}.mat"
        if (not raw_path.is_file() or Path(record["raw_mat"]).resolve() != raw_path.resolve()
            or record["job_sha256"] != job["job_sha256"]
            or record["grid_sha256"] != job["grid_sha256"]):
            failures.append(f"A2 PH raw identity/hash mismatch: {name}")
            continue
        raw = loadmat(raw_path, simplify_cells=True)
        check(raw["job_sha256"] == job["job_sha256"] and
              raw["grid_sha256"] == job["grid_sha256"],
              f"A2 PH embedded hashes mismatch: {name}")
        state, path = raw["ray_state_vec"], raw["ray_path_data"]
        if not state_fields.issubset(state) or not path_fields.issubset(path):
            failures.append(f"A2 PH missing required state/path fields: {name}")
            continue
        length = len(np.asarray(path["height"]).reshape(-1))
        check(length == int(record["point_count"]) and all(
            len(np.asarray(state[key]).reshape(-1)) == length for key in state_fields),
            f"A2 PH state/path length mismatch: {name}")
        check(np.max(np.abs(np.asarray(path["group_range"], float).reshape(-1)
                            - np.asarray(state["group_path"], float).reshape(-1) / 1000))
              < 1e-7, f"A2 PH group-path unit mismatch: {name}")
        check(all(np.all(np.asarray(path[key], float) == 0) for key in
                  ("geomag_x", "geomag_y", "geomag_z", "collision_frequency")),
              f"A2 PH reported nonzero B/collision: {name}")
        if int(job["explicit_state"]):
            initial = raw["initial"]
            check(isinstance(initial, dict) and
                  {"pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z",
                   "group_path", "geom_path", "phase_path", "absorption",
                   "indep_var", "ODE_step_size"}.issubset(initial),
                  f"A2 explicit initial state incomplete: {name}")
            if isinstance(initial, dict):
                grid = grid_by_hash.get(job["grid_sha256"])
                if grid is not None:
                    h = np.asarray(grid["height_km"], float).reshape(-1)
                    profile = np.asarray(grid["ne_profile_m3"], float).reshape(-1)
                    ne0 = float(np.interp(float(job["origin_height_km"]), h, profile))
                    n2 = 1 - float(job["launch_k2"]) * ne0 / (
                        float(job["frequency_MHz"]) * 1e6) ** 2
                    direction = experiment.launch_direction_ecef(
                        30.0, 120.0, float(job["elevation_deg"]), 0.0, experiment.R_KM)
                    expected_pos = experiment.geocentric_to_ecef(
                        30.0, 120.0, float(job["origin_height_km"]),
                        experiment.R_KM) * 1000
                    pos = np.array([float(initial[key]) for key in
                                    ("pos_x", "pos_y", "pos_z")])
                    wave = np.array([float(initial[key]) for key in
                                     ("dir_x", "dir_y", "dir_z")])
                    first_pos = np.array([np.asarray(state[key], float).reshape(-1)[0]
                                          for key in ("pos_x", "pos_y", "pos_z")])
                    first_wave = np.array([np.asarray(state[key], float).reshape(-1)[0]
                                           for key in ("dir_x", "dir_y", "dir_z")])
                    check(n2 > 0 and np.linalg.norm(pos - expected_pos) < 1e-3 and
                          np.linalg.norm(wave - np.sqrt(n2) * direction) < 1e-10 and
                          np.linalg.norm(first_pos - pos) < 1e-3 and
                          np.linalg.norm(first_wave - wave) < 1e-10,
                          f"A2 explicit launch position/wave normal mismatch: {name}")
        if job["stage"] == "boundary":
            grid = grid_by_hash.get(job["grid_sha256"])
            if grid is not None:
                below = np.asarray(path["height"], float) < float(np.min(grid["height_km"]))
                if np.any(below):
                    ne_below = np.asarray(path["electron_density"], float)[below]
                    check(np.all(ne_below == 0), f"A2 PH below-hmin medium not vacuum: {name}")
    identified = root / "identified"
    id_jobs: list[dict] = []
    if identified.is_dir():
        freeze_path = identified / "coefficient_freeze.json"
        freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
        check(freeze.get("schema") == "a4-c1m-a2-identified-k2-v1" and
              float(freeze["delta_Ne_min_m3"]) == 1e8,
              "A2 identified coefficient freeze/schema mismatch")
        selected = []
        for evidence in freeze["training"]:
            raw_path = Path(evidence["raw_mat"])
            name = raw_path.stem
            check(name.startswith("affine_") and "_d0p25_f10_" in name and
                  name.endswith("_t3") and
                  raw_path.parent.resolve() == (root / "pharlap_raw").resolve(),
                  f"A2 coefficient used non-training raw: {name}")
            check(experiment.sha(raw_path) == evidence["raw_sha256"],
                  f"A2 coefficient training raw drift: {name}")
            selected.append(name)
            raw = loadmat(raw_path, simplify_cells=True)
            grid_path = Path(raw["job"]["grid_mat"])
            check(experiment.sha(grid_path) == evidence["grid_sha256"],
                  f"A2 coefficient training grid drift: {name}")
            grid = loadmat(grid_path, simplify_cells=True)
            h = np.asarray(raw["ray_path_data"]["height"], float).reshape(-1)
            axis = np.asarray(grid["height_km"], float).reshape(-1)
            profile = np.asarray(grid["ne_profile_m3"], float).reshape(-1)
            independent_ne = np.interp(h, axis, profile, left=np.nan, right=np.nan)
            state = raw["ray_state_vec"]
            p2 = sum(np.asarray(state[key], float).reshape(-1) ** 2
                     for key in ("dir_x", "dir_y", "dir_z"))
            delta_ne = independent_ne - independent_ne[0]
            mask = np.isfinite(independent_ne) & (np.abs(delta_ne) >= 1e8)
            check(np.count_nonzero(mask) == evidence["sample_count"],
                  f"A2 coefficient sample selection drift: {name}")
            ceff = -(10e6) ** 2 * (p2[mask] - p2[0]) / delta_ne[mask]
            check(abs(float(np.median(ceff)) - evidence["C_eff_median"]) < 1e-8,
                  f"A2 coefficient median does not reproduce: {name}")
        check(len(selected) == 8 and len(set(selected)) == 8,
              "A2 coefficient must use eight distinct 10MHz training rays")
        check(abs(float(np.median([item["C_eff_median"] for item in freeze["training"]]))
                  - float(freeze["identified_K2"])) < 1e-8,
              "A2 identified coefficient selection mismatch")
        id_manifest = json.loads((identified / "input_manifest.json").read_text(encoding="utf-8"))
        check(experiment.sha(freeze_path) == id_manifest["coefficient_freeze_sha256"] and
              experiment.sha(identified / "input.mat") == id_manifest["input_mat_sha256"],
              "A2 identified input/freeze drift")
        id_jobs = _a2_jobs(identified / "input.mat")
        check(len(id_jobs) == 6 and all(float(job["frequency_MHz"]) == 10.0 and
              float(job["launch_k2"]) == float(freeze["identified_K2"])
              for job in id_jobs), "A2 identified launch used holdout frequency or wrong K2")
        id_summary = json.loads((identified / "pharlap_summary.json").read_text(encoding="utf-8"))
        check(id_summary.get("status") == "completed_a2_export" and
              len(id_summary["records"]) == 6 and
              {record["job_id"] for record in id_summary["records"]} ==
              {job["job_id"] for job in id_jobs},
              "A2 identified PH export incomplete")
        for job in id_jobs:
            raw_path = identified / "pharlap_raw" / f"{job['job_id']}.mat"
            raw = loadmat(raw_path, simplify_cells=True)
            check(raw["job_sha256"] == job["job_sha256"] and
                  raw["grid_sha256"] == job["grid_sha256"],
                  f"A2 identified PH raw hash mismatch: {job['job_id']}")
            state = raw["ray_state_vec"]
            check(state_fields.issubset(state),
                  f"A2 identified PH state incomplete: {job['job_id']}")
            initial = raw["initial"]
            if isinstance(initial, dict):
                first = np.array([np.asarray(state[key], float).reshape(-1)[0]
                                  for key in ("dir_x", "dir_y", "dir_z")])
                injected = np.array([float(initial[key]) for key in
                                     ("dir_x", "dir_y", "dir_z")])
                check(np.linalg.norm(first - injected) < 1e-10,
                      f"A2 identified PH first state changed: {job['job_id']}")
        internal_grid = Path(id_jobs[0]["grid_mat"])
        id_readback = identified / "pharlap_raw" / (
            f"grid_readback_{id_jobs[0]['grid_sha256'][:16]}.mat")
        if id_readback.is_file():
            original = loadmat(internal_grid, simplify_cells=True)
            readback = loadmat(id_readback, simplify_cells=True)
            check(all(key in readback and np.array_equal(value, readback[key])
                      for key, value in original.items() if not key.startswith("__")),
                  "A2 identified MATLAB grid readback mismatch")
        else:
            failures.append("A2 identified MATLAB grid readback missing")
    holdout_jobs, holdout_failures = _verify_a2_identified_holdout(
        root, jobs, freeze)
    failures.extend(holdout_failures)
    failures.extend(_verify_a2_ah_files(root, jobs, id_jobs, holdout_jobs))
    failures.extend(_verify_a2_formal_summaries(root))
    spatial_failures, _ = _verify_a2_spatial_fd(root)
    failures.extend(spatial_failures)
    failures.extend(_verify_a2_cause_matrix(root))
    failures.extend(_verify_a2_production_identity(root))
    return failures


def verify_p0(run: Path) -> list[str]:
    """Read-only freeze/fixture check; no PH or AH solver is executed."""
    failures = []

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    manifest = experiment.read_json(run / "p0_manifest.json")
    result = experiment.read_json(run / "p0_checks.json")
    inputs = experiment.read_json(run / "p0_inputs.json")
    parent_audit = experiment.read_json(run / "p0_parent_audit.json")
    check(parent_audit == experiment.p0_a2_audit(Path(manifest["parent_a2"])),
          "P0 parent raw-state audit not reproduced")
    for key, path in (("runner", Path(experiment.__file__)), ("verifier", Path(__file__))):
        check(experiment.sha(path) == manifest["sources"][key]["sha256"],
              f"P0 executing tool differs from frozen source: {key}")
    # Independently reproduce saved A2 increments with the existing raw-data verifier.
    parent = Path(manifest["parent_a2"])
    for folder, filename, ah_dir in (
        (parent, "a2_base_formal_identified_hermite_smooth_fd0p125.json",
         parent / "ah/formal/identified/hermite/fd0.125"),
        (parent / "identified_holdout", "a2_identified_holdout_formal_identified_hermite_affine.json",
         parent / "identified_holdout/ah/formal/identified/hermite")):
        summary = experiment.read_json(folder / "analysis" / filename)
        for gate in summary["gates"]:
            rows = sorted([row for row in summary["rows"]
                           if (f"{row['case_id']}|f{row['frequency_MHz']:g}|"
                               f"e{row['elevation_deg']:.14g}") == gate["case"]],
                          key=lambda row: row["level"])
            traces = [_a2_native_trace(folder, ah_dir, row) for row in rows]
            for row, trace in zip(rows, traces):
                fixed = row["stage"] in ("internal", "affine")
                end = 100.0 if fixed else min(trace["ph_s"][-1], trace["ah_s"][-1])
                target = np.linspace(0.0, end, 201 if fixed else 501)
                delta = np.linalg.norm(_a2_sample(trace["ph_s"], trace["ph_xyz"], target)
                    - _a2_sample(trace["ah_s"], trace["ah_xyz"], target), axis=1) * 1000
                ph_end = (_a2_sample(trace["ph_s"], trace["ph_xyz"], np.array([100.0]))[0]
                          if fixed else trace["ph_xyz"][-1])
                ph_group = (float(_a2_sample(trace["ph_s"], trace["ph_group"], np.array([100.0]))[0])
                            if fixed else float(trace["ph_group"][-1]))
                ah_group = float(trace["ah_record"]["group_path_km_native"])
                actual = {"position_rms_m": float(np.sqrt(np.mean(delta**2))),
                          "position_max_m": float(np.max(delta)),
                          "position_endpoint_m": float(np.linalg.norm(ph_end - trace["ah_xyz"][-1]) * 1000),
                          "group_relative_error": abs(ah_group - ph_group) / abs(ph_group)}
                for metric, value in actual.items():
                    check(np.isclose(row[metric], value, rtol=1e-8, atol=1e-20),
                          f"P0 raw residual differs: {row['job_id']}/{metric}")
            recalculated = _a2_gate_from_traces(rows, traces, rows)
            for metric, values in gate["numeric"]["checks"].items():
                for field in ("delta12_ph", "delta12_ah", "delta23_ph", "delta23_ah"):
                    check(np.isclose(values[field], recalculated["checks"][metric][field],
                                     rtol=1e-8, atol=1e-20),
                          f"P0 raw increment differs: {gate['case']}/{metric}/{field}")
    for key, record in manifest["sources"].items():
        relative = record["snapshot"]
        path = run / relative if relative else Path(record["path"])
        if relative and not path.resolve().is_relative_to(run.resolve()):
            failures.append(f"P0 snapshot escapes run: {key}")
            continue
        check(path.is_file() and experiment.sha(path) == record["sha256"],
              f"P0 frozen source hash differs: {key}")
    for key in ("runner", "verifier", "matlab_source"):
        check(bool(manifest["sources"][key]["snapshot"]), f"P0 snapshot absent: {key}")
    check(manifest["preparation_snapshot_before_inputs"] is True,
          "P0 preparation not declared frozen first")
    check((run / "p0_manifest.json").stat().st_mtime_ns <= (run / "p0_inputs.json").stat().st_mtime_ns,
          "P0 freeze manifest written after inputs")
    check(experiment.sha(run / "p0_inputs.json") == result["input_sha256"], "P0 input hash differs")
    for name, curved in (("straight", False), ("curved", True)):
        expected = experiment._p0_fixture(curved)
        for key, value in expected.items():
            check(np.array_equal(np.asarray(inputs[name][key]), np.asarray(value))
                  if isinstance(value, np.ndarray) else inputs[name][key] == value,
                  f"P0 analytic fixture changed: {name}/{key}")
    fresh = experiment.p0_checks()
    check(result["checks"] == fresh["checks"], "P0 fixture checks not reproduced")
    check(result["evidence"] == fresh["evidence"], "P0 numeric evidence not reproduced")
    check(all(isinstance(v, bool) and v for v in result["checks"].values()),
          "P0 check flags must be actual true JSON booleans")
    for key in ("passed", "evidence_integrity_passed", "synthetic_numeric_checks_passed"):
        check(result[key] is True, f"P0 failed: {key}")
    check(result["physical_rays_validated_in_p0"] is False, "P0 claims physical ray validation")
    check(result["parent_a2_global_passed"] is False, "P0 conceals historical A2 failure")
    check(manifest["production_modified"] is False and manifest["pharlap_executed"] is False,
          "P0 stage boundary changed")
    # Independently apply the per-solver rule to the stored counterexample.
    example = result["evidence"]["single_solver_nonmonotonic"]
    check(example["d23"]["ah"]["rms_m"] > example["d12"]["ah"]["rms_m"],
          "P0 counterexample no longer contains AH nonconvergence")
    check(example["decision"]["passed"] is False, "P0 accepts one-sided nonconvergence")
    plateau = result["evidence"]["method_difference_plateau"]
    check(plateau["metrics"]["rms_m"] > 0 and plateau["decision"]["passed"] is False,
          "P0 accepts an unexplained plateau")
    check(plateau["decision"]["checks"]["rms_m"]["three_times_budget"] == 0.0,
          "P0 includes an unverified interpolation-method budget")
    return failures


def _p1_ph_observations(raw: dict) -> dict:
    """Keep both native coordinates; equal array lengths do not establish co-location."""
    state, path = raw["ray_state_vec"], raw["ray_path_data"]
    vector = lambda obj, key: np.asarray(obj[key], dtype=float).reshape(-1)
    xyz = np.column_stack([vector(state, key) for key in ("pos_x", "pos_y", "pos_z")]) / 1000.0
    path_xyz = np.asarray(experiment.geocentric_to_ecef(
        vector(path, "lat"), vector(path, "lon"), vector(path, "height"), experiment.R_KM))
    state_s, path_s = vector(state, "geometrical_path") / 1000.0, vector(path, "geometric_distance")
    state_g, path_g = vector(state, "group_path") / 1000.0, vector(path, "group_range")
    p = np.column_stack([vector(state, key) for key in ("dir_x", "dir_y", "dir_z")])
    observed = {"state_points": len(state_s), "path_points": len(path_s),
                "state_end_height_m": float((np.linalg.norm(xyz[-1]) - experiment.R_KM) * 1000),
                "path_end_height_m": float(vector(path, "height")[-1] * 1000),
                "state_s_end_km": float(state_s[-1]), "path_s_end_km": float(path_s[-1]),
                "state_group_end_km": float(state_g[-1]), "path_group_end_km": float(path_g[-1]),
                "initial_xyz_km": xyz[0].tolist(), "initial_p": p[0].tolist(),
                "ray_label": int(np.asarray(raw["ray_data"]["ray_label"]).reshape(-1)[-1]),
                "finite": bool(all(np.all(np.isfinite(a)) for a in
                                   (xyz, path_xyz, state_s, path_s, state_g, path_g, p)))}
    paired = xyz.shape == path_xyz.shape and state_s.shape == path_s.shape
    observed["equal_lengths"] = bool(paired)
    if paired:
        observed.update({"path_state_max_m": float(np.max(np.linalg.norm(path_xyz - xyz, axis=1)) * 1000),
                         "path_state_s_max_m": float(np.max(np.abs(path_s - state_s)) * 1000),
                         "path_state_group_max_m": float(np.max(np.abs(path_g - state_g)) * 1000)})
    observed["native_coincident_5cm"] = bool(paired and observed["finite"] and all(
        observed[key] <= .05 for key in ("path_state_max_m", "path_state_s_max_m", "path_state_group_max_m")))
    return observed


def verify_p1(run: Path, report: dict | None = None) -> list[str]:
    """Audit immutable inputs and raw provenance, without executing either solver.

    Provenance success is independent of physical acceptance and output co-location.
    """
    failures, observations, physics, ah_status, repeat_checks = [], [], [], [], []
    result = report if report is not None else {}

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    def hashed(path: Path, digest: str, label: str) -> bool:
        ok = path.is_file() and experiment.sha(path) == digest
        check(ok, "P1 SHA mismatch/missing: " + label)
        return ok

    def same_job(actual: dict, expected: dict, label: str) -> None:
        check(set(actual) == set(expected), "P1 job fields differ: " + label)
        for key in expected:
            check(key in actual and np.array_equal(np.asarray(actual[key]).reshape(-1),
                  np.asarray(expected[key]).reshape(-1)), f"P1 job differs: {label}:{key}")

    try:
        manifest = experiment.read_json(run / "p1_manifest.json")
        check(manifest["schema"] == "c1m-p1-v1", "P1 manifest schema mismatch")
        check(manifest["production_modified"] is False, "P1 production-modified declaration differs")
        root = run if (run / "sources").is_dir() else run.parent
        sources = manifest["sources"]
        required = {"runner", "verifier", "matlab_source", "ah_solver", "ah_dispersion",
                    "ah_provider", "config", "pharlap_mex", "p0_verification", "p0_manifest"}
        check(required <= set(sources), "P1 required source inventory incomplete")
        for key, record in sources.items():
            snapshot = record.get("snapshot")
            target = root / snapshot if snapshot else Path(record["path"])
            if snapshot:
                check(target.resolve().is_relative_to(root.resolve()), "P1 snapshot outside run: " + key)
            hashed(target, record["sha256"], key)
        # The executing tools and production implementation must still match the freeze.
        for key in ("runner", "verifier", "matlab_source", "ah_solver", "ah_dispersion",
                    "ah_provider", "config", "utils", "pharlap_mex"):
            if key in sources:
                target = {"runner": Path(experiment.__file__), "verifier": Path(__file__)}.get(
                    key, Path(sources[key]["path"]))
                hashed(target, sources[key]["sha256"], "executing/current " + key)
        parent = experiment.OUT / manifest["parent"]
        old = experiment.read_json(parent / "p0_manifest.json")
        old_result = experiment.read_json(parent / "verification.json")
        check(old_result.get("p0_passed") is True and old_result.get("physical_consistency_validated") is False,
              "P1 parent P0 verification semantics differ")
        for key, record in old["sources"].items():
            target = parent / record["snapshot"] if record.get("snapshot") else Path(record["path"])
            hashed(target, record["sha256"], "parent frozen " + key)
            if key in ("ah_solver", "ah_dispersion", "ah_provider", "config", "pharlap_mex", "matlab_executable"):
                check(key in sources and sources[key]["sha256"] == record["sha256"],
                      "P1 production changed since P0: " + key)
        job_record = experiment.read_json(run / "jobs.json")
        hashed(run / "input.mat", job_record["input_sha256"], "input.mat")
        jobs = job_record["jobs"]
        by_id = {job["job_id"]: job for job in jobs}
        check(len(by_id) == len(jobs) and bool(jobs), "P1 job IDs empty/duplicated")
        input_jobs = _a2_jobs(run / "input.mat")
        check(len(input_jobs) == len(jobs), "P1 MAT job count differs")
        for expected, actual in zip(jobs, input_jobs):
            same_job(actual, expected, expected["job_id"])
            canonical = {key: value for key, value in expected.items() if key != "job_sha256"}
            digest = hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode("utf-8")).hexdigest()
            check(digest == expected["job_sha256"], "P1 job SHA differs: " + expected["job_id"])
        check((root / "p1_manifest.json").stat().st_mtime_ns <= (run / "input.mat").stat().st_mtime_ns,
              "P1 sources were not frozen before input generation")
        grids = experiment.read_json(root / "grids.json")["grids"]
        by_grid = {record["sha256"]: record for record in grids}
        check(len({record["case_id"] for record in grids}) == len(grids), "P1 grid case IDs duplicated")
        for record in grids:
            target = Path(record["path"])
            hashed(target, record["sha256"], "grid " + record["case_id"])
            grid = experiment.mat(target)
            for key, digest in record["array_sha256"].items():
                check(key in grid and experiment.array_sha(np.asarray(grid[key])) == digest,
                      "P1 grid array SHA differs: " + record["case_id"] + ":" + key)
            for key in ("Bx_T", "By_T", "Bz_T", "collision_freq_Hz"):
                check(np.all(np.asarray(grid[key]) == 0), "P1 nonzero B/collision: " + record["case_id"])
        raw_dir = run / "pharlap_raw"
        summary = experiment.read_json(run / "pharlap_summary.json")
        records = summary["records"]
        if isinstance(records, dict):
            records = [records]
        check(len(records) == len(jobs) and {r["job_id"] for r in records} == set(by_id),
              "P1 PH summary job IDs differ")
        raw_ph = {}
        for job in jobs:
            job_id = job["job_id"]
            check(job["grid_sha256"] in by_grid, "P1 job references unfrozen grid: " + job_id)
            check(any(record["sha256"] == job["grid_sha256"] and
                      Path(job["grid_mat"]).resolve() == Path(record["path"]).resolve() for record in grids),
                  "P1 job grid path differs: " + job_id)
            path = raw_dir / (job_id + ".mat")
            if not path.is_file():
                failures.append("P1 PH raw absent: " + job_id)
                continue
            raw = experiment.mat(path)
            raw_ph[job_id] = raw
            same_job(raw["job"], job, "PH " + job_id)
            for key, expected in (("job_sha256", job["job_sha256"]), ("grid_sha256", job["grid_sha256"]),
                                  ("exporter_sha256", sources["matlab_source"]["sha256"]),
                                  ("mex_sha256", sources["pharlap_mex"]["sha256"])):
                check(raw.get(key) == expected, "P1 PH raw identity differs: " + job_id + ":" + key)
            readback = raw_dir / ("grid_readback_" + job["grid_sha256"][:16] + ".mat")
            if job["grid_sha256"] not in {r.get("grid_sha256") for r in observations}:
                read = experiment.mat(readback)
                for key, digest in by_grid[job["grid_sha256"]]["array_sha256"].items():
                    check(key in read and experiment.array_sha(np.asarray(read[key])) == digest,
                          "P1 PH readback array differs: " + job_id + ":" + key)
            observation = _p1_ph_observations(raw)
            observation.update({"job_id": job_id, "grid_sha256": job["grid_sha256"],
                                "raw_sha256": experiment.sha(path)})
            observations.append(observation)
            if not observation["native_coincident_5cm"]:
                physics.append({"job_id": job_id, "finding": "PH path/state not verified co-located at 5cm diagnostic scale"})
            if "initial_p" in job:
                observation["initial_p_requested_delta"] = float(np.linalg.norm(
                    np.asarray(observation["initial_p"]) - np.asarray(job["initial_p"])))
                observation["initial_xyz_requested_delta_m"] = float(np.linalg.norm(
                    np.asarray(observation["initial_xyz_km"]) - np.asarray(job["initial_xyz_km"])) * 1000)
                if observation["initial_p_requested_delta"] > 1e-10 or observation["initial_xyz_requested_delta_m"] > 1e-5:
                    physics.append({"job_id": job_id, "finding": "PH did not retain requested restart state"})
            if "parent_raw" in job:
                hashed(Path(job["parent_raw"]), job["parent_raw_sha256"], "restart parent " + job_id)
        for job in jobs:
            if job.get("role") != "cache_control" or job["job_id"] not in raw_ph:
                continue
            keys = ("grid_sha256", "frequency_MHz", "elevation_deg", "origin_height_km", "level")
            candidates = [j for j in jobs if j.get("role") != "cache_control" and all(j[k] == job[k] for k in keys)]
            check(bool(candidates), "P1 cache reference absent: " + job["job_id"])
            if not candidates or candidates[0]["job_id"] not in raw_ph:
                continue
            first, second = raw_ph[candidates[0]["job_id"]], raw_ph[job["job_id"]]
            equal = all(set(first[layer]) == set(second[layer]) and all(np.array_equal(
                np.asarray(first[layer][key]), np.asarray(second[layer][key]), equal_nan=True)
                for key in first[layer]) for layer in ("ray_path_data", "ray_state_vec"))
            repeat_checks.append({"job_id": job["job_id"], "reference_job_id": candidates[0]["job_id"], "native_outputs_identical": bool(equal)})
            if not equal:
                physics.append({"job_id": job["job_id"], "finding": "PH cache/repeat native outputs changed"})
        ah_sources = {key: sources[value]["sha256"] for key, value in
                      {"runner": "runner", "config": "config", "dispersion": "ah_dispersion",
                       "tracer": "ah_solver", "provider": "ah_provider"}.items()}
        seen_ah = set()
        for path in sorted((run / "ah").glob("fd*/*.json")):
            raw = experiment.read_json(path)
            job_id = raw["job_id"]
            check(job_id in by_id, "P1 AH unexpected job: " + job_id)
            if job_id not in by_id:
                continue
            seen_ah.add(job_id)
            job = by_id[job_id]
            check(raw["job_sha256"] == job["job_sha256"] and raw["grid_sha256"] == job["grid_sha256"] and
                  raw["source_sha256"] == ah_sources, "P1 AH source/job/grid identity differs: " + str(path))
            check(raw["coefficient_K2"] == manifest["coefficient"] and raw["interpolation"] == "hermite",
                  "P1 AH coefficient/interpolation differs: " + str(path))
            if raw.get("raw_npz_sha256"):
                target = path.with_suffix(".npz")
                if hashed(target, raw["raw_npz_sha256"], "AH native " + str(path)):
                    with np.load(target, allow_pickle=False) as native:
                        check(native["xyz_km"].shape[0] == raw["point_count"], "P1 AH native point count differs: " + str(path))
            else:
                check(raw.get("status") == "failed", "P1 AH native SHA absent without recorded failure: " + str(path))
            ah_status.append({"job_id": job_id, "fd_km": raw["spatial_fd_base_km"], "status": raw["status"],
                              "termination_status": raw.get("termination_status"), "error": raw.get("error"),
                              "failure_reason": raw.get("failure_reason")})
            if raw.get("status") != "completed":
                physics.append({"job_id": job_id, "finding": "AH physical/integration result not completed",
                                "status": raw.get("status"), "reason": raw.get("failure_reason", raw.get("error"))})
        check(seen_ah == set(by_id), "P1 AH job inventory incomplete")
    except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        failures.append("P1 unreadable/inconsistent evidence: " + type(exc).__name__ + ": " + str(exc))
    result.update({"schema": "c1m-p1-readonly-verification-v1", "source_verified": not failures,
                   "source_failures": failures, "physics_passed": False,
                   "physics_acceptance": "not established by provenance verification",
                   "ph_output_observations": observations, "physics_findings": physics,
                   "cache_repeat_checks": repeat_checks, "ah_status": ah_status})
    print(f"P1 provenance: {len(failures)} error(s); PH output findings: {len(physics)}; physical acceptance pending")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id")
    parser.add_argument("--stage", choices=("a1", "a2", "p0", "p1"), default="a1")
    parser.add_argument("--p1-dataset", choices=("base", "restart"), default="base")
    parser.add_argument("--report-json", type=Path, help="Save the P1 read-only audit report")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.run_id:
        if Path(args.run_id).name != args.run_id or args.run_id in {'.','..'}:
            parser.error("run-id must be a directory name")
        if args.stage == "p1":
            report = {}
            target = experiment.OUT / args.run_id
            if args.p1_dataset == "restart":
                target = target / "restart"
            failures = verify_p1(target, report)
            if args.report_json:
                experiment.write_json(args.report_json, report)
        else:
            failures = {"a1": verify, "a2": verify_a2, "p0": verify_p0}[args.stage](
                experiment.OUT / args.run_id)
        for failure in failures:
            print("FAIL: " + failure)
        print(f"C1M {args.stage} read-only verification: {len(failures)} failure(s)")
        return int(bool(failures))
    if not args.self_check:
        parser.error("use --run-id or --self-check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
