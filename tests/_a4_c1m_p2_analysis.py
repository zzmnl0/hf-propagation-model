"""P2 preparation and version-bound dispatch; preserve the historical P1/A2 data."""
from __future__ import annotations
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
sys.dont_write_bytecode = True
import numpy as np
from scipy.io import loadmat
import _run_a4_c1m_attribution as common

CODE = Path(__file__).resolve().parents[1]
BASE = CODE / "output2/a4_pharlap/c1/c1m"
A2 = BASE / "a2_20260928_v1/a2"
P1 = BASE / "p1_20261007_v1"
ODE_SOURCE = P1 / "analysis_source_versions/p1_ode_mechanism_v1/independent_mechanism.py"
CORE_SOURCE = CODE / "tests/_a4_c1m_p2_analysis_core_v1.py"
LIMITS = {"rms_m": 10.0, "max_m": 50.0, "landing_m": 50.0, "group_rel": 1e-5}
ARC_KM = 100.0


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _records(path: Path) -> dict[str, dict]:
    value = _json(path)
    records = value.get("records", value.get("files", []))
    if isinstance(records, dict):
        records = [{"path": k, "sha256": v if isinstance(v, str) else v["sha256"]} for k, v in records.items()]
    indexed = {}
    for record in records:
        target = Path(record["path"])
        if not target.is_absolute():
            target = path.parent / target
        indexed[str(target.resolve()).casefold()] = record
    return indexed


def _historical_record(path: Path, records: dict[str, dict]) -> dict:
    key = str(path.resolve()).casefold()
    if key not in records:
        raise ValueError("Historical delivery has no record: " + str(path))
    expected = records[key]["sha256"]
    if common.sha(path) != expected:
        raise ValueError("Historical input changed: " + str(path))
    return {"path": str(path.resolve()), "sha256": expected}


def _snapshot(run: Path, key: str, path: Path) -> dict:
    record = {"path": str(path.resolve()), "sha256": common.sha(path), "snapshot": None}
    if path.suffix.lower() not in (".exe", ".mexw64", ".mat"):
        destination = run / "sources" / (key + path.suffix)
        if destination.exists():
            raise FileExistsError(destination)
        shutil.copy2(path, destination)
        record["snapshot"] = str(destination.relative_to(run))
        if common.sha(destination) != record["sha256"]:
            raise ValueError("Source snapshot mismatch")
    return record


def read_ph(path: Path, job: dict) -> dict:
    """Augment the existing A2 native reader; old A2 has no exporter hash field."""
    raw = common.mat(path)
    trace = common._a2_ph_trace(path, job)
    state, output = raw["ray_state_vec"], raw["ray_path_data"]
    trace["p"] = np.column_stack([common.vector(state, k) for k in ("dir_x", "dir_y", "dir_z")])
    trace["ne_out"] = trace["ne_m3"]
    trace["path_height"] = trace["height"].copy()
    trace["height"] = np.linalg.norm(trace["xyz"], axis=1) - common.R_KM
    trace["n_out"] = common.vector(output, "refractive_index")
    trace["ng_out"] = common.vector(output, "group_refractive_index")
    trace["path_xyz"] = np.asarray(common.geocentric_to_ecef(common.vector(output, "lat"), common.vector(output, "lon"), trace["path_height"], common.R_KM))
    trace["termination"] = "ground" if trace["ray_label"] == 1. else "other"
    trace["exporter_identity_available"] = "exporter_sha256" in raw
    for key in ("geomag_x", "geomag_y", "geomag_z", "collision_frequency"):
        if np.any(common.vector(output, key) != 0.):
            raise ValueError("Historical internal ray is not B0/collisionless")
    if any(np.any(~np.isfinite(trace[key])) for key in ("s", "xyz", "p", "group", "ne_out", "n_out", "ng_out")):
        raise ValueError("Invalid native internal array")
    return trace


def _source_paths() -> dict[str, Path]:
    paths = {"runner": CODE / "tests/_run_a4_c1m_attribution.py", "p2_prepare_dispatcher": Path(__file__),
        "p2_plan": CODE / "word/a4_c1m_p2_implementation_plan_20261008.txt",
        "matlab_exporter": CODE / "tests/matlab/a4_c1m_uniform_export.m",
        "ah_solver": CODE / "models/ray_tracer_ah_3d_ecef.py", "ah_dispersion": CODE / "models/dispersion_ah.py",
        "ah_provider": CODE / "models/ne_provider_geocentric.py", "config": CODE / "config.py", "utils": CODE / "utils.py",
        "python_executable": Path(sys.executable), "matlab_executable": common.MATLAB_EXE, "pharlap_mex": common.MEX,
        "historical_P1_ODE": ODE_SOURCE, "historical_A2_radial_check": A2 / "spherical_check.py"}
    for key, path in {"p2_exit": CODE / "tests/_a4_c1m_p2_exit.py",
                      "p2_exit_matlab": CODE / "tests/matlab/a4_c1m_p2_exit_export.m"}.items():
        if path.is_file():
            paths[key] = path
    return paths


def prepare(run: Path) -> None:
    run = Path(run).resolve()
    if (run / "p2_manifest.json").exists():
        raise FileExistsError("P2 manifest already exists")
    if run.exists():
        allowed = {"exit", "exit_prepare_failed_v1"}
        unexpected = [p.name for p in run.iterdir() if p.name not in allowed]
        if unexpected:
            raise ValueError("P2 contains unrelated or internal outputs: " + repr(unexpected))
        if (run / "exit").exists() and not (run / "exit/manifest.json").is_file():
            raise ValueError("Existing exit component has no valid manifest")
        if (run / "exit_prepare_failed_v1").exists() and not (run / "exit_prepare_failed_v1/failure.json").is_file():
            raise ValueError("Unexpected failed-component directory")
    run.mkdir(parents=True, exist_ok=True)
    (run / "sources").mkdir()
    (run / "inputs").mkdir()
    a2_delivery = _records(A2 / "sha256_manifest.json")
    p1_delivery = _records(P1 / "delivery_sha256.json")
    p1_manifest = _json(P1 / "p1_manifest.json")
    history = {}
    for path in (P1 / "p1_manifest.json", P1 / "sources/runner.py", ODE_SOURCE):
        history[str(path)] = _historical_record(path, p1_delivery)
    if common.sha(P1 / "sources/runner.py") != p1_manifest["sources"]["runner"]["sha256"]:
        raise ValueError("P1 historical runner snapshot differs from its freeze")
    history[str(P1 / "delivery_sha256.json")] = {"path": str(P1 / "delivery_sha256.json"), "sha256": common.sha(P1 / "delivery_sha256.json")}
    sources = {key: _snapshot(run, key, path) for key, path in _source_paths().items()}
    for key in ("ah_solver", "ah_dispersion", "ah_provider", "config", "utils", "pharlap_mex"):
        if sources[key]["sha256"] != p1_manifest["sources"][key]["sha256"]:
            raise ValueError("Production/MEX changed since P1: " + key)
    selected, auxiliary, inputs = [], [], {}
    for dataset, folder, formal in (("identified", A2 / "identified", True), ("common", A2, False)):
        input_manifest = _json(folder / "input_manifest.json")
        jobs = [j for j in input_manifest["jobs"] if j["stage"] == "internal"]
        if len(jobs) != 6 or sorted({j["elevation_deg"] for j in jobs}) != [22.0, 39.5]:
            raise ValueError("Expected both historical internal angles at all three levels")
        for path in (folder / "input_manifest.json", folder / "input.mat"):
            inputs[str(path)] = _historical_record(path, a2_delivery)
        input_sha = common.sha(folder / "input.mat")
        if input_sha != input_manifest["input_mat_sha256"]:
            raise ValueError("Historical input MAT changed")
        for job in jobs:
            grid_path = Path(job["grid_mat"])
            inputs[str(grid_path)] = _historical_record(grid_path, a2_delivery)
            if common.sha(grid_path) != job["grid_sha256"]:
                raise ValueError("A2 internal grid differs job identity")
            grid = loadmat(grid_path, variable_names=["height_km", "ne_profile_m3"], simplify_cells=True)
            ph_path = folder / "pharlap_raw" / (job["job_id"] + ".mat")
            inputs[str(ph_path)] = _historical_record(ph_path, a2_delivery)
            ph = read_ph(ph_path, job)
            raw = common.mat(ph_path)
            if any(str(raw[k]) != str(job[k]) for k in ("job_sha256", "grid_sha256")):
                raise ValueError("Original PH metadata mismatch")
            encoded = {k: v for k, v in job.items() if k != "job_sha256"}
            canonical = hashlib.sha256(json.dumps(encoded, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
            if canonical != job["job_sha256"]:
                raise ValueError("Original canonical job identity mismatch")
            ne0 = float(np.interp(100., grid["height_km"], grid["ne_profile_m3"]))
            expected_xyz = np.asarray(common.geocentric_to_ecef(30., 120., 100., common.R_KM))
            expected_p = np.sqrt(1. - job["launch_k2"] * ne0 / 1e14) * common.launch_direction_ecef(30., 120., job["elevation_deg"], 0., common.R_KM)
            if np.linalg.norm(ph["xyz"][0] - expected_xyz) > 1e-8 or np.linalg.norm(ph["p"][0] - expected_p) > 1e-12:
                raise ValueError("Original common initial state mismatch")
            if ph["s"][-1] < ARC_KM:
                raise ValueError("Original PH missing the complete 100-km interval")
            ah_folder = (folder / "ah/formal/identified/hermite") if formal else (folder / "ah/formal/common/hermite")
            ah_meta = ah_folder / (job["job_id"] + ".json")
            ah_npz = ah_meta.with_suffix(".npz")
            for path in (ah_meta, ah_npz):
                inputs[str(path)] = _historical_record(path, a2_delivery)
            ah = _json(ah_meta)
            if ah["job_sha256"] != job["job_sha256"] or ah["grid_sha256"] != job["grid_sha256"] or ah["status"] != "completed":
                raise ValueError("Original AH identity/completion mismatch")
            if ah["raw_npz_sha256"] != common.sha(ah_npz) or abs(ah["geom_path_km"] - ARC_KM) > 1e-10:
                raise ValueError("Original AH raw or 100-km endpoint mismatch")
            item = {"dataset": dataset, "formal": formal, "original_job": job, "PH_raw": str(ph_path), "AH_raw": str(ah_npz),
                "AH_meta": str(ah_meta), "original_fd_km": .5,
                "initial_xyz_km": ph["xyz"][0].tolist(), "initial_p": ph["p"][0].tolist(),
                "input_Ne100_m3": ne0, "output_Ne0_m3": float(ph["ne_out"][0]), "common_arc_km": ARC_KM}
            (selected if formal else auxiliary).append(item)
    for path in (A2 / "identified/coefficient_freeze.json", A2 / "identified/analysis/a2_identified_formal_identified_hermite_internal.json",
                 A2 / "sources/ah/5ab173b72dba2b16/runner.py", A2 / "sources/ah/5ab173b72dba2b16/tracer.py"):
        inputs[str(path)] = _historical_record(path, a2_delivery)
    coefficient = float(_json(A2 / "identified/coefficient_freeze.json")["identified_K2"])
    _write(run / "inputs/internal_jobs.json", {"primary": selected, "auxiliary": auxiliary})
    shutil.copy2(A2 / "grids/internal_common.mat", run / "inputs/internal_common.mat")
    if common.sha(run / "inputs/internal_common.mat") != selected[0]["original_job"]["grid_sha256"]:
        raise ValueError("Copied internal grid differs")
    if not (run / "exit/manifest.json").exists():
        import _a4_c1m_p2_exit as exit_component
        exit_component.prepare(run)
    exit_manifest = run / "exit/manifest.json"
    exit_reference = {"manifest": str(exit_manifest), "sha256": common.sha(exit_manifest), "already_prepared": True,
        "component_evidence": {str(p): common.sha(p) for p in sorted((run / "exit").rglob("*.json"))}}
    failed = run / "exit_prepare_failed_v1/failure.json"
    manifest = {"schema": "a4-c1m-p2-prepare-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": sources, "historical_inputs": inputs, "P1_history": history,
        "prepared_inputs": {str(run / "inputs/internal_jobs.json"): common.sha(run / "inputs/internal_jobs.json"),
                            str(run / "inputs/internal_common.mat"): common.sha(run / "inputs/internal_common.mat")},
        "constants": {"C": coefficient, "coefficient_source": str(A2 / "identified/coefficient_freeze.json"), "R_km": common.R_KM,
                      "C_km_s": 299792.458}, "arc_km": ARC_KM, "limits": LIMITS,
        "PH_profiles": common.TOL[:3].tolist(), "AH_profiles": common.AH_PROFILES[:3],
        "primary_dataset": "internal_identified_22_and_39p5", "original_fd_km": .5, "candidate_fd_km": .125,
        "original_fd_evidence": "Frozen 5ab runner constructor has no FD override; same snapshot tracer defaults to 0.5 km",
        "new_AH_runs": "6 candidate FD0.125 identified jobs, isolated workers<=2 with BLAS threads=1",
        "new_PH_runs": "No internal PH rerun; original 3 legal levels reused. Exit component has its own 12 new controls.",
        "auxiliary_common_is_formal": False,
        "original_A2_provenance_complete": False,
        "original_A2_provenance_note": "P2 verifies selected historical data against archived delivery but does not close the older preparation-history gap",
        "analysis_implementation": {"path": str(CORE_SOURCE), "binding": "analysis_manifest_written_before_execution", "prepared_at_this_stage": False},
        "verifier_implementation": {"path": str(CODE / "tests/_a4_c1m_p2_verify.py"), "binding": "verifier_freeze_written_before_execution", "prepared_at_this_stage": False},
        "exit_component": exit_reference, "preserved_failed_exit": {"path": str(failed), "sha256": common.sha(failed)} if failed.exists() else None,
        "environment": {"python_executable": sys.executable, "python_version": sys.version, "numpy_version": np.__version__,
            "scipy_version": __import__("scipy").__version__, "BLAS_threads": 1},
        "production_modified": False, "physics_passed": False, "global_passed": False}
    _write(run / "p2_manifest.json", manifest)
    print("P2 prepared original identified/common inputs; primary=6 auxiliary=6; " + str(run), flush=True)


def validate_prepare(run: Path) -> dict:
    manifest = _json(run / "p2_manifest.json")
    for record in manifest["sources"].values():
        if common.sha(Path(record["path"])) != record["sha256"]:
            raise ValueError("Prepared execution source changed: " + record["path"])
        if record["snapshot"] and common.sha(run / record["snapshot"]) != record["sha256"]:
            raise ValueError("Prepared snapshot changed")
    for category in ("historical_inputs", "P1_history"):
        for record in manifest[category].values():
            if common.sha(Path(record["path"])) != record["sha256"]:
                raise ValueError("Historical input/snapshot changed: " + record["path"])
    for path, expected in manifest["prepared_inputs"].items():
        if common.sha(Path(path)) != expected:
            raise ValueError("Prepared input changed: " + path)
    if common.sha(Path(manifest["exit_component"]["manifest"])) != manifest["exit_component"]["sha256"]:
        raise ValueError("Prepared exit manifest changed")
    for path, expected in manifest["exit_component"].get("component_evidence", {}).items():
        if common.sha(Path(path)) != expected:
            raise ValueError("Frozen exit component evidence changed: " + path)
    return manifest


def analyze(run: Path) -> None:
    run = Path(run).resolve()
    manifest = validate_prepare(run)
    if not CORE_SOURCE.is_file():
        raise FileNotFoundError("Analysis core has not yet been delivered: " + str(CORE_SOURCE))
    output = run / "analysis/p2_independent"
    if output.exists():
        raise FileExistsError("Preserve previous analysis; bind a distinct analysis version")
    (output / "sources").mkdir(parents=True)
    core_snapshot = output / "sources" / CORE_SOURCE.name
    shutil.copy2(CORE_SOURCE, core_snapshot)
    sha = common.sha(CORE_SOURCE)
    _write(output / "analysis_manifest.json", {"schema": "a4-c1m-p2-analysis-source-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(), "prepared_manifest_sha256": common.sha(run / "p2_manifest.json"),
        "sources": {"core": {"path": str(CORE_SOURCE), "sha256": sha, "snapshot": str(core_snapshot.relative_to(output))}},
        "C": manifest["constants"]["C"], "production_modified": False, "fitted_parameters": False,
        "reconstruction_method_difference_in_budget": False})
    spec = importlib.util.spec_from_file_location("_a4_c1m_p2_bound_core", core_snapshot)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    try:
        module.analyze(run, output, manifest)
    except Exception as exc:
        _write(output / "execution_failure.json", {"type": type(exc).__name__, "message": str(exc), "outputs_preserved": True})
        raise
    if common.sha(CORE_SOURCE) != sha or common.sha(core_snapshot) != sha:
        raise ValueError("Analysis source changed during execution")
