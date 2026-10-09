"""P2 exit semantics: frozen explicit-state PH controls, no production changes."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

for _thread_var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_var] = "1"

import numpy as np
from scipy.io import loadmat, savemat

MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
PH_ROOT = Path(r"D:\pharlap\pharlap_4.7.4")
MEX = PH_ROOT / "mex" / "raytrace_3d_sp.mexw64"
PYTHON = Path(r"C:\Users\12238\.conda\envs\pytorch_cpu\python.exe")
REFERENCE = "c1_h60_d2_f10_e39p5425_t1"
TOLS = [[1e-8, .01, 10.], [1e-10, .001, 5.], [1e-12, .001, 1.]]
STATE_INPUT = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z",
               "group_path", "geom_path", "phase_path", "absorption", "indep_var", "ODE_step_size")


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def vec(struct: dict, key: str) -> np.ndarray:
    return np.asarray(struct[key], dtype=float).reshape(-1)


def read_raw(path: Path) -> dict:
    raw = loadmat(path, simplify_cells=True, mat_dtype=True)
    state, output = raw["ray_state_vec"], raw["ray_path_data"]
    xyz = np.column_stack([vec(state, key) for key in STATE_INPUT[:3]]) / 1000.
    p = np.column_stack([vec(state, key) for key in STATE_INPUT[3:6]])
    lat, lon = np.deg2rad(vec(output, "lat")), np.deg2rad(vec(output, "lon"))
    path_h = vec(output, "height")
    path_xyz = (6371. + path_h)[:, None] * np.column_stack(
        [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)])
    arrays = {"xyz": xyz, "p": p, "s": vec(state, "geometrical_path") / 1000.,
              "G": vec(state, "group_path") / 1000., "height": np.linalg.norm(xyz, axis=1)-6371.,
              "path_xyz": path_xyz, "path_h": path_h, "path_s": vec(output, "geometric_distance"),
              "path_G": vec(output, "group_range"), "Ne": vec(output, "electron_density")*1e6,
              "n": vec(output, "refractive_index"), "ng": vec(output, "group_refractive_index"),
              "indep_var": vec(state, "indep_var"), "ODE_step_size": vec(state, "ODE_step_size")}
    assert all(v.shape[0] == xyz.shape[0] for v in arrays.values()), "path/state length mismatch"
    assert all(np.all(np.isfinite(v)) for v in arrays.values()), "Nonfinite ray coordinate or state"
    for key in ("geomag_x", "geomag_y", "geomag_z", "collision_frequency"):
        assert np.all(vec(output, key)==0.), "B0/collisionless contract changed"
    arrays["raw"] = raw
    return arrays


def exit_index(data: dict) -> int:
    descending = np.sum(data["xyz"] * data["p"], axis=1) < 0.
    indices = np.flatnonzero((data["height"] < 60.) & descending)
    assert indices.size and indices[0]>0 and indices[0]+2 < len(data["s"])
    return int(indices[0])


def verify(manifest: dict, out: Path) -> dict:
    checks = {"inputs_unchanged": all(sha(Path(k))==v for k,v in manifest["inputs_sha256"].items()),
              "external_sources_unchanged": all(sha(Path(k))==v for k,v in manifest["external_sha256"].items()),
              "source_originals_unchanged": all(sha(Path(v["path"]))==v["sha256"] for v in manifest["sources"]),
              "source_snapshots_unchanged": all(sha(out/v["snapshot"])==v["sha256"] for v in manifest["sources"])}
    assert all(checks.values()), "Exit freeze changed: " + str(checks)
    return checks


def prepare(run: Path) -> None:
    run = run.resolve()
    out = run / "exit"
    out.mkdir(parents=True, exist_ok=False)
    sources = out / "sources"
    sources.mkdir()
    module = Path(__file__).resolve()
    exporter = module.parent / "matlab" / "a4_c1m_p2_exit_export.m"
    if not exporter.is_file():
        exporter = module.with_name("a4_c1m_p2_exit_export.m")
    records = []
    for source in (module, exporter):
        snapshot = sources / source.name
        shutil.copy2(source, snapshot)
        records.append({"path": str(source), "sha256": sha(source), "snapshot": str(snapshot.relative_to(out))})
    # Both preparation sources are frozen before opening any P1 input.
    write_json(out / "source_freeze.json", {"created_utc": datetime.now(timezone.utc).isoformat(), "sources": records})
    p1 = run.with_name("p1_20261007_v1")
    natural = [j for j in read_json(p1/"jobs.json")["jobs"] if j["scope"]=="full"]
    assert len(natural)==12
    reference = next(j for j in natural if j["job_id"]==REFERENCE)
    ref_raw = p1 / "pharlap_raw" / (REFERENCE+".mat")
    data = read_raw(ref_raw)
    first = exit_index(data)
    state = data["raw"]["ray_state_vec"]
    fixed = {key: float(vec(state, "geometrical_path" if key=="geom_path" else key)[first]) for key in STATE_INPUT}
    fixed["ODE_step_size"] = 1000.
    p_fixed = data["p"][first]
    assert abs(float(np.dot(p_fixed,p_fixed))-1.)>1e-6, "Reference defect must be resolved above roundoff"
    jobs = []
    for position, index, normalize in (("last_inside",first-1,False), ("first_outside",first,False),
                                        ("next_outside",first+1,False), ("first_outside_unit",first,True)):
        initial = dict(fixed)
        initial.update({key: float(vec(state,key)[index]) for key in STATE_INPUT[:3]})
        if normalize:
            initial.update({key: float(value) for key,value in zip(STATE_INPUT[3:6],p_fixed/np.linalg.norm(p_fixed))})
        for level,tol in enumerate(TOLS,1):
            job = {"job_id": f"exit_{position}_t{level}", "position_role": position,
                   "normalize_input_p": normalize, "source_job_id": REFERENCE,
                   "source_raw_mat": str(ref_raw), "source_raw_sha256": sha(ref_raw),
                   "source_position_index": index, "source_magnitude_index": first,
                   "source_elevation_deg": float(reference["elevation_deg"]),
                   "source_natural_p": data["p"][index].tolist(), "diagnostic_only": True,
                   "common_p_inside_is_not_a_physical_initial_state": position=="last_inside",
                   "common_initial_ODE_step_m": 1000.,
                   "initial": initial, "level": level, "tol": tol, "frequency_MHz": 10.,
                   "h_boundary_km": 60., "diagnostic_span_km": 10.,
                   "grid_mat": reference["grid_mat"], "grid_sha256": reference["grid_sha256"],
                   "exporter_sha256": sha(exporter), "mex_sha256": sha(MEX)}
            job["job_sha256"] = hashlib.sha256(json.dumps(job,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
            jobs.append(job)
    input_mat = out / "input.mat"
    savemat(input_mat, {"exit_jobs": np.array(jobs,dtype=object)}, do_compression=True, oned_as="row", long_field_names=True)
    write_json(out / "jobs.json", {"jobs": jobs, "natural_jobs": natural})
    roundtrip = loadmat(input_mat, simplify_cells=True, mat_dtype=True)["exit_jobs"]
    for job,loaded in zip(jobs,roundtrip):
        assert tuple(loaded["initial"])==STATE_INPUT
        assert all(float(loaded["initial"][k])==job["initial"][k] for k in STATE_INPUT)
        assert np.array_equal(loaded["tol"],job["tol"])
    inputs = [p1/"jobs.json", p1/"p1_manifest.json", input_mat, out/"jobs.json"]
    inputs.extend(p1/"pharlap_raw"/(j["job_id"]+".mat") for j in natural)
    inputs.extend(Path(j["grid_mat"]) for j in natural)
    ext = [MATLAB,MEX,PH_ROOT/"src"/"matlab"/"raytrace_3d_sp.m", PH_ROOT/"src"/"C"/"raytrace-3d-sp_matlab_wrapper.c"]
    manifest = {"schema": "c1m-p2-exit-freeze-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
                "sources": records, "inputs_sha256": {str(p):sha(p) for p in inputs},
                "external_sha256": {str(p):sha(p) for p in ext}, "matrix_count": len(jobs),
                "natural_reference_count": len(natural), "source_freeze_precedes_input_generation": True,
                "physical_acceptance_passed": False, "closed_core_location_identified": False}
    write_json(out/"manifest.json",manifest)
    expression = matlab_expression(out)
    ps_quote = lambda v: "'"+str(v).replace("'","''")+"'"
    script = "\n".join(["$ErrorActionPreference = 'Stop'", "$env:OPENBLAS_NUM_THREADS='1'",
                        "$env:OMP_NUM_THREADS='1'", "$env:MKL_NUM_THREADS='1'",
                        "& "+ps_quote(MATLAB)+" -singleCompThread -batch "+ps_quote(expression),
                        "if ($LASTEXITCODE -ne 0) { throw 'MATLAB exit export failed; partial raw records preserved' }"])+"\n"
    (out/"run_matlab.ps1").write_text(script,encoding="utf-8")
    print("P2 EXIT PREPARED 12 controls; sources frozen before input generation")


def matlab_expression(out: Path) -> str:
    quote = lambda p: "'"+str(p).replace("\\","/").replace("'","''")+"'"
    return (f"addpath({quote(out/'sources')}); addpath({quote(PH_ROOT/'mex')},'-begin'); "
            f"addpath(genpath({quote(PH_ROOT/'src'/'matlab')}),'-end'); "
            f"diary({quote(out/'matlab_diary.log')}); "
            f"a4_c1m_p2_exit_export({quote(out/'input.mat')},{quote(out/'pharlap_raw')},{quote(out/'execution_checkpoint.json')}); diary off;")


def point_row(data: dict, index: int, job_id: str, role: str) -> dict:
    pnorm = float(np.linalg.norm(data["p"][index]))
    row = {"job_id": job_id, "role": role, "index_zero_based": index,
           "s_absolute_km": float(data["s"][index]), "G_absolute_km": float(data["G"][index]),
           "state_height_km": float(data["height"][index]), "path_height_km": float(data["path_h"][index]),
           "p2": pnorm**2, "inverse_p": 1./pnorm, "Ne_m3": float(data["Ne"][index]),
           "n_path": float(data["n"][index]), "ng_path": float(data["ng"][index]),
           "path_state_position_error_m": float(np.linalg.norm(data["xyz"][index]-data["path_xyz"][index])*1000),
           "path_state_s_error_m": float((data["path_s"][index]-data["s"][index])*1000),
           "path_state_G_error_m": float((data["path_G"][index]-data["G"][index])*1000),
           "indep_var": float(data["indep_var"][index]), "ODE_step_size": float(data["ODE_step_size"][index])}
    row.update({key:float(value) for key,value in zip(STATE_INPUT[:6],np.r_[data["xyz"][index]*1000,data["p"][index]])})
    return row


def interval_row(data: dict, index: int, job_id: str, role: str) -> dict:
    ds = float(data["s"][index+1]-data["s"][index])
    dG = float(data["G"][index+1]-data["G"][index])
    chord = float(np.linalg.norm(data["xyz"][index+1]-data["xyz"][index]))
    ds_path = float(data["path_s"][index+1]-data["path_s"][index])
    dG_path = float(data["path_G"][index+1]-data["path_G"][index])
    assert ds>0 and chord>0 and ds_path>0, "Non-increasing native geometrical path"
    left, right = (1./float(np.linalg.norm(data["p"][j])) for j in (index,index+1))
    return {"job_id":job_id,"role":role,"left_index_zero_based":index,
            "left_height_km":float(data["height"][index]),"right_height_km":float(data["height"][index+1]),
            "delta_s_km":ds,"delta_G_km":dG,"chord_km":chord,"Gprime_state":dG/ds,
            "Gprime_path":dG_path/ds_path,"Gprime_chord":dG/chord,"chord_over_s":chord/ds,
            "left_inverse_p":left,"right_inverse_p":right,"Gprime_minus_one":dG/ds-1.,
            "Gprime_minus_left_inverse_p":dG/ds-left,"Gprime_minus_right_inverse_p":dG/ds-right,
            "p_change":float(np.linalg.norm(data["p"][index+1]-data["p"][index])),
            "direction_change":float(np.linalg.norm(data["p"][index+1]*right-data["p"][index]*left)),
            "delta_indep_var":float(data["indep_var"][index+1]-data["indep_var"][index])}


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(run: Path, analysis_id: str = "analysis_v1") -> None:
    out = run.resolve()/"exit"
    manifest = read_json(out/"manifest.json")
    checks = verify(manifest,out)
    analysis = out/analysis_id
    analysis.mkdir(exist_ok=False)
    matrix = read_json(out/"jobs.json")
    natural_points, natural_intervals, natural_cases = [], [], []
    for job in matrix["natural_jobs"]:
        data = read_raw(run.with_name("p1_20261007_v1")/"pharlap_raw"/(job["job_id"]+".mat"))
        for key in ("job_sha256","grid_sha256","exporter_sha256"):
            assert str(data["raw"][key])==job[key], "P1 raw identity changed: "+key
        first = exit_index(data)
        for index,role in ((first-1,"last_inside"),(first,"first_outside"),(first+1,"next_outside"),(len(data["s"])-1,"terminal")):
            natural_points.append(point_row(data,index,job["job_id"],role))
        for index in range(first-1,min(first+4,len(data["s"])-1)):
            natural_intervals.append(interval_row(data,index,job["job_id"],"natural_exit"))
        normalized = data["p"][first]/np.linalg.norm(data["p"][first])
        natural_cases.append({"job_id":job["job_id"],"level":job["level"],"first_outside_p2":float(np.dot(data["p"][first],data["p"][first])),
                              "later_unit_p2_error":float(np.max(abs(np.sum(data["p"][first+1:]**2,axis=1)-1.))),
                              "later_direction_from_first_unit_error":float(np.max(np.linalg.norm(data["p"][first+1:]-normalized,axis=1))),
                              "state_path_nonterminal_max_m":float(np.max(np.linalg.norm(data["xyz"][:-1]-data["path_xyz"][:-1],axis=1))*1000),
                              "terminal_state_path_error_m":float(np.linalg.norm(data["xyz"][-1]-data["path_xyz"][-1])*1000)})
    points, intervals, cases = [], [], []
    for job in matrix["jobs"]:
        path = out/"pharlap_raw"/(job["job_id"]+".mat")
        if not path.is_file():
            cases.append({"job_id":job["job_id"],"status":"raw_missing","diagnostic_coverage_passed":False})
            continue
        data = read_raw(path)
        raw = data["raw"]
        for key in ("job_sha256","grid_sha256","exporter_sha256","mex_sha256"):
            assert str(raw[key])==job[key], "Raw identity mismatch: "+key
        initial_errors = {key:float(vec(raw["ray_state_vec"],"geometrical_path" if key=="geom_path" else key)[0]-job["initial"][key]) for key in STATE_INPUT}
        native_s = data["s"]-data["s"][0]
        stop = min(int(np.searchsorted(native_s,10.,side="left")),len(native_s)-1)
        selected = range(stop+1)
        points.extend(point_row(data,index,job["job_id"],"restart_native_10km") for index in selected)
        intervals.extend(interval_row(data,index,job["job_id"],"restart_native_10km") for index in range(stop))
        exterior = [index for index in range(stop) if data["height"][index]<60. and data["height"][index+1]<60.]
        ratios = [interval_row(data,index,job["job_id"],"outside") for index in exterior]
        outside_p = data["p"][1:stop+1][data["height"][1:stop+1]<60.]
        unit = np.array([job["initial"][k] for k in STATE_INPUT[3:6]])
        unit /= np.linalg.norm(unit)
        cases.append({"job_id":job["job_id"],"status":"completed","level":job["level"],"position_role":job["position_role"],
                      "normalize_input_p":job["normalize_input_p"],"diagnostic_only":True,
                      "common_p_inside_is_not_a_physical_initial_state":job["common_p_inside_is_not_a_physical_initial_state"],
                      "raw_sha256":sha(path),"raw_point_count":len(native_s),"diagnostic_point_count":stop+1,
                      "covered_span_km":float(native_s[stop]),"diagnostic_coverage_passed":bool(native_s[-1]>=10.),
                      "initial_field_errors":initial_errors,"initial_p2":float(np.sum(data["p"][0]**2)),
                      "next_p2":float(np.sum(data["p"][1]**2)),
                      "initial_s_preserved":bool(initial_errors["geom_path"]==0.),"initial_G_preserved":bool(initial_errors["group_path"]==0.),
                      "first_interval":interval_row(data,0,job["job_id"],"first_progressed_interval"),
                      "outside_interval_count":len(ratios),
                      "outside_Gprime_minus_one_max":max([abs(v["Gprime_minus_one"]) for v in ratios],default=None),
                      "outside_p2_minus_one_after_initial_max":float(np.max(abs(np.sum(outside_p**2,axis=1)-1.))) if len(outside_p) else None,
                      "outside_p_minus_input_unit_max":float(np.max(np.linalg.norm(outside_p-unit,axis=1))) if len(outside_p) else None,
                      "state_path_nonterminal_max_m":float(np.max(np.linalg.norm(data["xyz"][:-1]-data["path_xyz"][:-1],axis=1))*1000),
                      "terminal_state_path_error_m":float(np.linalg.norm(data["xyz"][-1]-data["path_xyz"][-1])*1000)})
    for name,rows in (("natural_points",natural_points),("natural_intervals",natural_intervals),("restart_points",points),("restart_intervals",intervals)):
        if rows: write_csv(analysis/(name+".csv"),rows)
    write_csv(analysis/"natural_cases.csv",natural_cases)
    flat_cases = []
    for case in cases:
        flat = {k:v for k,v in case.items() if not isinstance(v,dict)}
        flat.update({"initial_error_"+k:v for k,v in case.get("initial_field_errors",{}).items()})
        flat.update({"first_interval_"+k:v for k,v in case.get("first_interval",{}).items() if isinstance(v,(int,float))})
        flat_cases.append(flat)
    if flat_cases:
        fields = list(dict.fromkeys(k for row in flat_cases for k in row))
        with (analysis/"restart_cases.csv").open("w",encoding="utf-8",newline="") as handle:
            writer = csv.DictWriter(handle,fieldnames=fields)
            writer.writeheader(); writer.writerows(flat_cases)
    original = loadmat(Path(matrix["jobs"][0]["grid_mat"]),simplify_cells=True,mat_dtype=True)
    back_path = out/"grid_readback.mat"
    back = loadmat(back_path,simplify_cells=True,mat_dtype=True) if back_path.is_file() else {}
    checks["grid_readback_arrays_equal"] = bool(back and all(k in back and np.array_equal(original[k],back[k]) for k in original if not k.startswith("__")))
    checks["all_12_controls_completed"] = bool(len(cases)==12 and all(c["status"]=="completed" for c in cases))
    checks["all_10km_spans_covered"] = all(c["diagnostic_coverage_passed"] for c in cases)
    summary = {"schema":"c1m-p2-exit-semantics-v1","checks":checks,"natural_cases":natural_cases,"restart_cases":cases,
               "sources_and_input_integrity_passed":bool(all(checks.values())),
               "physical_acceptance_passed":False,"global_passed":False,"closed_core_location_identified":False,
               "independent_ODE_modified_to_fit_PH":False,
               "interpretation_scope":"Public explicit-state input, native state/path output and increments only. Closed-core integration/save/reconstruction location is not identified.",
               "error_notes":"Final ground path projection is separate from nonterminal state/path comparisons. No truncated failed trace substitutes for 10 km coverage."}
    write_json(analysis/"summary.json",summary)
    verify(manifest,out)
    hashes = {str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file() and p.name!="sha256_manifest.json"}
    write_json(analysis/"sha256_manifest.json",{"records":hashes})
    print("P2 EXIT ANALYZED integrity="+str(summary["sources_and_input_integrity_passed"])+" controls="+str(len(cases)))


def run(run: Path) -> None:
    out = run.resolve()/"exit"
    verify(read_json(out/"manifest.json"),out)
    assert not (out/"execution.json").exists(), "Execution already exists; use a new run"
    command = [str(MATLAB),"-singleCompThread","-batch",matlab_expression(out)]
    write_json(out/"launch.json",{"command":command,"threads":1,"start_utc":datetime.now(timezone.utc).isoformat()})
    with (out/"matlab_stdout.log").open("w",encoding="utf-8") as log:
        process = subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy(),
                                 creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0),check=False)
    write_json(out/"execution.json",{"returncode":process.returncode,"end_utc":datetime.now(timezone.utc).isoformat()})
    verify(read_json(out/"manifest.json"),out)
    analyze(run)
    assert process.returncode==0, "MATLAB failed; execution and raw records preserved"


def selfcheck() -> None:
    assert len(STATE_INPUT)==12 and len(TOLS)==3 and all(t[2]>=1. for t in TOLS)
    data = {"s":np.array([100.,101.]),"G":np.array([200.,201.]),"path_s":np.array([100.,101.]),
            "path_G":np.array([200.,201.]),"xyz":np.array([[6430.,0.,0.],[6429.,0.,0.]]),
            "height":np.array([59.,58.]),"p":np.array([[-.9,0.,0.],[-1.,0.,0.]]),"indep_var":np.array([5.,6.])}
    row = interval_row(data,0,"synthetic","test")
    assert row["Gprime_state"]==1. and abs(row["Gprime_minus_left_inverse_p"]-(1.-1./.9))<1e-15
    assert row["direction_change"]==0.
    assert json.loads(json.dumps({"passed":True}))["passed"] is True
    print("P2 EXIT SELFCHECK 0 fail(s)")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run",type=Path)
    parser.add_argument("--stage",choices=("prepare","run","analyze","selfcheck"),required=True)
    parser.add_argument("--analysis-id",default="analysis_v1")
    args = parser.parse_args()
    if args.stage=="selfcheck": selfcheck()
    else:
        assert args.run is not None, "--run required"
        if args.stage=="prepare": prepare(args.run)
        elif args.stage=="run": run(args.run)
        else: analyze(args.run,args.analysis_id)
