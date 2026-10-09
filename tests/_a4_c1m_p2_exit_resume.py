"""Read-only completion lookup for the original batch plus identical-job retry."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil


def sha(path: Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle,"sha256").hexdigest()


def read(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def finalize(run: Path) -> dict:
    out = run.resolve()/"exit"
    index_path = out/"completion_index_v3.json"
    assert not index_path.exists(), "Completion index already frozen"
    source = Path(__file__).resolve()
    snapshots = out/"completion_sources_v3"
    snapshots.mkdir(exist_ok=False)
    snapshot = snapshots/source.name
    shutil.copy2(source,snapshot)
    write(snapshots/"source_freeze.json",{"source":str(source),"sha256":sha(source),"snapshot":str(snapshot)})
    summary_path = out/"retry_nomem_v1"/"analysis_v3"/"summary.json"
    summary = read(summary_path)
    combined = read(out/"retry_nomem_v1"/"analysis_v2"/"summary.json")
    assert summary["diagnostic_execution_complete"] and all(summary["checks"].values())
    assert all(combined["checks"].values()) and len(summary["raw_sources"])==12
    assert summary["original_batch_returncode"]==1 and summary["retry_returncode"]==0
    inputs = {}
    for path,groups in ((out/"manifest.json",("inputs_sha256","external_sha256")),
                        (out/"retry_nomem_v1"/"manifest.json",("input_sha256",)),
                        (out/"retry_nomem_v1"/"sequence_manifest_v3.json",("inputs_sha256","sources_sha256"))):
        manifest = read(path)
        for group in groups:
            for key,value in manifest[group].items():
                assert key not in inputs or inputs[key]==value
                inputs[key] = value
    for record in read(out/"manifest.json")["sources"]:
        inputs[record["path"]] = record["sha256"]
    inputs[str(source)] = sha(source)
    failed = run.resolve()/"exit_prepare_failed_v1"
    if failed.is_dir():
        inputs.update({str(p):sha(p) for p in failed.rglob("*") if p.is_file() and "__pycache__" not in p.parts})
    assert all(sha(Path(k))==v for k,v in inputs.items())
    artifacts = {str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    index = {"schema":"c1m-p2-exit-completion-index-v3","summary":str(summary_path.relative_to(out)),
             "inputs_sha256":inputs,"artifacts_sha256":artifacts,"resume_source_sha256":sha(source),
             "successful_job_count":12,"raw_sources":summary["raw_sources"],
             "original_attempt_returncode":1,"identical_job_retry_returncode":0,
             "component_diagnostic_execution_complete":True,"physical_acceptance_passed":False,
             "global_passed":False,"closed_core_location_identified":False}
    write(index_path,index)
    return verify_completion(run)


def verify_completion(run: Path) -> dict:
    out = run.resolve()/"exit"
    index = read(out/"completion_index_v3.json")
    assert sha(Path(__file__))==index["resume_source_sha256"], "Completion-reader source changed"
    assert all(sha(Path(k))==v for k,v in index["inputs_sha256"].items()), "Component input/source changed"
    assert all(sha(out/k)==v for k,v in index["artifacts_sha256"].items()), "Component artifact changed"
    assert index["successful_job_count"]==12 and index["component_diagnostic_execution_complete"] is True
    for record in index["raw_sources"]:
        assert sha(Path(record["raw_mat"]))==record["sha256"]
    summary = read(out/index["summary"])
    print("P2 EXIT READ-ONLY VERIFIED 12 controls; frozen summary="+str(out/index["summary"]))
    return summary


def run(run: Path) -> dict | None:
    out = run.resolve()/"exit"
    if (out/"completion_index_v3.json").is_file():
        return verify_completion(run)
    assert not (out/"execution.json").exists(), "Partial run preserved; explicit versioned recovery required"
    source = out/"sources"/"_a4_c1m_p2_exit.py"
    expected = next(v["sha256"] for v in read(out/"manifest.json")["sources"] if Path(v["snapshot"]).name==source.name)
    assert sha(source)==expected
    spec = importlib.util.spec_from_file_location("frozen_exit",source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run(run)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--stage",choices=("finalize","run"),required=True)
    args = parser.parse_args()
    {"finalize":finalize,"run":run}[args.stage](args.run)
