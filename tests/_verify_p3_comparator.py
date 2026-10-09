"""Fault injection into the real P3 comparator, preserving source raw states."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import _a4_c1m_p3_compare as compare


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    data = run / "experiments"
    output = args.output.resolve() if args.output else run / "comparator_negative_controls"
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    shutil.copy2(source, output / source.name)
    manifest = json.loads((data / "manifest.json").read_text(encoding="utf-8"))
    jobs = [job for job in manifest["jobs"]
            if job["role"] == "R1" and job["job_id"].startswith("internal_identified_f10_e22_t")]
    assert len(jobs) == 3
    jobs.sort(key=lambda job: job["level"])
    original_read = compare.read_trace
    original_reference = compare.p2.RadialReference
    originals = {(side, job["level"]): original_read(data, job, side)
                 for side in ("ah", "ph") for job in jobs}
    source_hashes = {trace["raw_path"]: sha(trace["raw_path"]) for trace in originals.values()}
    reference_cache = {}

    def cached_reference(*pos, **kw):
        # All faults affect AH only; the identical PH reference is reused exactly.
        key = (pos[5], pos[6])
        if key not in reference_cache:
            reference_cache[key] = original_reference(*pos, **kw)
        return reference_cache[key]

    def evaluate(name, injection=None):
        directory = output / name
        directory.mkdir()
        def reader(_data, job, side):
            trace = copy.deepcopy(originals[(side, job["level"])])
            if injection and side == "ah":
                injection(trace, job["level"])
            return trace
        compare.read_trace = reader
        row = compare.compare_case(data, jobs, directory)
        (directory / "result.json").write_text(json.dumps(row, indent=2, allow_nan=False), encoding="utf-8")
        return row

    results = []
    try:
        compare.p2.RadialReference = cached_reference
        baseline = evaluate("baseline")
        assert baseline["checks"]["thresholds"] and baseline["checks"]["initial_state"]
        assert baseline["checks"]["scope"]

        def wrong_units(trace, level):
            trace["xyz"] *= 1000.0
        units = evaluate("wrong_position_units", wrong_units)
        assert not units["checks"]["thresholds"] and not units["strict_passed"]
        results.append({"name": "wrong_position_units", "fault": "AH km positions interpreted as metres",
                        "rejected_gate": "thresholds", "passed": True})

        def wrong_initial(trace, level):
            if level == 1:
                trace["p"][0] *= 1.001
        initial = evaluate("coarse_level_initial_mismatch", wrong_initial)
        assert not initial["checks"]["initial_state"] and not initial["strict_passed"]
        results.append({"name": "coarse_level_initial_mismatch", "fault": "only AH level 1 initial p scaled by 1.001",
                        "rejected_gate": "initial_state", "passed": True})

        def truncate(trace, level):
            if level == 2:
                keep = trace["s"] < 90.0
                for key in ("s", "xyz", "p", "group", "tangent"):
                    trace[key] = trace[key][keep]
                trace["termination"] = "grid_exit"
        try:
            evaluate("premature_termination", truncate)
        except ValueError as error:
            assert "Early termination" in str(error)
            results.append({"name": "premature_termination", "fault": "AH level 2 truncated before 90 km",
                            "rejected_gate": "scope", "error": str(error), "passed": True})
        else:
            raise AssertionError("premature termination was accepted")

        def nonmonotonic(trace, level):
            displacement = {1: 0.0, 2: 0.001, 3: 0.1}[level]
            trace["xyz"][:, 0] += displacement * trace["s"] / 100.0
        convergence = evaluate("nonmonotonic_refinement", nonmonotonic)
        control = convergence["integration"]["ah"]
        assert control["delta23"]["rms_m"] > control["delta12"]["rms_m"]
        assert not convergence["checks"]["integration"] and not convergence["strict_passed"]
        results.append({"name": "nonmonotonic_refinement", "fault": "AH terminal x offsets 0,1,100 m over three levels",
                        "rejected_gate": "integration", "passed": True})
    finally:
        compare.read_trace = original_read
        compare.p2.RadialReference = original_reference
        assert all(sha(path) == digest for path, digest in source_hashes.items())
    summary = {"passed": len(results) == 4 and all(row["passed"] for row in results),
               "tests": results, "source_sha256": sha(source),
               "comparator_sha256": sha(compare.__file__),
               "reference_source_sha256": sha(compare.p2.__file__),
               "fixture_sha256": source_hashes, "raw_inputs_unchanged": True,
               "reference_cache": "real PH references reused; no reference or acceptance formulas mocked"}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"passed": summary["passed"], "checks": len(results)}), flush=True)


if __name__ == "__main__":
    main()
