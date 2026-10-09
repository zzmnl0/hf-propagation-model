"""Exact lightweight production suite; historical tests are never discovered."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
TESTS = (
    "_verify_production_boundary.py",
    "_verify_a4_c0r1_geometry.py",
    "_verify_a4_c0r3_ecef_hamiltonian.py",
    "_verify_absorption_consistency.py",
    "_verify_a4_c0r6_entry_migration.py",
    "_verify_legacy_source_separation.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--test", choices=TESTS, action="append")
    parser.add_argument("--child", choices=TESTS, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.list:
        print("\n".join(TESTS))
        return 0
    if args.output is None:
        parser.error("--output must name a new run directory")
    out = args.output.resolve()
    sys.path.insert(0, str(ROOT))
    from _verify_production_boundary import install_guard, check_loaded, static_violations
    if args.child:
        # NumPy consults cached platform data; Windows' initial query invokes cmd /c ver.
        platform.uname()
        iri_driver = None
        if args.child == "_verify_a4_c0r6_entry_migration.py":
            iri_spec = importlib.util.find_spec("iri2016")
            iri_driver = Path(iri_spec.origin).parent / "iri2016_driver.exe"
        allowed_calls = install_guard(iri_driver)
        spec = importlib.util.spec_from_file_location("production_case", ROOT / "tests" / args.child)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # Existing acceptance scripts retain their assertions but write only into this run.
        for key in ("OUT", "REPORT", "FIGURE", "CSV_DIR"):
            if hasattr(module, key):
                old = getattr(module, key)
                setattr(module, key, out if key == "OUT" else out / old.name)
        try:
            result = (module.run if hasattr(module, "run") else module.main)()
        finally:
            (out / "subprocess_audit.json").write_text(json.dumps(allowed_calls, indent=2) + "\n", encoding="utf-8")
            check_loaded()
        return result
    out.mkdir(parents=True, exist_ok=False)
    selected = args.test or list(TESTS)
    paths = list(ROOT.glob("*.py"))
    for folder in ("models", "viz", "analysis"):
        paths.extend((ROOT / folder).rglob("*.py"))
    paths.extend(ROOT / "tests" / test for test in selected if test not in (
        "_verify_production_boundary.py", "_verify_legacy_source_separation.py"))
    violations = static_violations(paths)
    rows = []
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", MPLBACKEND="Agg")
    if not violations:
        for test in selected:
            case_out = out / Path(test).stem
            case_out.mkdir()
            process = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                      "--child", test, "--output", str(case_out)],
                                     cwd=ROOT, env=env, capture_output=True, text=True,
                                     encoding="utf-8", errors="replace")
            (case_out / "run.log").write_text(process.stdout + process.stderr, encoding="utf-8")
            rows.append({"test": test, "returncode": process.returncode,
                         "passed": process.returncode == 0})
            print(f"{test}: {'PASS' if process.returncode == 0 else 'FAIL'}", flush=True)
    report = {"scope": "exact lightweight suite; only pinned IRI native driver subprocess allowed for C0R6",
              "guard_scope": "project code/legacy imports, reads and execution; independent historical hash verification runs outside this suite",
              "static_fixture_exclusion": "boundary rejection fixtures and legacy_source_separation directory inventory are excluded from literal scanning; both execute under runtime guard",
              "selected_tests": selected, "excluded_tests": sorted(path.name for path in (ROOT / "tests").glob("_verify*.py") if path.name not in selected),
              "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths + [Path(__file__)] + [ROOT / "tests" / test for test in selected]},
              "static_violations": violations, "cases": rows,
              "passed": not violations and len(rows) == len(selected) and all(row["passed"] for row in rows)}
    report["full_suite_passed"] = report["passed"] and set(selected) == set(TESTS)
    report["full_suite_tests"] = list(TESTS)
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
