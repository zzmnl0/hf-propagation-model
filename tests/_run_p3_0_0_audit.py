"""Record workspace identity and exact production physics before/after cleanup."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True, encoding="utf-8")


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def capture(run, stage):
    import numpy as np
    from models.ne_provider_geocentric import ChapmanGeocentricNeProvider3D, VacuumGeocentricNeProvider3D
    from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF, integrate_fixed_path_ah3d_ecef
    import main
    import main2
    import config

    out = run / stage
    out.mkdir(parents=True, exist_ok=False)
    arrays, summaries = {}, []
    for label, provider, branch, gyro in (
        ("vacuum", VacuumGeocentricNeProvider3D(6371.), "O", 0.),
        ("chapman_b0", ChapmanGeocentricNeProvider3D(6371.), "O", 0.),
        ("chapman_O", ChapmanGeocentricNeProvider3D(6371.), "O", 1.197),
        ("chapman_X", ChapmanGeocentricNeProvider3D(6371.), "X", 1.197),
    ):
        model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371., freq_MHz=10.,
                                       wave_mode=branch, gyro_MHz=gyro,
                                       b_ecef_hat=(.2, -.3, .9327379053))
        result = integrate_fixed_path_ah3d_ecef(model, tx_lat_deg=30., tx_lon_deg=120.,
            tx_height_km=100., elevation_deg=35., bearing_deg=37., path_length_km=10., profile="smoke")
        assert result["termination_status"] == "fixed_geometric_path_completed", result["termination_status"]
        for key in ("trajectory_ecef_km", "wavevector_ecef_km^-1", "trajectory_geom_path_km",
                    "group_path_km", "geom_path_km", "tau_ms"):
            arrays[label + ":" + key] = np.asarray(result[key])
        summaries.append({"case": label, "termination": result["termination_status"],
                          "accepted_steps": result["accepted_steps"]})
    import copy
    original = copy.deepcopy(config.UPDATE2)
    for entry in (main, main2):
        try:
            with entry._temporary_update2({"OX_MODE": "O"}):
                with entry._temporary_update2({"OX_MODE": "X"}):
                    raise RuntimeError("fixture")
        except RuntimeError as exc:
            assert str(exc) == "fixture"
        assert config.UPDATE2 == original
    assert not any(n == "legacy" or n.startswith("legacy.") for n in sys.modules)
    np.savez(out / "physics.npz", **arrays)
    write(out / "physics.json", {"cases": summaries, "configuration_restored": True,
                                 "legacy_loaded": False})
    if stage == "after":
        previous = np.load(run / "before/physics.npz", allow_pickle=False)
        checks = {k: bool(np.array_equal(v, previous[k])) for k, v in arrays.items()}
        checks["same_fields"] = set(previous.files) == set(arrays)
        checks["same_statuses"] = json.loads((run / "before/physics.json").read_text()) == json.loads((out / "physics.json").read_text())
        write(out / "equivalence.json", {"checks": checks, "passed": all(checks.values())})
        assert all(checks.values()), checks
        # Compare only the old public power helper, without importing a historical solver.
        source = git("show", "54fb763:main.py").lstrip("\ufeff")
        function = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "_add_power_3d")
        namespace = {"cfg": config, "np": np,
                     "_tag_communication_modes": main._tag_communication_modes,
                     "_tag_point_target_radar_modes": main._tag_point_target_radar_modes}
        exec(compile(ast.Module(body=[function], type_ignores=[]), "git-baseline-main-power", "exec"), namespace)
        power_checks = {}
        for enabled in (False, True):
            for radar in (False, True):
                previous_mode, current_mode = copy.deepcopy(result), copy.deepcopy(result)
                with main._temporary_update2({"ABSORPTION_COMPLEX": {"enable": enabled}}):
                    namespace["_add_power_3d"]([previous_mode], model.provider, radar_mode=radar)
                    main._add_power_3d([current_mode], model.provider, radar_mode=radar)
                added = set(previous_mode) - set(result)
                power_checks[f"complex{enabled}_radar{radar}"] = all(previous_mode[k] == current_mode[k] for k in added)
        write(out / "power_equivalence.json", {"checks": power_checks, "passed": all(power_checks.values())})
        assert all(power_checks.values()), power_checks
    print(stage + " physics captured; cases=" + str(len(summaries)), flush=True)


def inventory(run):
    run.mkdir(parents=True, exist_ok=False)
    files = git("ls-files", "--cached", "--others", "--exclude-standard", "-z").split("\0")
    records = []
    for rel in sorted(set(filter(None, files))):
        path = ROOT / rel
        if not path.is_file():
            records.append({"path": rel, "exists": False})
            continue
        deps = []
        if path.suffix == ".py":
            tree = ast.parse(path.read_text(encoding="utf-8-sig"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    deps.extend(a.name for a in node.names)
                elif isinstance(node, ast.ImportFrom):
                    deps.append("." * node.level + (node.module or ""))
        records.append({"path": rel, "sha256": sha(path), "bytes": path.stat().st_size,
                        "imports": sorted(set(deps)), "tracked": bool(git("ls-files", "--", rel).strip())})
    write(run / "inventory_before.json", {"root": str(ROOT), "head": git("rev-parse", "HEAD").strip(),
        "branch": git("branch", "--show-current").strip(), "status": git("status", "--porcelain=v1"),
        "staged": git("diff", "--cached", "--name-only"), "records": records,
        "ignored": git("ls-files", "--others", "--ignored", "--exclude-standard", "--directory")})
    for phase in ("p1_20261007_v1", "p2_20261008_v1"):
        folder = ROOT / "output2/a4_pharlap/c1/c1m" / phase
        manifest = folder / "delivery_sha256.json"
        data = json.loads(manifest.read_text(encoding="utf-8"))
        rows = data.get("records", data.get("files", []))
        assert rows, manifest
        for row in rows:
            path = folder / row["path"]
            assert sha(path) == row["sha256"], path
        write(run / (phase + "_identity.json"), {"manifest": str(manifest), "sha256": sha(manifest), "records": len(rows), "verified": True})
    print("inventory and historical identities recorded", flush=True)


def finalize(run, suite):
    baseline = json.loads((run / "inventory_before.json").read_text(encoding="utf-8"))
    tests = json.loads((suite / "summary.json").read_text(encoding="utf-8"))
    assert tests["passed"] and tests["full_suite_passed"] and not tests["static_violations"]
    expected = {"_verify_production_boundary.py", "_verify_a4_c0r1_geometry.py",
                "_verify_a4_c0r3_ecef_hamiltonian.py", "_verify_absorption_consistency.py",
                "_verify_a4_c0r6_entry_migration.py", "_verify_legacy_source_separation.py"}
    assert set(tests["selected_tests"]) == expected and len(tests["cases"]) == len(expected)
    for rel, digest in tests["source_sha256"].items():
        assert sha(ROOT / rel) == digest, rel
    for name in ("equivalence.json", "power_equivalence.json"):
        assert json.loads((run / "after" / name).read_text())["passed"]
    unchanged = []
    for row in baseline["records"]:
        rel = row["path"]
        if rel.startswith("legacy/") or rel in (
            "models/dispersion_ah.py", "models/ray_tracer_ah_3d_ecef.py", "models/ne_provider_geocentric.py"):
            assert sha(ROOT / rel) == row["sha256"], rel
            unchanged.append(rel)
    histories = []
    for identity in run.glob("*_identity.json"):
        record = json.loads(identity.read_text())
        manifest = Path(record["manifest"])
        assert sha(manifest) == record["sha256"]
        data = json.loads(manifest.read_text(encoding="utf-8"))
        for row in data.get("records", data.get("files", [])):
            assert sha(manifest.parent / row["path"]) == row["sha256"]
        histories.append(record)
    owned = ["README.md", "config.py", "main.py", "main2.py", "models/unified_hybrid_pipeline.py",
             "tests/_run_p3_0_0_audit.py", "tests/_run_production_checks.py", "tests/_verify_production_boundary.py",
             "word/a4_c1m_p3_0_0_cleanup_report_20261009.txt"]
    owned += [r["path"] for r in baseline["records"] if not r.get("tracked", True)
              and Path(r["path"]).suffix in (".py", ".m", ".txt")]
    owned = sorted(set(owned))
    assert not git("status", "--porcelain", "--", *owned).strip(), "Task-owned changes must be committed"
    for rel in owned:
        assert git("ls-files", "--", rel).strip(), rel
    current = []
    original_records = {r["path"]: r for r in baseline["records"]}
    for rel in filter(None, git("ls-files", "--cached", "--others", "--exclude-standard", "-z").split("\0")):
        p = ROOT / rel
        if not p.is_file():
            continue
        category = ("historical_reference" if rel.startswith("legacy/") else
                    "active_production_or_research" if rel.startswith("models/") else
                    "formal_test" if Path(rel).name in tests["selected_tests"] else
                    "unselected_diagnostic_or_historical_tool" if rel.startswith("tests/") else
                    "documentation" if rel.startswith("word/") else "support_or_entrypoint")
        imports = []
        if p.suffix == ".py":
            for node in ast.walk(ast.parse(p.read_text(encoding="utf-8-sig"))):
                if isinstance(node, ast.Import):
                    imports.extend(a.name for a in node.names)
                elif isinstance(node, ast.ImportFrom):
                    imports.append("." * node.level + (node.module or ""))
        current.append({"path": rel, "category": category, "sha256": sha(p),
                        "tracked": bool(git("ls-files", "--", rel).strip()), "imports": sorted(set(imports)),
                        "changed_from_inventory": sha(p) != original_records.get(rel, {}).get("sha256")})
    write(run / "workspace_final.json", {"head": git("rev-parse", "HEAD").strip(),
        "branch": git("branch", "--show-current").strip(), "status": git("status", "--porcelain"),
        "owned_paths": owned, "records": current,
        "remaining_changes_policy": "Pre-existing other-task documentation edits and DOCX files retained, not committed by cleanup",
        "ignored_policy": "Existing output/output2/cache/work ignore rules retained; generated evidence linked from tracked report"})
    completion = {"phase": "P3-0-0", "solver_physics_modified": False, "production_code_modified": True, "backup_created": False,
        "inventory_complete": True, "production_legacy_isolated": True, "formal_tests_legacy_free": True,
        "legacy_guard_verified": True, "old_config_cleanup_complete": True, "active_entrypoints_preserved": True,
        "default_behavior_unchanged": True, "historical_evidence_unchanged": True,
        "git_tracking_complete": True, "cleanup_tests_passed": True, "p3_0_0_passed": True,
        "ready_for_p3_0": True, "p3_physics_accepted": False, "recommend_p4": False,
        "head": git("rev-parse", "HEAD").strip(), "baseline_head": baseline["head"],
        "suite": str(suite), "suite_summary_sha256": sha(suite / "summary.json"),
        "unchanged_legacy_and_physics_files": unchanged, "historical_manifests": histories,
        "scope": "exact selected formal suite plus four short forward rays and power-helper parity; not all full scenarios",
        "retained": ["legacy references", "archived config rejection keys", "active B1", "shared helpers with callers"],
        "excluded": "No historical legacy numerical tests, no P3 physical compatibility patch, no remote push"}
    write(run / "completion.json", completion)
    write(run / "delivery_sha256.json", {"records": [
        {"path": str(p.relative_to(run)), "sha256": sha(p), "bytes": p.stat().st_size}
        for p in sorted(run.rglob("*")) if p.is_file() and p.name != "delivery_sha256.json"]})
    print("P3-0-0 complete; task Git clean; physics unchanged; ready for P3-0", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("inventory", "before", "after", "finalize"))
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--suite", type=Path)
    args = parser.parse_args()
    if args.stage == "inventory":
        inventory(args.run.resolve())
    elif args.stage == "finalize":
        if args.suite is None:
            parser.error("finalize requires --suite")
        finalize(args.run.resolve(), args.suite.resolve())
    else:
        capture(args.run.resolve(), args.stage)
