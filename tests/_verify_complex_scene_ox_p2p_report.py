"""Verify five-scene O/X ray fans and seeded local P2P report products."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2"
ASSETS = OUT / "report_assets" / "reflection_route_2026-07-26"
SCENES = ("baseline", "tid", "spreadf", "es", "bubble")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path: Path):
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def run() -> int:
    frozen = _load(ASSETS / "formal_complex_scene_manifest.json")["inputs"]
    checks = {}
    metrics = {}
    artifacts = []

    for scene in SCENES:
        o_png = OUT / "tx_fan_unified_{}.png".format(scene)
        o_json = OUT / "ray_fan_unified_{}.json".format(scene)
        ox_png = OUT / "tx_fan_ox_unified_{}.png".format(scene)
        ox_json = OUT / "ray_fan_ox_unified_{}.json".format(scene)
        p2p_png = OUT / "p2p_modes_unified_{}.png".format(scene)
        p2p_json = OUT / "p2p_unified_{}.json".format(scene)
        p2p_csv = OUT / "modes_unified_{}.csv".format(scene)
        artifacts.extend((o_png, o_json, ox_png, ox_json, p2p_png,
                          p2p_json, p2p_csv))

        frozen_png = frozen[
            "output2/tx_fan_unified_{}.png".format(scene)]["sha256"]
        frozen_json = frozen[
            "output2/ray_fan_unified_{}.json".format(scene)]["sha256"]
        checks["{}_o_png_frozen".format(scene)] = _sha256(o_png) == frozen_png
        checks["{}_o_json_frozen".format(scene)] = (
            _sha256(o_json) == frozen_json)

        ox = _load(ox_json)
        p2p = _load(p2p_json)
        diagnostics = ox["ray_fan_diagnostics"]
        checks["{}_ox_modes".format(scene)] = set(ox["wave_modes"]) == {"O", "X"}
        for branch in ("O", "X"):
            diag = diagnostics[branch]
            accounted = int(diag["returned"]) + sum(
                int(value) for value in diag["exception_counts"].values())
            checks["{}_{}_attempted_33".format(scene, branch)] = (
                int(diag["attempted"]) == 33)
            checks["{}_{}_attempts_accounted".format(scene, branch)] = (
                accounted == 33)
        checks["{}_provider_hash_consistent".format(scene)] = (
            ox["provider_hash"] == p2p["provider_hash"]
            == _load(o_json)["provider_hash"])
        checks["{}_p2p_scope".format(scene)] = (
            p2p["search_scope"]
            == "seeded_local_roots_not_full_domain_mode_audit"
            and p2p["completeness_claim"] is False
            and int(p2p["candidate_nearest_per_branch"]) == 1
            and p2p["solve_profile"] == "smoke"
            and p2p["final_retrace_profile"] == "research"
            and int(p2p["max_nfev"]) == 5)
        checks["{}_images_valid".format(scene)] = (
            ox_png.stat().st_size > 10_000
            and p2p_png.stat().st_size > 10_000)

        modes = p2p["modes"]
        valid_modes = all(
            mode.get("termination_status") == "ground"
            and mode.get("p2p_trace_profile") == "research"
            and mode.get("p2p_optimization_trace_profile") == "smoke"
            and mode.get("p2p_search_scope")
            == "seeded_local_roots_not_full_domain_mode_audit"
            and np.isfinite(float(mode["p2p_landing_error_km"]))
            and float(mode["p2p_landing_error_km"]) < 5.0
            and np.isfinite(float(mode["apex_height_km"]))
            and np.isfinite(float(mode["group_path_km"]))
            and np.isfinite(float(mode["tau_ms"]))
            for mode in modes)
        checks["{}_p2p_roots_valid".format(scene)] = valid_modes
        with p2p_csv.open("r", encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        checks["{}_csv_matches_roots".format(scene)] = len(rows) == len(modes)
        metrics[scene] = {
            "provider_hash": ox["provider_hash"],
            "ray_fan": {
                branch: diagnostics[branch] for branch in ("O", "X")
            },
            "p2p_root_counts": p2p["root_counts"],
            "p2p_max_landing_error_km": (
                max(float(mode["p2p_landing_error_km"]) for mode in modes)
                if modes else None),
            "p2p_status": (
                "seeded_roots_found" if modes
                else "no_converged_local_seeded_root_not_no_mode_evidence"),
        }

    checks = {name: bool(value) for name, value in checks.items()}
    report = {
        "schema_version": "complex_scene_ox_p2p_report_acceptance_v1",
        "status": "passed" if all(checks.values()) else "failed",
        "search_scope": "seeded_local_roots_not_full_domain_mode_audit",
        "complex_scene_mode_completeness": "not_verified_p2_6h_deferred",
        "checks": checks,
        "metrics": metrics,
    }
    acceptance = ASSETS / "formal_complex_scene_ox_p2p_acceptance.json"
    acceptance.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    artifacts.append(acceptance)
    manifest = {
        "schema_version": "complex_scene_ox_p2p_report_manifest_v1",
        "artifacts": {
            str(path.relative_to(ROOT)).replace("\\", "/"): {
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in artifacts
        },
        "code": {
            str(path.relative_to(ROOT)).replace("\\", "/"): _sha256(path)
            for path in (
                ROOT / "main2.py",
                ROOT / "models" / "unified_hybrid_pipeline.py",
                ROOT / "models" / "point_to_point_ah_3d_ecef.py",
                ROOT / "viz" / "plot_utils.py",
                Path(__file__),
            )
        },
    }
    (ASSETS / "formal_complex_scene_ox_p2p_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    failures = [name for name, value in checks.items() if not value]
    print("complex_scene_ox_p2p_report: {} fail(s)".format(len(failures)))
    for name in failures:
        print(" FAIL {}".format(name))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(run())
