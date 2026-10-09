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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("inventory", "before", "after"))
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "inventory":
        inventory(args.run.resolve())
    else:
        capture(args.run.resolve(), args.stage)
