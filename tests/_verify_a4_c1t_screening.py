"""Verify A4-C1T screening artifacts without running external solvers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1t"
MANIFEST = OUT / "c1t_freeze_manifest.json"
ALLOWED_STATUS = {
    "qualified_hamiltonian_solver_found",
    "no_qualified_solver_use_farr_local_fullwave",
    "validation_blocked_by_build_or_compute_resources",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    failures = []
    if not MANIFEST.is_file():
        failures.append("missing C1T manifest")
    else:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        if manifest.get("status") not in ALLOWED_STATUS:
            failures.append("invalid final status")
        if not all(manifest.get("checks", {}).values()):
            failures.append("one or more frozen checks failed")
        for raw_path, expected in manifest.get("source_sha256", {}).items():
            path = Path(raw_path)
            if not path.is_file() or sha256(path) != expected:
                failures.append(f"source drift: {path}")
        for item in manifest.get("frozen_artifacts", []):
            path = Path(item["path"])
            if not path.is_file() or sha256(path) != item["sha256"]:
                failures.append(f"artifact drift: {path}")
            elif path.suffix.lower() == ".png":
                pixels = np.asarray(mpimg.imread(path), dtype=float)
                if (pixels.size == 0 or not np.isfinite(pixels).all()
                        or min(pixels.shape[:2]) < 200 or np.var(pixels) <= 1e-8):
                    failures.append(f"invalid visualization: {path}")

    matrix_path = OUT / "c1t_candidate_admission_matrix.json"
    farr_path = OUT / "c1t_farr_reproduction_config.json"
    status_path = OUT / "c1t_status.json"
    for path in (matrix_path, farr_path, status_path):
        if not path.is_file():
            failures.append(f"missing artifact: {path}")
    if not failures:
        matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
        candidates = matrix.get("candidates", [])
        if matrix.get("qualified_candidates"):
            failures.append("an unverified third solver was admitted")
        if matrix.get("screening_conclusion") != (
                "no_qualified_solver_use_farr_local_fullwave"):
            failures.append("unexpected screening conclusion")
        au_ray = next((row for row in candidates
                       if row["candidate"] == "AU-Ray public demo"), None)
        farr = next((row for row in candidates
                     if row["candidate"] == "FARR"), None)
        if au_ray is None or au_ray["full_solver_source"]:
            failures.append("AU-Ray source-audit semantics changed")
        if farr is None or farr["independent_hamiltonian"]:
            failures.append("FARR Hamiltonian semantics changed")
        config = json.loads(farr_path.read_text(encoding="utf-8"))
        if (config["paper_case"]["frequency_mhz"] != 1.0
                or config["paper_case"]["domain_km"] != [150.0, 37.5, 120.0]
                or config["paper_case"]["mpi_processes_reported"] != 3000
                or len(config.get("levels", [])) != 3
                or config["acceptance"]["ray_rms_thresholds_applicable"]):
            failures.append("FARR paper-case or semantic contract drifted")
        status = json.loads(status_path.read_text(encoding="utf-8"))
        if (status["formal_c1_status"]
                != "failed_c1_chapman_high_low_external_comparison"
                or status["c1q_attribution"]
                != "pharlap_effective_dynamics_not_recoverable"
                or status["production_modified"]
                or status["a4_c2_admitted"]):
            failures.append("C1/C1Q/manual-gate semantics changed")

    print(json.dumps({"failures": failures}, indent=2, ensure_ascii=True))
    print(f"A4-C1T verification: {len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
