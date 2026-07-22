"""Freeze pre-C0R evidence and inventory legacy-coordinate production callers."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "output2" / "a4_pharlap" / "c0r0_freeze"
MANIFEST = FREEZE / "c0r0_freeze_manifest.json"
LEGACY = "legacy_coordinate_reference_only"

C0_FILES = (
    "output2/a4_pharlap/a4_c0_validation.json",
    "output2/a4_pharlap/run.log",
    "output2/a4_pharlap/figures/trajectory_c0_ah3d.png",
    "output2/a4_pharlap/figures/trajectory_c0_comparison.png",
    "output2/a4_pharlap/figures/trajectory_c0_pharlap.png",
    "output2/a4_pharlap/inputs/a4_c0_vacuum.mat",
    "output2/a4_pharlap/preflight/pharlap_preflight.json",
    "output2/a4_pharlap/preflight/pharlap_preflight.mat",
    "output2/a4_pharlap/raw/a4_c0_vacuum_ah3d.json",
    "output2/a4_pharlap/raw/a4_c0_vacuum_pharlap.mat",
    "tests/_run_a4_pharlap_c0.py",
    "tests/_verify_a4_pharlap_c0.py",
    "word/a4_c0_implementation_report.txt",
)

SCAN_TOKENS = (
    "trajectory_3d",
    "earth_flat",
    "TX_POS",
    "RX_POS",
    "z_cart",
    "legacy_parabolic",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_once(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if _sha256(source) != _sha256(destination):
            raise RuntimeError(f"frozen artifact differs: {destination}")
        return
    shutil.copy2(source, destination)


def _record(source: Path, destination: Path, artifact_class: str) -> dict:
    _copy_once(source, destination)
    return {
        "source_path": str(source),
        "frozen_path": str(destination),
        "sha256": _sha256(source),
        "bytes": source.stat().st_size,
        "artifact_class": artifact_class,
        "coordinate_semantics": LEGACY,
    }


def _legacy_sources() -> list[tuple[Path, str]]:
    items: list[tuple[Path, str]] = []
    groups = (
        (ROOT / "output2" / "a3rnew_p3", "a3rnew_p3_internal_evidence"),
        (ROOT / "output2" / "a3rnew_p5", "a3rnew_p5_internal_evidence"),
    )
    for directory, artifact_class in groups:
        items.extend((path, artifact_class) for path in directory.rglob("*") if path.is_file())
    benchmark = ROOT / "output2" / "a3r_benchmark"
    for pattern in ("*p0c*", "*high_low*"):
        items.extend((path, "coleman_p0c_root_evidence")
                     for path in benchmark.glob(pattern) if path.is_file())
    for directory in (ROOT / "output", ROOT / "output2"):
        items.extend((path, "main_main2_figure") for path in directory.glob("*.png"))
    for name, artifact_class in (
        ("a3rnew_p3_implementation_report.txt", "a3rnew_p3_internal_report"),
        ("a3rnew_p5_accuracy_acceptance_report.txt", "a3rnew_p5_internal_report"),
    ):
        path = ROOT / "word" / name
        if path.is_file():
            items.append((path, artifact_class))
    return sorted(set(items), key=lambda item: item[0].as_posix().lower())


def _production_files() -> list[Path]:
    files = [ROOT / name for name in ("config.py", "main.py", "main2.py")]
    for directory in (ROOT / "models", ROOT / "viz", ROOT / "analysis"):
        files.extend(directory.rglob("*.py"))
    return sorted(path for path in set(files) if path.is_file())


def _migration_inventory() -> dict:
    occurrences = []
    for path in _production_files():
        relative = path.relative_to(ROOT).as_posix()
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for token in SCAN_TOKENS:
                if token in line:
                    occurrences.append({
                        "token": token,
                        "path": relative,
                        "line": line_number,
                        "text": line.strip(),
                    })
    callers = sorted({item["path"] for item in occurrences})
    return {
        "scan_scope": ["config.py", "main.py", "main2.py", "models/**/*.py",
                       "viz/**/*.py", "analysis/**/*.py"],
        "tokens": list(SCAN_TOKENS),
        "occurrence_count": len(occurrences),
        "caller_count": len(callers),
        "callers": callers,
        "occurrences": occurrences,
        "migration_rule": "all listed production callers migrate before entrypoint switch",
        "partial_migration_allowed": False,
    }


def _git_head() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
        capture_output=True, check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def run() -> int:
    missing = [relative for relative in C0_FILES if not (ROOT / relative).is_file()]
    matlab_script = Path(r"D:\Matlab\PHaRLAP\a4_pharlap_preflight.m")
    if not matlab_script.is_file():
        missing.append(str(matlab_script))
    if missing:
        print("A4-C0R-0: FAIL missing=" + ",".join(missing))
        return 1

    if MANIFEST.exists():
        print(f"A4-C0R-0: frozen manifest already exists: {MANIFEST}")
        return 0

    frozen = []
    for relative in C0_FILES:
        source = ROOT / relative
        destination = FREEZE / "c0_evidence" / relative
        frozen.append(_record(source, destination, "a4_c0_historical_evidence"))
    frozen.append(_record(
        matlab_script,
        FREEZE / "c0_evidence" / "external_matlab" / matlab_script.name,
        "a4_c0_matlab_adapter",
    ))

    legacy = []
    for source, artifact_class in _legacy_sources():
        relative = source.relative_to(ROOT)
        destination = FREEZE / LEGACY / relative
        legacy.append(_record(source, destination, artifact_class))

    manifest = {
        "schema_version": "a4-c0r0-freeze-v1",
        "stage": "A4-C0R-0",
        "status": "passed_legacy_coordinate_artifacts_frozen",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_head_at_freeze": _git_head(),
        "coordinate_semantics": LEGACY,
        "claim_boundary": (
            "historical internal/interface evidence only; not PHaRLAP external "
            "physics validation and not valid for the new ECEF production contract"
        ),
        "source_overwrite_allowed": False,
        "frozen_c0_artifacts": frozen,
        "legacy_coordinate_artifacts": legacy,
        "migration_inventory": _migration_inventory(),
        "required_new_trajectory_fields": [
            "trajectory_ecef_km",
            "trajectory_geocentric_lat_lon_height",
        ],
        "forbidden_entrypoint_switch_until_inventory_closed": True,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    temporary = MANIFEST.with_suffix(".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                         encoding="utf-8")
    temporary.replace(MANIFEST)
    print(
        "A4-C0R-0: PASS "
        f"c0={len(frozen)} legacy={len(legacy)} "
        f"callers={manifest['migration_inventory']['caller_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
