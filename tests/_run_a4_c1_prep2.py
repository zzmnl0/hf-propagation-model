"""Run and freeze the PHaRLAP full-elevation Chapman Prep-2 scan."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MATLAB_DIR = Path(r"D:\Matlab\PHaRLAP")
MATLAB_SCRIPT = MATLAB_DIR / "a4_c1_prep2_pharlap.m"
MEX = PHARLAP / "mex" / "raytrace_3d_sp.mexw64"
PREP1 = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep1"
PREP1_REPORT = PREP1 / "a4_c1_prep1_validation.json"
PREP1_MANIFEST = PREP1 / "a4_c1_prep1_freeze_manifest.json"
OUT = ROOT / "output2" / "a4_pharlap" / "c1_prep" / "prep2"
RAW = OUT / "a4_c1_prep2_pharlap_raw.mat"
MATLAB_JSON = OUT / "a4_c1_prep2_pharlap_summary.json"
LOG = OUT / "a4_c1_prep2_matlab.log"
CSV = OUT / "a4_c1_prep2_scan.csv"
ROOT_CSV = OUT / "a4_c1_prep2_roots.csv"
REPORT = OUT / "a4_c1_prep2_validation.json"
MANIFEST = OUT / "a4_c1_prep2_freeze_manifest.json"
FIGURE = OUT / "a4_c1_prep2_full_domain_scan.png"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def _verify_prep1() -> tuple[Path, str]:
    manifest = json.loads(PREP1_MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "passed_input_contract_frozen":
        raise RuntimeError("Prep-1 is not frozen")
    for item in manifest["frozen_artifacts"]:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"Prep-1 frozen artifact changed: {path}")
    report = json.loads(PREP1_REPORT.read_text(encoding="utf-8"))
    main = next(item for item in report["inputs"] if item["variant"] == "main_2km")
    path = Path(main["path"])
    if _sha256(path) != main["sha256"]:
        raise RuntimeError("Prep-1 main input hash mismatch")
    return path, main["sha256"]


def _run_matlab(input_path: Path) -> subprocess.CompletedProcess:
    quote = lambda path: str(path).replace("'", "''")
    expression = (
        f"addpath('{quote(MATLAB_DIR)}'); "
        f"a4_c1_prep2_pharlap('{quote(input_path)}','{quote(RAW)}',"
        f"'{quote(MATLAB_JSON)}')")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run([str(MATLAB), "-batch", expression], env=env,
                            text=True, capture_output=True, timeout=7200,
                            errors="replace")
    LOG.write_text((result.stdout or "") + "\n" + (result.stderr or ""),
                   encoding="utf-8")
    return result


def _write_csv(path: Path, records: list[dict]) -> None:
    fields = [
        "scan_stage", "elevation_deg", "ray_label", "status",
        "ground_range_km", "along_track_residual_km", "cross_track_km",
        "apex_km", "geometric_path_km", "group_path_km", "tau_ms",
        "landing_lat_deg", "landing_lon_deg", "path_point_count",
        "bracket_low_deg", "bracket_high_deg", "exitflag", "root_solver",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def _plot_scan(summary: dict, source_hashes: dict[str, str]) -> dict:
    coarse = summary.get("coarse_records", [])
    adaptive = summary.get("adaptive_records", [])
    roots = summary.get("roots", [])
    colors = {"ground": "#009E73", "penetrated": "#E69F00",
              "grid_exit": "#D55E00", "max_points": "#CC79A7",
              "catastrophic_error": "#000000", "unknown": "#777777"}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), constrained_layout=True)
    for status, color in colors.items():
        selected = [record for record in coarse if record.get("status") == status]
        if selected:
            elevation = [record["elevation_deg"] for record in selected]
            axes[0].scatter(elevation, [record["ground_range_km"] for record in selected],
                            s=15, color=color, label=status, zorder=2)
            axes[1].scatter(elevation, [record["apex_km"] for record in selected],
                            s=15, color=color, label=status, zorder=2)
    if adaptive:
        axes[0].scatter([record["elevation_deg"] for record in adaptive],
                        [record["ground_range_km"] for record in adaptive],
                        s=7, marker="x", color="#0072B2", alpha=0.55,
                        label="adaptive 0.1 deg", zorder=1)
        axes[1].scatter([record["elevation_deg"] for record in adaptive],
                        [record["apex_km"] for record in adaptive],
                        s=7, marker="x", color="#0072B2", alpha=0.55,
                        label="adaptive 0.1 deg", zorder=1)
    axes[0].axhline(1185.848260273658, color="black", linestyle="--", linewidth=1.0,
                    label="RX target range")
    for index, root in enumerate(roots, start=1):
        axes[0].scatter(root["elevation_deg"], root["ground_range_km"], s=90,
                        marker="*", color="#D55E00", edgecolor="black", zorder=4,
                        label="bounded ground roots" if index == 1 else None)
        axes[1].scatter(root["elevation_deg"], root["apex_km"], s=90,
                        marker="*", color="#D55E00", edgecolor="black", zorder=4,
                        label="bounded ground roots" if index == 1 else None)
    axes[0].set(xlabel="Launch elevation (deg)", ylabel="Terminal ground range (km)",
                xlim=(1.0, 89.0))
    axes[1].set(xlabel="Launch elevation (deg)", ylabel="Apex height (km)",
                xlim=(1.0, 89.0))
    for label, axis in zip(("A", "B"), axes):
        axis.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=7, loc="best")
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top",
                  fontweight="bold")
    fig.suptitle(
        "A4-C1 Prep-2 | PHaRLAP full elevation forward-map scan\n"
        "one-way | Chapman F2 | 10 MHz | OX=0 | B=0 | collision=0 | "
        "not AH3D comparison",
        fontsize=9.5)
    fig.savefig(FIGURE, dpi=300, facecolor="white", metadata={
        "Title": "A4-C1 Prep-2 PHaRLAP full-domain scan",
        "Description": "PHaRLAP-only forward-map validation; not A4-C1 external validation",
    })
    plt.close(fig)
    image = mpimg.imread(FIGURE)
    checks = {
        "exists": FIGURE.is_file() and FIGURE.stat().st_size > 0,
        "finite_pixels": bool(np.all(np.isfinite(image))),
        "minimum_dimensions": bool(image.shape[0] >= 900 and image.shape[1] >= 1800),
        "nonblank_pixel_variance": bool(float(np.std(image[..., :3])) >= 1.0e-3),
        "source_hashes_recorded": all(len(value) == 64 for value in source_hashes.values()),
    }
    return {"path": str(FIGURE), "sha256": _sha256(FIGURE),
            "source_sha256": source_hashes, "pixel_shape": list(image.shape),
            "pixel_std": float(np.std(image[..., :3])), "checks": checks}


def _verify_existing(allow_rebuild: bool) -> int | None:
    if not MANIFEST.is_file():
        return None
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    ok = (manifest.get("status") == "passed_pharlap_full_domain_scan"
          and all(Path(item["path"]).is_file()
                  and _sha256(Path(item["path"])) == item["sha256"]
                  for item in manifest.get("frozen_artifacts", [])))
    if ok:
        print("A4-C1-Prep-2 frozen: OK")
        return 0
    if allow_rebuild:
        print("A4-C1-Prep-2 frozen artifacts changed; explicit rebuild enabled")
        return None
    print("A4-C1-Prep-2 frozen artifacts changed; rebuild refused")
    return 1


def run(allow_rebuild: bool = False) -> int:
    existing = _verify_existing(allow_rebuild)
    if existing is not None:
        return existing
    OUT.mkdir(parents=True, exist_ok=True)
    prerequisites = {
        "matlab": MATLAB.is_file(), "pharlap_mex": MEX.is_file(),
        "matlab_script": MATLAB_SCRIPT.is_file(),
        "prep1_report": PREP1_REPORT.is_file(),
        "prep1_manifest": PREP1_MANIFEST.is_file(),
    }
    if not all(prerequisites.values()):
        _atomic_json(REPORT, {"status": "blocked_prerequisite", "checks": prerequisites})
        return 2
    try:
        input_path, input_hash = _verify_prep1()
    except Exception as error:
        _atomic_json(REPORT, {"status": "blocked_prep1_hash", "error": str(error)})
        return 2

    process = _run_matlab(input_path)
    summary = (json.loads(MATLAB_JSON.read_text(encoding="utf-8"))
               if MATLAB_JSON.is_file() else {})
    records = summary.get("combined_records", [])
    roots = summary.get("roots", [])
    figure = (_plot_scan(summary, {"raw": _sha256(RAW),
                                   "matlab_summary": _sha256(MATLAB_JSON)})
              if records and RAW.is_file() and MATLAB_JSON.is_file() else {})
    if records:
        _write_csv(CSV, records)
    if roots:
        _write_csv(ROOT_CSV, [{**root, "scan_stage": "bounded_root"} for root in roots])
    expected = np.arange(1.0, 89.0 + 0.25, 0.5)
    coarse = summary.get("coarse_records", [])
    coarse_angles = np.array([item.get("elevation_deg", np.nan) for item in coarse])
    root_checks = [
        root.get("status") == "ground"
        and abs(float(root.get("along_track_residual_km", np.inf))) <= 1.0
        and abs(float(root.get("cross_track_km", np.inf))) <= 1.0
        and int(root.get("exitflag", 0)) > 0
        for root in roots
    ]
    checks = {
        **prerequisites,
        "prep1_frozen_hash": True,
        "matlab_returncode": process.returncode == 0,
        "matlab_status": summary.get("status") == "passed_pharlap_full_domain_scan",
        "complete_coarse_domain": len(coarse) == expected.size
                                  and np.allclose(coarse_angles, expected, atol=1e-12),
        "no_catastrophic_error": int(summary.get("status_counts", {}).get(
            "catastrophic_error", 0)) == 0,
        "no_grid_exit_after_domain_expansion": int(summary.get(
            "status_counts", {}).get("grid_exit", 0)) == 0,
        "adaptive_scan_present": len(summary.get("adaptive_records", [])) > 0,
        "ground_only_brackets": bool(summary.get("checks", {}).get(
            "ground_only_brackets", False)),
        "at_least_two_roots": len(roots) >= 2,
        "all_roots_ground_and_within_1km": bool(root_checks) and all(root_checks),
        "raw_output": RAW.is_file() and RAW.stat().st_size > 0,
        "legacy_angles_not_used": bool(summary.get("checks", {}).get(
            "legacy_angles_not_used", False)),
        "visualization": bool(figure) and all(figure["checks"].values()),
    }
    status = "passed_pharlap_full_domain_scan" if all(checks.values()) else "failed"
    report = {
        "schema_version": "a4-c1-prep2-validation-v1", "status": status,
        "case_semantics": "pharlap_full_domain_scan_not_ah3d_or_c1_external_validation",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks, "prep1_input": {"path": str(input_path), "sha256": input_hash},
        "hashes": {
            "python": _sha256(Path(__file__)), "matlab": _sha256(MATLAB_SCRIPT),
            "pharlap_mex": _sha256(MEX),
            "raw": _sha256(RAW) if RAW.is_file() else "",
            "matlab_summary": _sha256(MATLAB_JSON) if MATLAB_JSON.is_file() else "",
        },
        "summary": summary,
        "artifacts": {"raw_mat": str(RAW), "matlab_summary": str(MATLAB_JSON),
                      "scan_csv": str(CSV), "root_csv": str(ROOT_CSV),
                      "matlab_log": str(LOG), "figure": figure},
        "limitations": [
            "Prep-2 is PHaRLAP-only and does not compare AH3D trajectories",
            "Root candidates are full-domain forward-map roots, not proof of complex-field mode completeness",
            "High/low variational semantics are assigned only after later branch-signature audit",
        ],
    }
    _atomic_json(REPORT, report)
    if status != "passed_pharlap_full_domain_scan":
        print(f"A4-C1-Prep-2: FAIL report={REPORT}")
        return 1
    artifacts = [RAW, MATLAB_JSON, CSV, ROOT_CSV, FIGURE, MATLAB_SCRIPT, MEX,
                 Path(__file__), input_path]
    _atomic_json(MANIFEST, {
        "schema_version": "a4-c1-prep2-freeze-v1", "status": status,
        "case_semantics": report["case_semantics"],
        "validation_report": {"path": str(REPORT), "sha256": _sha256(REPORT)},
        "frozen_artifacts": [
            {"path": str(path), "sha256": _sha256(path)} for path in artifacts],
    })
    print(f"A4-C1-Prep-2: OK roots={len(roots)} report={REPORT}")
    return 0


if __name__ == "__main__":
    arguments = sys.argv[1:]
    if any(argument != "--rebuild-frozen" for argument in arguments):
        raise SystemExit("usage: _run_a4_c1_prep2.py [--rebuild-frozen]")
    raise SystemExit(run(allow_rebuild="--rebuild-frozen" in arguments))
