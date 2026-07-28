"""Freeze A4-C1T solver screening and FARR fallback readiness artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1t"
FARR_ROOT = Path(
    r"D:\0文件\0文献\电离层短波传播\reference-code"
    r"\farr-Release\farr-Release"
)
AU_RAY_FILES = [
    ROOT / "sources" / "au_ray_README.pdf",
    ROOT / "sources" / "research_au_ray_download.txt",
]
PLAN = ROOT / "word" / "a4_c1t_third_solver_screening_and_farr_fallback_plan.txt"
REPORT = ROOT / "word" / "a4_c1t_third_solver_screening_report.txt"
FINAL_C1 = (
    ROOT / "output2" / "a4_pharlap" / "c1" / "final"
    / "a4_c1_validation.json"
)
C1Q_FINAL = (
    ROOT / "output2" / "a4_pharlap" / "c1" / "c1q" / "c1q_6"
    / "c1q_6_final_attribution.json"
)

GATES = [
    "full_solver_source",
    "usable_license",
    "independent_hamiltonian",
    "canonical_equations_auditable",
    "b0_collision0_10mhz",
    "external_frozen_ne",
    "exact_c1q_geometry_angles",
    "full_trajectory_state",
    "controlled_build",
]

CANDIDATES = [
    {
        "candidate": "AU-Ray public demo",
        "source_ref": "Zenodo AU-RAY DEMO; Hirvonen et al. (2025)",
        "full_solver_source": False,
        "usable_license": False,
        "independent_hamiltonian": True,
        "canonical_equations_auditable": False,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": False,
        "exact_c1q_geometry_angles": False,
        "full_trajectory_state": True,
        "controlled_build": False,
        "reason": "Compiled restricted Linux demo; no C++ core source and fixed inputs.",
    },
    {
        "candidate": "FARR",
        "source_ref": "Local GPLv3 release; Green et al. (2025)",
        "full_solver_source": True,
        "usable_license": True,
        "independent_hamiltonian": False,
        "canonical_equations_auditable": False,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": True,
        "exact_c1q_geometry_angles": False,
        "full_trajectory_state": False,
        "controlled_build": False,
        "reason": "Auditable 3D FDTD full-wave solver, not a Hamiltonian ray solver.",
    },
    {
        "candidate": "IONORT public archive",
        "source_ref": "github.com/blair3sat/ionosphere-rt@6c2f73c",
        "full_solver_source": True,
        "usable_license": False,
        "independent_hamiltonian": True,
        "canonical_equations_auditable": True,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": True,
        "exact_c1q_geometry_angles": True,
        "full_trajectory_state": True,
        "controlled_build": False,
        "reason": "Punch-card-era Fortran archive; no explicit license; documented as non-buildable.",
    },
    {
        "candidate": "AutoNVIS IONORT-style",
        "source_ref": "github.com/n4hy/AutoNVIS@4bfda5a",
        "full_solver_source": True,
        "usable_license": False,
        "independent_hamiltonian": True,
        "canonical_equations_auditable": False,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": False,
        "exact_c1q_geometry_angles": True,
        "full_trajectory_state": True,
        "controlled_build": True,
        "reason": "License is TBD; documented simplified dr/ds and dk/ds semantics are not C1Q canonical equations.",
    },
    {
        "candidate": "HamSCI PyLap",
        "source_ref": "github.com/HamSCI/PyLap@4b32d0a",
        "full_solver_source": False,
        "usable_license": True,
        "independent_hamiltonian": False,
        "canonical_equations_auditable": False,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": True,
        "exact_c1q_geometry_angles": True,
        "full_trajectory_state": True,
        "controlled_build": False,
        "reason": "C/Python wrapper calls PHaRLAP raytrace_3d_ binary and is not independent.",
    },
    {
        "candidate": "andersontodds ray-trace",
        "source_ref": "github.com/andersontodds/ray-trace@e8c50b0",
        "full_solver_source": True,
        "usable_license": False,
        "independent_hamiltonian": True,
        "canonical_equations_auditable": False,
        "b0_collision0_10mhz": False,
        "external_frozen_ne": False,
        "exact_c1q_geometry_angles": False,
        "full_trajectory_state": True,
        "controlled_build": False,
        "reason": "Early-stage whistler/plasmasphere code with meridional-plane assumptions.",
    },
    {
        "candidate": "Jones-Stephenson OT 75-76 listing",
        "source_ref": "Jones and Stephenson (1975), OT Report 75-76",
        "full_solver_source": False,
        "usable_license": False,
        "independent_hamiltonian": True,
        "canonical_equations_auditable": True,
        "b0_collision0_10mhz": True,
        "external_frozen_ne": True,
        "exact_c1q_geometry_angles": True,
        "full_trajectory_state": True,
        "controlled_build": False,
        "reason": "Report listing is not a maintained, licensed, reproducibly buildable source release.",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_signature(root: Path) -> dict:
    digest = hashlib.sha256()
    count = 0
    total_bytes = 0
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root)
        if any(part in {".git", "build", "__pycache__", "googletest"}
               for part in rel.parts):
            continue
        file_hash = sha256(path)
        digest.update(rel.as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        count += 1
        total_bytes += path.stat().st_size
    return {
        "root": str(root),
        "sha256": digest.hexdigest(),
        "file_count": count,
        "total_bytes": total_bytes,
        "excluded": [".git", "build", "__pycache__", "googletest"],
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_admission(rows: list[dict], path: Path) -> None:
    values = np.array([[int(row[gate]) for gate in GATES] for row in rows])
    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    image = ax.imshow(values, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(range(len(GATES)), [gate.replace("_", "\n") for gate in GATES],
                  fontsize=8)
    ax.set_yticks(range(len(rows)), [row["candidate"] for row in rows], fontsize=9)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, "PASS" if values[i, j] else "FAIL",
                    ha="center", va="center", fontsize=7,
                    color="black" if values[i, j] else "white")
    ax.set_title("A4-C1T source-auditable Hamiltonian solver admission matrix")
    ax.set_xlabel("All gates are mandatory; no candidate passed every gate")
    fig.colorbar(image, ax=ax, ticks=[0, 1], label="Gate result")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def farr_levels() -> list[dict]:
    return [
        {
            "level": "pilot_coarse",
            "dx_m": 45.0, "dy_m": 37.5, "dz_m": 45.0,
            "dt_s": 6.51e-8, "nt": 16667, "production_acceptance": False,
        },
        {
            "level": "paper_baseline",
            "dx_m": 30.0, "dy_m": 25.0, "dz_m": 30.0,
            "dt_s": 4.34e-8, "nt": 25000, "production_acceptance": True,
        },
        {
            "level": "selected_refined",
            "dx_m": 24.0, "dy_m": 20.0, "dz_m": 24.0,
            "dt_s": 3.472e-8, "nt": 31250, "production_acceptance": True,
        },
    ]


def add_grid_cost(rows: list[dict]) -> None:
    lengths_m = (150_000.0, 37_500.0, 120_000.0)
    for row in rows:
        shape = [
            int(np.ceil(lengths_m[0] / row["dx_m"])),
            int(np.ceil(lengths_m[1] / row["dy_m"])),
            int(np.ceil(lengths_m[2] / row["dz_m"])),
        ]
        row["nx"], row["ny"], row["nz"] = shape
        row["grid_cells"] = int(np.prod(shape, dtype=np.int64))
        row["simulated_time_ms"] = row["dt_s"] * row["nt"] * 1e3


def plot_farr_cost(rows: list[dict], path: Path) -> None:
    labels = [row["level"] for row in rows]
    cells = [row["grid_cells"] / 1e9 for row in rows]
    steps = [row["nt"] / 1000 for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].bar(x, cells, color=["#4c78a8", "#f58518", "#e45756"])
    axes[0].set_xticks(x, labels, rotation=18)
    axes[0].set_ylabel("Grid cells (billion)")
    axes[0].set_title("Spatial cost")
    axes[1].bar(x, steps, color=["#4c78a8", "#f58518", "#e45756"])
    axes[1].set_xticks(x, labels, rotation=18)
    axes[1].set_ylabel("Time steps (thousand)")
    axes[1].set_title("Temporal cost")
    fig.suptitle("Planned FARR paper-case convergence levels (not executed)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    missing = [str(path) for path in [*AU_RAY_FILES, FARR_ROOT, PLAN, REPORT,
                                      FINAL_C1, C1Q_FINAL] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing frozen input(s): " + ", ".join(missing))

    rows = []
    for source in CANDIDATES:
        row = dict(source)
        row["qualified"] = all(bool(row[gate]) for gate in GATES)
        rows.append(row)
    qualified = [row["candidate"] for row in rows if row["qualified"]]

    levels = farr_levels()
    add_grid_cost(levels)
    local_build_ready = (
        os.name != "nt"
        and shutil.which("mpirun") is not None
        and shutil.which("cmake") is not None
    )
    farr_config = {
        "semantics": "local_full_wave_fallback_not_hamiltonian_c1_closure",
        "execution_status": (
            "ready_for_hpc_execution" if local_build_ready
            else "blocked_by_linux_mpi_hpc_environment"
        ),
        "paper_case": {
            "frequency_mhz": 1.0,
            "launch_elevation_deg": 60.0,
            "array_elements": 16,
            "half_power_beamwidth_deg": 7.16,
            "source_bandwidth_khz": 100.0,
            "domain_km": [150.0, 37.5, 120.0],
            "mpi_processes_reported": 3000,
            "pharlap_model": "2D, B=0, common IRI-2020 density",
            "paper_reference": {
                "time_of_flight_ms": 0.84,
                "first_hop_distance_km": 127.0,
                "arrival_angle_deg": 58.0,
            },
        },
        "levels": levels,
        "acceptance": {
            "paper_reproduction": {
                "time_of_flight_abs_error_ms_max": 0.02,
                "first_hop_distance_abs_error_km_max": 5.0,
                "arrival_angle_abs_error_deg_max": 2.0,
            },
            "baseline_to_refined": {
                "time_of_flight_change_ms_max": 0.01,
                "first_hop_distance_change_km_max": 3.0,
                "arrival_angle_change_deg_max": 1.0,
            },
            "main_lobe_definition": "maximum electromagnetic energy density per time step",
            "ray_rms_thresholds_applicable": False,
        },
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "mpirun": shutil.which("mpirun"),
            "cmake": shutil.which("cmake"),
        },
    }

    matrix_csv = OUT / "c1t_candidate_admission_matrix.csv"
    matrix_json = OUT / "c1t_candidate_admission_matrix.json"
    farr_csv = OUT / "c1t_farr_reproduction_levels.csv"
    farr_json = OUT / "c1t_farr_reproduction_config.json"
    matrix_png = OUT / "c1t_candidate_admission_matrix.png"
    farr_png = OUT / "c1t_farr_fallback_readiness.png"
    write_csv(matrix_csv, rows, ["candidate", *GATES, "qualified",
                                 "source_ref", "reason"])
    matrix_json.write_text(json.dumps({
        "gates": GATES,
        "candidates": rows,
        "qualified_candidates": qualified,
        "screening_conclusion": (
            "qualified_hamiltonian_solver_found" if qualified
            else "no_qualified_solver_use_farr_local_fullwave"
        ),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(farr_csv, levels, list(levels[0].keys()))
    farr_json.write_text(json.dumps(farr_config, indent=2, ensure_ascii=False),
                         encoding="utf-8")
    plot_admission(rows, matrix_png)
    plot_farr_cost(levels, farr_png)

    final_c1 = json.loads(FINAL_C1.read_text(encoding="utf-8"))
    c1q = json.loads(C1Q_FINAL.read_text(encoding="utf-8"))
    status = {
        "stage": "A4-C1T",
        "screening_conclusion": (
            "qualified_hamiltonian_solver_found" if qualified
            else "no_qualified_solver_use_farr_local_fullwave"
        ),
        "final_status": (
            "qualified_hamiltonian_solver_found" if qualified
            else (
                "no_qualified_solver_use_farr_local_fullwave"
                if local_build_ready
                else "validation_blocked_by_build_or_compute_resources"
            )
        ),
        "formal_c1_status": final_c1.get("status"),
        "c1q_attribution": c1q.get("final_attribution"),
        "production_modified": False,
        "a4_c2_admitted": False,
        "next_gate": "manual_decision_required_before_farr_hpc_run_or_a4_c2",
    }
    status_path = OUT / "c1t_status.json"
    status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False),
                           encoding="utf-8")

    local_sources = {
        str(path): sha256(path) for path in AU_RAY_FILES
    }
    local_sources["farr_tree"] = tree_signature(FARR_ROOT)
    source_paths = [
        Path(__file__), PLAN, REPORT, matrix_csv, matrix_json,
        farr_csv, farr_json, matrix_png, farr_png, status_path,
    ]
    manifest = {
        "stage": "A4-C1T",
        "status": status["final_status"],
        "source_sha256": {
            str(path): sha256(path) for path in [Path(__file__), PLAN, REPORT]
        },
        "local_third_party_sources": local_sources,
        "frozen_artifacts": [
            {"path": str(path), "sha256": sha256(path)}
            for path in source_paths[3:]
        ],
        "checks": {
            "no_unqualified_solver_admitted": not qualified,
            "au_ray_demo_not_source_auditable": not rows[0]["full_solver_source"],
            "farr_not_hamiltonian": not rows[1]["independent_hamiltonian"],
            "formal_c1_failure_preserved": (
                final_c1.get("status")
                == "failed_c1_chapman_high_low_external_comparison"
            ),
            "c1q_attribution_preserved": (
                c1q.get("final_attribution")
                == "pharlap_effective_dynamics_not_recoverable"
            ),
            "production_unchanged": not status["production_modified"],
            "c2_blocked": not status["a4_c2_admitted"],
        },
    }
    manifest_path = OUT / "c1t_freeze_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")
    (OUT / "run.log").write_text(
        f"screened_candidates={len(rows)}\n"
        f"qualified_candidates={len(qualified)}\n"
        f"screening_conclusion={status['screening_conclusion']}\n"
        f"final_status={status['final_status']}\n"
        f"farr_execution_status={farr_config['execution_status']}\n",
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
