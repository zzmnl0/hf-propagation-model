"""A3Rnew-P3 literature and analytical validation.

This runner does not tune production inputs or modify the P2P solver.  It
combines frozen P0C evidence with independent analytical AH checks and emits a
machine-readable record whose reproduction semantics are explicit.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from models.dispersion_ah import C_KM_S, appleton_hartree_roots, dispersion_D, select_root
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.ray_tracer_ah_3d import (
    RefractiveIndexAH3D,
    group_velocity_ah3d,
    trace_single_ray_ah_3d,
)
from tests._benchmark_p2p_ah3d import (
    CASES,
    _ColemanChapmanProvider,
    _build_case,
)
from models.scenario_provider_factory import ScenarioProviderFactory


OUT = ROOT / "output2" / "a3rnew_p3"
P0C = ROOT / "output2" / "a3r_benchmark"
SCHEMA = "a3rnew-p3-v1"
FREQUENCIES_MHZ = (8.0, 9.5, 11.0, 12.5, 14.0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _finite(value) -> bool:
    return value is not None and np.isfinite(float(value))


def _analytical_checks() -> dict:
    omega = 2.0 * np.pi * cfg.FREQ_MHZ * 1e6
    direction = np.array([1.0, 0.0, 0.2])
    b_hat = np.array([0.0, 0.0, 1.0])
    vacuum_roots = appleton_hartree_roots(omega, 0.0, direction, b_hat, 1.2)

    x = np.array([0.0, 10.0])
    z = np.array([100.0, 110.0])
    vacuum_n3d = RefractiveIndex3D(np.zeros((2, 2)), x, z, cfg.FREQ_MHZ,
                                   earth_flat=False)
    vacuum = RefractiveIndexAH3D(
        vacuum_n3d, cfg.FREQ_MHZ, "O", {**cfg.GEOMAG, "fH_MHz": 0.0},
        use_dispersion_v2=True,
    )
    k_vac = np.array([omega / C_KM_S, 0.0, 0.0])
    vg_vac = group_velocity_ah3d(5.0, 0.0, 105.0, k_vac, vacuum)

    ne_iso = 0.35 * (cfg.FREQ_MHZ * 1e6 / cfg.K_FP) ** 2
    iso_n3d = RefractiveIndex3D(np.full((2, 2), ne_iso), x, z, cfg.FREQ_MHZ,
                                earth_flat=False)
    iso_models = {
        mode: RefractiveIndexAH3D(
            iso_n3d, cfg.FREQ_MHZ, mode, {**cfg.GEOMAG, "fH_MHz": 0.0},
            use_dispersion_v2=True,
        ) for mode in ("O", "X")
    }
    iso_n2 = {mode: model.n2_raw(5.0, 0.0, 105.0, direction)
              for mode, model in iso_models.items()}
    n_iso = np.sqrt(iso_n2["O"])
    k_iso = np.array([n_iso * omega / C_KM_S, 0.0, 0.0])
    vg_iso = group_velocity_ah3d(5.0, 0.0, 105.0, k_iso, iso_models["O"])

    cutoff_ne = (cfg.FREQ_MHZ * 1e6 / cfg.K_FP) ** 2
    cutoff_states = {}
    for label, scale in (("below", 0.995), ("above", 1.005)):
        field = RefractiveIndex3D(
            np.full((2, 2), scale * cutoff_ne), x, z, cfg.FREQ_MHZ,
            earth_flat=False,
        )
        model = RefractiveIndexAH3D(
            field, cfg.FREQ_MHZ, "O", {**cfg.GEOMAG, "fH_MHz": 0.0},
            use_dispersion_v2=True,
        )
        cutoff_states[label] = model.branch_status(5.0, 0.0, 105.0, k_vac)

    previous = {"O": None, "X": None}
    scan = {"O": [], "X": []}
    for freq in np.linspace(8.0, 14.0, 121):
        roots = appleton_hartree_roots(
            2.0 * np.pi * freq * 1e6, 4.0e11, [1.0, 0.2, 0.3],
            [0.45, -0.15, -0.88], 1.2,
        )
        for mode in ("O", "X"):
            root = select_root(roots, mode, previous[mode])
            scan[mode].append(root)
            previous[mode] = root.n2

    checks = {
        "vacuum_n2_unity": all(abs(root.real - 1.0) <= 1e-12
                                for root in vacuum_roots.values()),
        "vacuum_group_speed_c": abs(np.linalg.norm(vg_vac) / C_KM_S - 1.0) <= 1e-7,
        "b0_ox_degenerate": abs(iso_n2["O"] - iso_n2["X"]) <= 1e-12,
        "isotropic_group_speed": abs(np.linalg.norm(vg_iso) / (C_KM_S * n_iso) - 1.0) <= 2e-6,
        "cutoff_classification": cutoff_states == {"below": "near_cutoff", "above": "evanescent"},
        "frequency_branch_identity": all(
            root.branch_id == mode and root.status != "branch_lost"
            for mode, roots in scan.items() for root in roots
        ),
        "frequency_positive_branch_gap": all(
            root.branch_gap > 0.0 for root in scan["O"] + scan["X"]
        ),
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "vacuum_group_speed_ratio": float(np.linalg.norm(vg_vac) / C_KM_S),
        "isotropic_group_speed_ratio": float(np.linalg.norm(vg_iso) / (C_KM_S * n_iso)),
        "b0_n2": iso_n2,
        "cutoff_states": cutoff_states,
        "dispersion_frequency_scan": {
            "range_MHz": [8.0, 14.0], "point_count": 121,
            "minimum_branch_gap": float(min(
                root.branch_gap for root in scan["O"] + scan["X"])),
        },
    }


def _coleman_evidence() -> tuple[dict, list[dict]]:
    prep = _json(P0C / "a3r_p0c_reference_preparation.json")
    closure = _json(P0C / "a3r_p0c_closure.json")
    repeated = closure["repeats"]
    matches = [match for run in repeated for match in run["matches"]]
    path_errors = []
    formal_rays = []
    for run in repeated:
        independent = {item["ray_label"]: item for item in run["independent"]}
        for match in run["matches"]:
            reference = independent[match["ray_label"]]["reference"]
            path_errors.append(abs(
                float(match["production"]["geom_path_km"])
                - float(reference["geom_path_km"])))
        for item in run["independent"]:
            formal_rays.extend((item["research"], item["reference"]))
    checks = {
        "two_modes_each_repeat": all(len(run["production"]) == 2 for run in repeated),
        "all_modes_matched": all(match["matched"] for match in matches),
        "elevation_tolerance_deg": max(match["beta_error_deg"] for match in matches) <= 0.05,
        "height_tolerance_km": max(match["height_error_km"] for match in matches) <= 1.0,
        "delay_tolerance_ms": max(match["tau_error_ms"] for match in matches) <= 0.01,
        "geometric_path_tolerance_km": max(path_errors) <= 1.0,
        "formal_integration_convergence": all(
            ray["branch_status"] == "ground"
            and ray["convergence_status"] == "passed"
            and ray["hamiltonian_residual_max"] <= 1e-6
            for ray in formal_rays
        ),
        "ground_and_hit_tolerance": all(
            mode["branch_status"] == "ground" and mode["p2p_landing_error_km"] < 5.0
            for run in repeated for mode in run["production"]
        ),
    }

    spec = next(item for item in CASES if item["name"] == "coleman_chapman_high_low")
    provider, n3d, _ = _build_case(ScenarioProviderFactory(cfg), spec)
    model = RefractiveIndexAH3D(
        n3d, 10.0, "O", geomag={**cfg.GEOMAG, "fH_MHz": 0.0},
        integration_profile="research",
    )
    trajectories = []
    for label, beta in zip(("low_ray", "high_ray"), prep["independent_roots_deg"]):
        ray = trace_single_ray_ah_3d(cfg.TX_POS, float(beta), 0.0, model, 10.0, t=0.0)
        trajectory = np.asarray(ray["trajectory_3d"], dtype=float)
        trajectories.append({"label": label, "beta_deg": float(beta), "trajectory": trajectory})

    return {
        "reproduction_semantics": "literature_guided_parameterized_reproduction",
        "strict_reproduction": False,
        "reason_not_strict": (
            "Coleman reports an approximately 1100 km case; the frozen receiver "
            "distance 1185.848260 km was selected independently inside the computed two-root interval"
        ),
        "frozen_inputs": {
            "foF2_MHz": 7.0, "hmF2_km": 350.0, "ymF2_km": 100.0,
            "frequency_MHz": 10.0, "geomagnetic_field": "B=0",
            "collision_frequency_Hz": 0.0,
            "target_distance_km": float(prep["selected_target_distance_km"]),
            "provider_hash": provider.get_metadata()["provider_hash"],
        },
        "reference_roots_deg": prep["independent_roots_deg"],
        "repeat_count": len(repeated),
        "checks": checks,
        "passed": all(checks.values()),
        "maximum_errors": {
            "elevation_deg": float(max(match["beta_error_deg"] for match in matches)),
            "height_km": float(max(match["height_error_km"] for match in matches)),
            "delay_ms": float(max(match["tau_error_ms"] for match in matches)),
            "geometric_path_km": float(max(path_errors)),
            "hamiltonian_residual": float(max(
                ray["hamiltonian_residual_max"] for ray in formal_rays)),
        },
    }, trajectories


def _nosikov_evidence() -> dict:
    return {
        "reproduction_semantics": "partial_reproduction_with_assumptions",
        "comparison_status": "partial_not_matched",
        "strict_reproduction": False,
        "available_inputs": {
            "N1_cm-3": 0.2e6, "N2_cm-3": 1.0e6,
            "z1_km": 110.0, "z2_km": 300.0,
            "delta_z1_km": 30.0, "delta_z2_km": 75.0,
            "frequency_MHz": 12.0,
        },
        "missing_for_unique_reproduction": [
            "complete transmitter and receiver endpoint geometry",
            "unambiguous Earth geometry and coordinate convention",
            "supplementary three-dimensional profile inputs",
            "complete boundary and initialization settings",
        ],
        "assumptions_used": [
            "published two-layer density parameters are recorded without inverse fitting",
            "no endpoint geometry is invented and no numerical match is claimed",
        ],
        "passed": True,
    }


def _frequency_evidence() -> tuple[dict, list[dict]]:
    provider = _ColemanChapmanProvider()
    ne = provider.get_Ne_grid(cfg.BG_X, cfg.BG_Z)
    records = []
    for freq in FREQUENCIES_MHZ:
        n3d = RefractiveIndex3D(
            ne, cfg.BG_X, cfg.BG_Z, freq, earth_flat=True, ne_provider=provider,
        )
        for mode in ("O", "X"):
            model = RefractiveIndexAH3D(
                n3d, freq, mode, geomag=cfg.GEOMAG,
                integration_profile="smoke", use_dispersion_v2=True,
            )
            ray = trace_single_ray_ah_3d(cfg.TX_POS, 35.339807513, 0.0,
                                         model, freq, t=0.0)
            record = {
                "frequency_MHz": freq, "branch_id": mode,
                "reported_wave_mode": ray.get("wave_mode"),
                "branch_status": ray.get("branch_status"),
                "branch_gap_min": ray.get("branch_gap_min"),
                "x_land_km": ray.get("x_land_km") if _finite(ray.get("x_land_km")) else None,
                "h_reflect_km": ray.get("h_reflect_km") if _finite(ray.get("h_reflect_km")) else None,
                "tau_ms": ray.get("tau_ms") if _finite(ray.get("tau_ms")) else None,
            }
            records.append(record)
    allowed = {"ground", "upper_boundary", "escaped", "near_cutoff", "evanescent", "turning_point",
               "branch_lost", "resonance_excluded", "step_failed", "max_steps"}
    checks = {
        "five_frequencies_both_branches": len(records) == 10,
        "no_silent_branch_switch": all(
            item["reported_wave_mode"] == item["branch_id"] for item in records),
        "explicit_physical_states": all(item["branch_status"] in allowed for item in records),
    }
    return {
        "launch_semantics": "fixed_launch_angle_ray_level_branch_audit_not_p2p_mode_continuation",
        "frequencies_MHz": list(FREQUENCIES_MHZ),
        "physical_state_transitions_allowed": True,
        "checks": checks, "passed": all(checks.values()), "records": records,
    }, records


def _plots(trajectories: list[dict], frequency_records: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for item in trajectories:
        tr = item["trajectory"]
        ax.plot(tr[:, 0], tr[:, 2], lw=1.8,
                label=f"{item['label']} ({item['beta_deg']:.3f} deg)")
    ax.set(xlabel="Ground-range coordinate x (km)", ylabel="Height z (km)",
           title="Coleman-guided Chapman F2 high/low rays; B=0; one-way AH3D")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "coleman_high_low_trajectories.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
    for mode, marker in (("O", "o"), ("X", "s")):
        rows = [row for row in frequency_records if row["branch_id"] == mode]
        freq = [row["frequency_MHz"] for row in rows]
        height = [np.nan if row["h_reflect_km"] is None else row["h_reflect_km"] for row in rows]
        delay = [np.nan if row["tau_ms"] is None else row["tau_ms"] for row in rows]
        axes[0].plot(freq, height, marker=marker, label=mode)
        axes[1].plot(freq, delay, marker=marker, label=mode)
    axes[0].set_ylabel("Apex/terminal height (km)")
    axes[1].set_ylabel("Accumulated one-way group delay (ms)")
    axes[1].set_xlabel("Frequency (MHz)")
    axes[0].set_title("Fixed-launch AH3D O/X frequency audit; state transitions retained")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "ah3d_ox_frequency_trend.png", dpi=180)
    plt.close(fig)


def main() -> int:
    required = [
        P0C / "a3r_p0c_reference_preparation.json",
        P0C / "a3r_p0c_closure.json",
        ROOT / "output2" / "a3rnew_p2l" / "p2_6l_reference_acceptance.json",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing frozen prerequisite: " + ", ".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)

    analytical = _analytical_checks()
    coleman, trajectories = _coleman_evidence()
    nosikov = _nosikov_evidence()
    frequency, frequency_records = _frequency_evidence()
    _plots(trajectories, frequency_records)

    p2 = _json(required[2])
    checks = {
        "analytical_limits": analytical["passed"],
        "coleman_parameterized_reproduction": coleman["passed"],
        "nosikov_semantics_complete": nosikov["passed"],
        "ray_frequency_branch_audit": frequency["passed"],
        "p2_prerequisite": p2.get("status") == "reference_case_limited_pass",
    }
    result = {
        "schema_version": SCHEMA,
        "stage": "A3Rnew-P3",
        "status": "passed" if all(checks.values()) else "failed",
        "evidence_scope": "internal_analytical_and_literature_guided_validation",
        "not_external_validation": "PHaRLAP comparison remains A4",
        "checks": checks,
        "analytical_limits": analytical,
        "coleman_2011": coleman,
        "nosikov_2020": nosikov,
        "ray_frequency_audit": frequency,
        "thresholds": {
            "landing_error_km_lt": 5.0,
            "hamiltonian_residual_le": 1e-6,
            "integration_landing_delta_km_le": 1.0,
            "integration_height_delta_km_le": 0.5,
            "integration_delay_relative_le": 1e-3,
            "coleman_elevation_error_deg_le": 0.05,
            "coleman_height_error_km_le": 1.0,
            "coleman_delay_error_ms_le": 0.01,
        },
        "provenance": {
            "prerequisite_hashes": {str(path.relative_to(ROOT)): _sha256(path) for path in required},
            "code_hashes": {
                path: _sha256(ROOT / path) for path in (
                    "models/dispersion_ah.py", "models/ray_tracer_ah_3d.py",
                    "models/ray_tracer_ah_3d_v2.py", "models/point_to_point_ah_3d.py",
                    "tests/_run_a3rnew_p3_validation.py",
                )
            },
            "configuration_semantics": "frozen inputs; one deterministic run; no parameter tuning",
        },
        "residual_risks": [
            "P2-6H complex-scene high-density completeness audit is deferred",
            "Nosikov 2020 cannot be uniquely reproduced from available inputs",
            "near-cutoff tunneling and mode conversion require wave methods",
            "AH3D remains Hamiltonian ray tracing rather than a full-wave solver",
        ],
        "next_gate": "P5new accuracy report before A4 PHaRLAP",
    }
    path = OUT / "a3rnew_p3_validation.json"
    path.write_text(json.dumps(
        result, indent=2, ensure_ascii=True,
        default=lambda value: value.item() if isinstance(value, np.generic) else str(value),
    ), encoding="utf-8")
    print(f"A3Rnew-P3 status={result['status']}")
    print(f"result={path}")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
