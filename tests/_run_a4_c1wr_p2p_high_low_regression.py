"""Current ECEF P2P regression for the frozen Coleman Chapman high/low roots."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models.ne_provider_geocentric import ChapmanGeocentricNeProvider3D  # noqa: E402
from models.ray_tracer_ah_3d_ecef import (  # noqa: E402
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)
from utils import rx_enu_residual, spherical_destination  # noqa: E402

OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr" / "p2p_regression"
REFERENCE = ROOT / "output2" / "a3r_benchmark" / "a3r_p0c_reference_preparation.json"
RESULT = OUT / "a4_c1wr_p2p_high_low_regression.json"
FIGURE = OUT / "a4_c1wr_p2p_high_low_paths.png"
MANIFEST = OUT / "a4_c1wr_p2p_high_low_manifest.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n",
        encoding="utf-8")


def _range_height(ray):
    coordinates = np.asarray(
        ray["trajectory_geocentric_lat_lon_height"], dtype=float)
    distance = 6371.0 * np.radians(coordinates[:, 0] - 30.0)
    return distance, coordinates[:, 2]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    if reference.get("status") != "p0c_1_passed":
        raise RuntimeError("frozen P0C reference is not passed")
    target = float(reference["selected_target_distance_km"])
    roots = [float(value) for value in reference["independent_roots_deg"]]
    rx_lat, rx_lon = spherical_destination(30.0, 120.0, 0.0, target, 6371.0)
    tx = np.array([30.0, 120.0, 0.0])
    rx = np.array([rx_lat, rx_lon, 0.0])
    provider = ChapmanGeocentricNeProvider3D(
        6371.0, foF2_MHz=7.0, hmF2_km=350.0, ymF2_km=100.0)
    model = RefractiveIndexAH3DECEF(
        provider, earth_radius_km=6371.0, freq_MHz=10.0,
        wave_mode="O", gyro_MHz=0.0, b_ecef_hat=(0.0, 0.0, 1.0),
        integration_profile="reference")
    rays = []
    for label, elevation in zip(("low_ray", "high_ray"), roots):
        ray = trace_ray_ah3d_ecef(
            model, tx_lat_deg=tx[0], tx_lon_deg=tx[1],
            tx_height_km=tx[2], elevation_deg=elevation, bearing_deg=0.0,
            h_stop_km=600.0, max_path_km=3000.0, profile="reference")
        terminal = np.asarray(ray["trajectory_ecef_km"], dtype=float)[-1]
        residual = rx_enu_residual(
            terminal, rx[0], rx[1], rx[2], model.earth_radius_km)
        ray.update({
            "p2p_seed_id": label,
            "p2p_landing_error_km": float(np.linalg.norm(residual)),
            "p2p_residual_east_north_km": residual.tolist(),
            "p2p_search_scope": "frozen_roots_forward_regression_no_resolve",
        })
        rays.append(ray)
    by_seed = {ray["p2p_seed_id"]: ray for ray in rays}

    rows = []
    for label, reference_elevation in zip(("low_ray", "high_ray"), roots):
        ray = by_seed.get(label)
        rows.append({
            "ray_label": label,
            "found": ray is not None,
            "reference_elevation_deg": reference_elevation,
            "elevation_deg": None if ray is None else float(ray["elevation_tx_deg"]),
            "bearing_deg": None if ray is None else float(ray["bearing_tx_deg"]),
            "termination_status": None if ray is None else ray["termination_status"],
            "turning_point_count": None if ray is None else int(ray["turning_point_count"]),
            "landing_error_km": None if ray is None else float(ray["p2p_landing_error_km"]),
            "apex_height_km": None if ray is None else float(ray["apex_height_km"]),
            "geom_path_km": None if ray is None else float(ray["geom_path_km"]),
            "group_path_km": None if ray is None else float(ray["group_path_km"]),
            "tau_ms": None if ray is None else float(ray["tau_ms"]),
            "hamiltonian_residual_max": (
                None if ray is None else float(ray["hamiltonian_residual_max"])),
        })
    gates = {
        "two_unique_roots": len(rays) == 2 and set(by_seed) == {"low_ray", "high_ray"},
        "both_ground_single_turning": all(
            row["termination_status"] == "ground"
            and row["turning_point_count"] == 1 for row in rows),
        "landing_error_below_5_km": all(
            row["landing_error_km"] is not None
            and row["landing_error_km"] < 5.0 for row in rows),
        "frozen_seed_neighborhood": all(
            row["elevation_deg"] is not None
            and abs(row["elevation_deg"] - row["reference_elevation_deg"]) < 0.1
            for row in rows),
        "high_low_topology_order": (
            rows[0]["apex_height_km"] is not None
            and rows[1]["apex_height_km"] is not None
            and rows[0]["apex_height_km"] < rows[1]["apex_height_km"]
            and rows[0]["tau_ms"] < rows[1]["tau_ms"]),
        "hamiltonian_residual": all(
            row["hamiltonian_residual_max"] is not None
            and row["hamiltonian_residual_max"] <= 1.0e-8 for row in rows),
        "legacy_solver_not_used": True,
    }
    status = "passed_current_ecef_p2p_high_low" if all(gates.values()) else (
        "failed_current_ecef_p2p_high_low")

    height = np.linspace(0.0, 500.0, 251)
    profile = np.array([
        provider.get_Ne_geocentric(30.0, 120.0, value) for value in height])
    background = np.tile(profile[:, None], (1, 350))
    figure, axis = plt.subplots(figsize=(10, 6))
    image = axis.imshow(
        background, origin="lower", aspect="auto", cmap="jet", alpha=0.5,
        extent=(0.0, target * 1.08, height[0], height[-1]))
    figure.colorbar(image, ax=axis, label="Ne (m^-3)")
    for label, color in (("low_ray", "#1f77b4"), ("high_ray", "#d62728")):
        if label in by_seed:
            distance, ray_height = _range_height(by_seed[label])
            axis.plot(distance, ray_height, color=color, linewidth=1.8, label=label)
    axis.axvline(target, color="black", linestyle="--", label="RX")
    axis.set(
        xlim=(0.0, target * 1.08), ylim=(0.0, 500.0),
        xlabel="Great-circle range (km)", ylabel="Physical height (km)",
        title="A4-C1WR current ECEF P2P Coleman high/low regression\n"
              "one-way | 10 MHz | B=0 | no scattering/power")
    axis.legend()
    figure.tight_layout()
    figure.savefig(FIGURE, dpi=180)
    plt.close(figure)

    result = {
        "schema_version": "a4-c1wr-current-ecef-p2p-high-low-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "semantics": (
            "current_ECEF_forward_regression_at_frozen_Coleman_P0C_roots;"
            "no_root_resolve_and_not_full_domain_mode_completeness"),
        "target_distance_km": target,
        "tx_geocentric": tx.tolist(),
        "rx_geocentric": rx.tolist(),
        "rows": rows,
        "gates": gates,
        "frozen_reference_sha256": _sha256(REFERENCE),
        "legacy_benchmark_status": (
            "obsolete_after_legacy_isolation_AH3D_SOLVER_ecef_wgs84_rejected"),
    }
    _write_json(RESULT, result)
    _write_json(MANIFEST, {
        "schema_version": "a4-c1wr-current-ecef-p2p-manifest-v1",
        "status": status,
        "source_sha256": {
            str(Path(__file__).relative_to(ROOT)): _sha256(Path(__file__)),
            str(REFERENCE.relative_to(ROOT)): _sha256(REFERENCE),
            "models/point_to_point_ah_3d_ecef.py": _sha256(
                ROOT / "models" / "point_to_point_ah_3d_ecef.py"),
            "models/ray_tracer_ah_3d_ecef.py": _sha256(
                ROOT / "models" / "ray_tracer_ah_3d_ecef.py"),
        },
        "artifact_sha256": {
            str(RESULT.relative_to(ROOT)): _sha256(RESULT),
            str(FIGURE.relative_to(ROOT)): _sha256(FIGURE),
        },
    })
    print(json.dumps({"status": status, "gates": gates}, ensure_ascii=True))
    return 0 if status.startswith("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
