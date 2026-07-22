"""Lightweight C0R-6 production-entry and output-semantics verification."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
import main2
from models.ne_provider_geocentric import (
    ChapmanGeocentricNeProvider3D,
    GeocentricNeProvider3D,
)
from models.unified_hybrid_pipeline import (
    UnifiedPropagationRequest,
    run_unified_pipeline,
)
from viz.plot_utils import plot_tx_ray_fan, save_modes_csv


OUT = ROOT / "output2" / "a4_pharlap"
REPORT = OUT / "c0r6_entry_validation.json"
FIGURE = OUT / "c0r6_ecef_entry_smoke.png"
CSV_DIR = OUT / "c0r6_csv"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    production_request = main2._request_for("tid")
    chapman = ChapmanGeocentricNeProvider3D(cfg.RE_KM)
    rt = {**cfg.RT_3D, "beta_min": 35.0, "beta_max": 35.0,
          "n_fan": 1, "max_path_km": 3000.0}
    request = UnifiedPropagationRequest(
        scenario="c0r6_ecef_smoke", provider=chapman,
        tx_geocentric=cfg.TX_GEOCENTRIC,
        rx_geocentric=cfg.RX_GEOCENTRIC,
        bearing_deg=cfg.LINK_BEARING_DEG,
        earth_radius_km=cfg.RE_KM, solver="ah3d_ecef",
        rt_config=rt, make_ray_fan=True, make_p2p=False, make_pd=False)
    result = run_unified_pipeline(request)
    ray = result.rays3d[0]

    x = np.linspace(0.0, 1500.0, 31)
    z = np.linspace(60.0, 600.0, 28)
    ne = chapman.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, x, z)
    link_info = {
        "TX_LAT": cfg.TX_LAT, "TX_LON": cfg.TX_LON,
        "freq_MHz": cfg.FREQ_MHZ, "Pt_W": cfg.PT_W,
        "RX_km": cfg.RX_RANGE, "IRI_DT": cfg.IRI_DT,
        "IRI_LAT": cfg.IRI_LAT, "IRI_LON": cfg.IRI_LON,
        "earth_model": "sphere", "earth_radius_km": cfg.RE_KM,
        "coordinate_model": "geocentric lat/lon/height -> ECEF state",
        "oneway_or_twoway": "one_way", "reflection_coupling": "on",
        "ray_pe_coupling": "off", "scatter_coupling": "diagnostic_only",
    }
    fig, _ = plot_tx_ray_fan(
        result.rays3d, ne, x, z, tx_pos=(0.0, 0.0),
        rx_pos=(cfg.RX_RANGE, 0.0), link_info=link_info,
        title="C0R-6 ECEF AH3D entry smoke")
    fig.savefig(FIGURE, dpi=100, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)

    mode = dict(ray)
    mode.update({
        "label": "O-smoke", "Pr_W": 1.0, "Pr_dBW": 0.0,
        "delta_tau_ms": 0.0, "h_reflect_km": ray.get("apex_height_km", 0.0),
        "provider_name": result.provider_metadata["name"],
        "provider_components": "background", "solver_chain": "ah3d_ecef",
        "observable_type": "communication_oneway", "delay_key": "tau_ms",
    })
    tau = np.array([mode["tau_ms"] - 0.01, mode["tau_ms"], mode["tau_ms"] + 0.01])
    csv_path = Path(save_modes_csv(
        "c0r6_smoke", [mode], tau, np.array([0.5, 1.0, 0.5]),
        cfg.FREQ_MHZ, str(CSV_DIR)))
    with csv_path.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))

    legacy_rejected = False
    try:
        run_unified_pipeline(UnifiedPropagationRequest(
            scenario="reject", provider=chapman, solver="iso3d_legacy",
            make_ray_fan=False, make_p2p=False, make_pd=False))
    except ValueError:
        legacy_rejected = True

    image = mpimg.imread(FIGURE)
    checks = {
        "formal_config_TX_RX_geocentric": (
            len(cfg.TX_GEOCENTRIC) == 3 and len(cfg.RX_GEOCENTRIC) == 3
            and cfg.EARTH_MODEL["coordinate_system"] == "geocentric_spherical"),
        "external_Ne_contract_explicit": all(
            key in cfg.UPDATE2["EXTERNAL_NE_PROVIDER"] for key in
            ("path", "coordinate_system", "axis_order", "units", "earth_radius_km")),
        "main2_uses_formal_provider": (
            isinstance(production_request.provider, GeocentricNeProvider3D)
            and production_request.solver == "ah3d_ecef"
            and production_request.tx_geocentric == cfg.TX_GEOCENTRIC
            and production_request.rx_geocentric == cfg.RX_GEOCENTRIC),
        "legacy_solver_has_no_fallback": legacy_rejected,
        "ray_has_only_explicit_ECEF_geocentric_paths": (
            "trajectory_ecef_km" in ray
            and "trajectory_geocentric_lat_lon_height" in ray
            and "trajectory_3d" not in ray),
        "traceability_semantics_complete": all(
            key in result.traceability for key in
            ("earth_model", "coordinate_model", "coordinate_semantics",
             "oneway_or_twoway", "ray_pe_coupling", "scatter_power_coupling")),
        "csv_semantics_complete": all(row.get(key) not in (None, "") for key in
            ("earth_model", "earth_radius_km", "coordinate_model",
             "coordinate_semantics", "plot_coordinate_semantics",
             "oneway_or_twoway")),
        "annotated_figure_nonblank": (
            FIGURE.stat().st_size > 10000 and np.isfinite(image).all()
            and float(np.std(image)) > 1.0e-3),
    }
    report = {
        "schema_version": "a4-c0r6-entry-v1",
        "stage": "A4-C0R-6",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "lightweight_production_entry_coordinate_and_output_contract",
        "checks": checks,
        "metrics": {
            "ray_count": len(result.rays3d),
            "termination_status": ray["termination_status"],
            "figure_bytes": FIGURE.stat().st_size,
            "figure_pixel_std": float(np.std(image)),
        },
        "artifacts": {"figure": str(FIGURE), "csv": str(csv_path)},
        "source_sha256": {
            name: _sha256(ROOT / path) for name, path in {
                "config": "config.py", "main": "main.py", "main2": "main2.py",
                "pipeline": "models/unified_hybrid_pipeline.py",
                "provider": "models/ne_provider_geocentric.py",
                "test": "tests/_verify_a4_c0r6_entry_migration.py",
            }.items()
        },
        "claim_boundary": (
            "Production entry, formal coordinate flow, ECEF ray smoke, and output "
            "semantics only; not PHaRLAP trajectory agreement or full-domain P2P."),
    }
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(" {}: {}".format(name, "OK" if passed else "FAIL"))
    failures = sum(not value for value in checks.values())
    print("A4-C0R-6: {} fail(s)".format(failures))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
