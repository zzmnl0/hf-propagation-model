"""Acceptance checks for the C0R-2 geocentric NeProvider contract."""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from legacy.models.ne_provider import GridNeProvider, ParameterizedIonosphere3DProvider
from models.ne_provider_geocentric import (
    AXIS_ORDER,
    COORDINATE_SYSTEM,
    ChapmanGeocentricNeProvider3D,
    ExternalGeocentricNeProvider3D,
    GeocentricGridNeProvider3D,
    GeocentricParameterizedIonosphere3DProvider,
    VacuumGeocentricNeProvider3D,
)
from utils import geocentric_to_ecef, spherical_destination


R_KM = 6371.0
REPORT = ROOT / "output2" / "a4_pharlap" / "c0r2_provider_validation.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _constant_geocentric_background(value: float = 1.0e11):
    lat = np.array([28.0, 30.0, 32.0])
    lon = np.array([118.0, 120.0, 122.0])
    height = np.array([100.0, 260.0, 380.0])
    Ne = np.full((lat.size, lon.size, height.size), value)
    return GeocentricGridNeProvider3D(Ne, lat, lon, height, R_KM)


def _external_checks() -> tuple[dict, dict]:
    lat = np.array([-1.0, 0.0, 1.0])
    lon = np.array([179.0, 180.0, 181.0])
    height = np.array([100.0, 200.0])
    Ne_cm3 = np.empty((3, 3, 2))
    for i in range(3):
        for j in range(3):
            for k in range(2):
                Ne_cm3[i, j, k] = 1.0e5 + 100.0 * i + 10.0 * j + k
    rejected = 0
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    directory = REPORT.parent
    created = []
    try:
        directory = Path(directory)
        valid = directory / "_c0r2_contract_valid.npz"
        created.append(valid)
        np.savez(
            valid, Ne=Ne_cm3, latitude_deg=lat, longitude_deg=lon,
            height_km=height, units="cm^-3",
            coordinate_system=COORDINATE_SYSTEM, axis_order=AXIS_ORDER,
            earth_radius_km=R_KM, longitude_periodic_complete=False,
        )
        provider = ExternalGeocentricNeProvider3D(valid)
        periodic_a = provider.get_Ne_geocentric(0.0, 181.0, 150.0)
        periodic_b = provider.get_Ne_geocentric(0.0, -179.0, 150.0)
        metadata = provider.get_metadata()

        invalid_files = (
            {"Ne": Ne_cm3},
            {
                "Ne": Ne_cm3, "latitude_deg": lat, "longitude_deg": lon,
                "height_km": height, "units": "m^-3",
                "coordinate_system": COORDINATE_SYSTEM,
                "axis_order": "longitude,latitude,height",
                "earth_radius_km": R_KM, "longitude_periodic_complete": False,
            },
            {
                "Ne": Ne_cm3, "latitude_deg": lat, "longitude_deg": lon,
                "height_km": height, "units": "unknown",
                "coordinate_system": COORDINATE_SYSTEM,
                "axis_order": AXIS_ORDER, "earth_radius_km": R_KM,
                "longitude_periodic_complete": False,
            },
            {
                "Ne": Ne_cm3, "latitude_deg": lat, "longitude_deg": lon,
                "height_km": height, "units": "m^-3",
                "coordinate_system": "local_cartesian",
                "axis_order": AXIS_ORDER, "earth_radius_km": R_KM,
                "longitude_periodic_complete": False,
            },
            {
                "Ne": Ne_cm3, "latitude_deg": lat, "longitude_deg": lon,
                "height_km": height, "units": "m^-3",
                "coordinate_system": COORDINATE_SYSTEM,
                "axis_order": AXIS_ORDER, "earth_radius_km": R_KM,
                "longitude_periodic_complete": "False",
            },
        )
        for index, fields in enumerate(invalid_files):
            path = directory / f"_c0r2_contract_invalid_{index}.npz"
            created.append(path)
            np.savez(path, **fields)
            try:
                ExternalGeocentricNeProvider3D(path)
            except ValueError:
                rejected += 1
    finally:
        for path in created:
            path.unlink(missing_ok=True)
    metrics = {
        "external_unit_value_m3": periodic_a,
        "external_periodic_difference_m3": abs(periodic_a - periodic_b),
        "invalid_external_contracts_rejected": rejected,
    }
    checks = {
        "external_cm3_conversion": abs(periodic_a - 1.001205e11) < 1.0,
        "external_longitude_periodic": periodic_a == periodic_b,
        "external_hashes_recorded": (
            len(metadata["grid_hash"]) == 64
            and len(metadata["source_metadata"]["source_sha256"]) == 64),
        "external_semantics_recorded": (
            metadata["coordinate_system"] == COORDINATE_SYSTEM
            and metadata["axis_order"] == AXIS_ORDER
            and metadata["longitude_periodic_complete"] is False
            and metadata["source_metadata"]["unit_conversion"]
            == "cm^-3_to_m^-3_x1e6"),
        "ambiguous_external_contracts_rejected": rejected == 5,
    }
    return metrics, checks


def _parameterized_checks() -> tuple[dict, dict]:
    background = _constant_geocentric_background()
    x = np.array([-100.0, 0.0, 100.0, 200.0])
    z = np.array([100.0, 260.0, 380.0])
    local_background = GridNeProvider(np.full((x.size, z.size), 1.0e11), x, z)
    tid = {
        "enable": True, "n_components": 1, "amplitude_list": [0.08],
        "az_deg_list": [0.0], "period_s_list": [1200.0],
        "lambda_h_km_list": [200.0], "z0_km": 260.0,
        "sigma_z_km": 60.0,
    }
    bubble = {
        "enable": True, "delta_max": 0.3, "x0_km": 100.0,
        "y0_km": 0.0, "z0_km": 260.0, "sigma_perp_km": 40.0,
        "sigma_para_km": 500.0,
    }
    spreadf = {
        "enable": True, "delta_Ne_rel_rms": 0.02, "L0_km": 80.0,
        "l0_m": 1000.0, "h_screen_km": 260.0, "sigma_z_km": 60.0,
        "anisotropy_ratio": 20.0, "n_modes": 16,
    }
    es = {
        **cfg.ES_BG, **cfg.ES, "enable": True, "foEs_MHz": 3.0,
        "h_Es_km": 260.0, "delta_h_m": 50000.0,
        "patch_sigma_x_km": 500.0, "patch_sigma_y_km": 500.0,
        "n_modes": 12,
    }
    geomag = {**cfg.GEOMAG, "link_bearing_deg": 0.0}
    new = GeocentricParameterizedIonosphere3DProvider(
        background, tx_lat_deg=30.0, tx_lon_deg=120.0,
        link_bearing_deg=0.0, earth_radius_km=R_KM,
        tid_params=tid, bubble_params=bubble, spreadf_params=spreadf,
        es_params=es, geomag=geomag, seed=17)
    legacy = ParameterizedIonosphere3DProvider(
        local_background, tid_params=tid, bubble_params=bubble,
        spreadf_params=spreadf, es_params=es, geomag=geomag, seed=17)
    target_lat, target_lon = spherical_destination(
        30.0, 120.0, 0.0, 100.0, R_KM)
    local = new.derive_link_coordinates(target_lat, target_lon, 260.0)
    new_value = new.get_Ne_geocentric(target_lat, target_lon, 260.0)
    legacy_value = legacy.get_Ne(local)
    cross_lat, cross_lon = spherical_destination(
        30.0, 120.0, 90.0, 50.0, R_KM)
    cross = new.derive_link_coordinates(cross_lat, cross_lon, 260.0)
    relative_difference = abs(new_value - legacy_value) / max(legacy_value, 1.0)
    metadata = new.get_metadata()
    metrics = {
        "derived_along_link_error_km": float(np.max(np.abs(
            local - np.array([100.0, 0.0, 260.0])))),
        "derived_cross_link_error_km": float(np.max(np.abs(
            cross - np.array([0.0, 50.0, 260.0])))),
        "parameterized_equation_relative_difference": relative_difference,
        "active_component_count": len(metadata["components"]),
    }
    checks = {
        "sphere_derived_along_link": metrics["derived_along_link_error_km"] <= 1.0e-9,
        "sphere_derived_cross_link": metrics["derived_cross_link_error_km"] <= 1.0e-9,
        "component_equations_unchanged": relative_difference <= 1.0e-12,
        "all_parameterized_components_active": len(metadata["components"]) == 4,
        "derived_coordinate_semantics": (
            metadata["component_coordinate_semantics"]
            == "derived_along_link_cross_link_physical_height_km"),
    }
    return metrics, checks


def run() -> int:
    vacuum = VacuumGeocentricNeProvider3D(R_KM)
    chapman = ChapmanGeocentricNeProvider3D(R_KM)
    peak = chapman.get_Ne_geocentric(30.0, 120.0, 350.0)
    peak_ecef = geocentric_to_ecef(30.0, 120.0, 350.0, R_KM)
    gradient = vacuum.get_grad_Ne_ecef(
        geocentric_to_ecef(30.0, 120.0, 300.0, R_KM), R_KM, step_km=0.5)

    lat = np.array([-1.0, 0.0, 1.0])
    lon = np.array([179.0, 180.0, 181.0])
    height = np.array([100.0, 200.0])
    ii, jj, kk = np.meshgrid(np.arange(3), np.arange(3), np.arange(2),
                             indexing="ij")
    Ne = 1.0e11 + ii * 1.0e8 + jj * 1.0e7 + kk * 1.0e6
    grid = GeocentricGridNeProvider3D(Ne, lat, lon, height, R_KM)
    point_value = grid.get_Ne_geocentric(0.0, -179.0, 150.0)
    ecef_value = grid.get_Ne_ecef(
        geocentric_to_ecef(0.0, -179.0, 150.0, R_KM), R_KM)
    global_lon = np.array([-180.0, -90.0, 0.0, 90.0])
    global_Ne = np.full((lat.size, global_lon.size, height.size), 1.0e11)
    global_grid = GeocentricGridNeProvider3D(
        global_Ne, lat, global_lon, height, R_KM,
        longitude_periodic_complete=True)
    global_wrap_a = global_grid.get_Ne_geocentric(0.0, 179.0, 150.0)
    global_wrap_b = global_grid.get_Ne_geocentric(0.0, -181.0, 150.0)
    shape_rejected = False
    try:
        GeocentricGridNeProvider3D(
            Ne[:2], lat, lon, height, R_KM)
    except ValueError:
        shape_rejected = True
    negative_density_rejected = False
    try:
        invalid_Ne = Ne.copy()
        invalid_Ne[0, 0, 0] = -1.0
        GeocentricGridNeProvider3D(invalid_Ne, lat, lon, height, R_KM)
    except ValueError:
        negative_density_rejected = True
    radius_mismatch_rejected = False
    try:
        vacuum.get_Ne_ecef(
            geocentric_to_ecef(30.0, 120.0, 300.0, R_KM), R_KM + 1.0)
    except ValueError:
        radius_mismatch_rejected = True

    external_metrics, external_checks = _external_checks()
    parameterized_metrics, parameterized_checks = _parameterized_checks()
    metrics = {
        "vacuum_gradient_max_m3_per_km": float(np.max(np.abs(gradient))),
        "chapman_peak_relative_error": abs(peak - chapman.peak_Ne_m3)
        / chapman.peak_Ne_m3,
        "chapman_ecef_relative_error": abs(
            chapman.get_Ne_ecef(peak_ecef, R_KM) - peak) / peak,
        "grid_geocentric_ecef_difference_m3": abs(point_value - ecef_value),
        **external_metrics,
        **parameterized_metrics,
    }
    checks = {
        "formal_interface_has_no_ambiguous_get_Ne": not hasattr(vacuum, "get_Ne"),
        "vacuum_provider": vacuum.get_Ne_geocentric(30.0, 120.0, 300.0) == 0.0,
        "vacuum_ecef_gradient": metrics["vacuum_gradient_max_m3_per_km"] == 0.0,
        "chapman_peak": metrics["chapman_peak_relative_error"] <= 1.0e-14,
        "chapman_same_geocentric_ecef_point": (
            metrics["chapman_ecef_relative_error"] <= 1.0e-14),
        "grid_same_geocentric_ecef_point": (
            metrics["grid_geocentric_ecef_difference_m3"] <= 1.0),
        "global_longitude_boundary_interpolated": (
            global_wrap_a == 1.0e11 and global_wrap_a == global_wrap_b),
        "grid_axis_swap_rejected": shape_rejected,
        "negative_density_rejected": negative_density_rejected,
        "earth_radius_mismatch_rejected": radius_mismatch_rejected,
        "grid_metadata_complete": (
            grid.get_metadata()["axis_order"] == AXIS_ORDER
            and len(grid.get_metadata()["grid_hash"]) == 64),
        **external_checks,
        **parameterized_checks,
    }
    report = {
        "schema_version": "a4-c0r2-provider-v1",
        "stage": "A4-C0R-2",
        "status": "passed" if all(checks.values()) else "failed",
        "case_semantics": "geocentric_ecef_ne_provider_contract_only",
        "coordinate_system": COORDINATE_SYSTEM,
        "axis_order": AXIS_ORDER,
        "internal_Ne_units": "m^-3",
        "earth_radius_km": R_KM,
        "metrics": metrics,
        "checks": checks,
        "source_sha256": {
            "provider": _sha256(ROOT / "models" / "ne_provider_geocentric.py"),
            "geometry": _sha256(ROOT / "utils.py"),
            "test": _sha256(Path(__file__)),
        },
        "claim_boundary": (
            "Provider coordinates, interpolation, units, hashes, and unchanged "
            "parameterized component equations only; AH3D/P2P/PE not migrated"
        ),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(REPORT)
    for name, passed in checks.items():
        print(f" {name}: {'OK' if passed else 'FAIL'}")
    failures = sum(not passed for passed in checks.values())
    print(f"A4-C0R-2: {failures} fail(s)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run())
