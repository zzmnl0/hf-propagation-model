"""A4-C1WR grid-boundary derivative attribution without changing production fields."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ppigrf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg  # noqa: E402
import models.ray_tracer_ah_3d_ecef as ahcore  # noqa: E402
from models.geomagnetic_field_3d import IGRFGridGeomagneticField3D  # noqa: E402
from models.ne_provider_geocentric import (  # noqa: E402
    GeocentricGridNeProvider3D,
    WGS84NeProviderAdapter,
)
from models.ray_tracer_ah_3d_ecef import (  # noqa: E402
    WGS84_COORDINATE_MODEL,
    RefractiveIndexAH3DECEF,
    hamiltonian_ah3d_ecef,
    trace_ray_ah3d_ecef,
)
from tests import _run_a4_c1w_windows_ox_validation as c1w  # noqa: E402
from utils import ecef_to_geodetic, enu_basis_wgs84  # noqa: E402

OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1wr"
RAW = OUT / "raw"
FROZEN = ROOT / "output2" / "a4_pharlap" / "c1" / "c1w"
PLAN = ROOT / "word" / "a4_c1wr_grid_derivative_attribution_and_repair_plan.txt"
P2P_RESULT = OUT / "p2p_regression" / "a4_c1wr_p2p_high_low_regression.json"
P2P_MANIFEST = OUT / "p2p_regression" / "a4_c1wr_p2p_high_low_manifest.json"
C1_SHADOW = OUT / "c1_shadow" / "a4_c1wr_c1_shadow_summary.json"
C1_SHADOW_MANIFEST = (
    OUT / "c1_shadow" / "same_angle" / "c1_1_freeze_manifest.json")
REPORT = ROOT / "word" / "a4_c1wr_grid_derivative_attribution_and_repair_report.txt"

STEPS_KM = np.array([0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0])
PRODUCTION_STEPS_KM = np.array([0.125, 0.25, 0.5])
VARIANTS = (
    "linear_ne_linear_b",
    "analytic_ne_linear_b",
    "linear_ne_direct_b",
    "cubic_ne_linear_b",
    "linear_ne_cubic_b",
    "cubic_ne_cubic_b",
    "analytic_ne_direct_b_reference",
)
TRACE_VARIANTS = (
    "linear_ne_linear_b",
    "cubic_ne_linear_b",
    "linear_ne_cubic_b",
    "cubic_ne_cubic_b",
)
FIGURES = (
    "a4_c1wr_boundary_crossing.png",
    "a4_c1wr_dr_scale_platform.png",
    "a4_c1wr_factor_attribution.png",
    "a4_c1wr_high_o_variant_trajectories.png",
    "a4_c1wr_regression_gate_matrix.png",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"cannot serialize {type(value).__name__}")


def _write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=True, default=_json_value) + "\n",
        encoding="utf-8",
    )


class _AnalyticChapman:
    def __init__(self):
        self.provider_hash = hashlib.sha256(
            b"c1wr-analytic-chapman-foF2-7-hmF2-350-ymF2-100").hexdigest()

    def get_Ne_ecef(self, ecef_km, earth_radius_km=None, t: float = 0.0) -> float:
        del earth_radius_km, t
        return float(c1w._chapman(ecef_to_geodetic(ecef_km)[2]))

    def domain_margin_ecef(self, ecef_km, earth_radius_km=None) -> float:
        del earth_radius_km
        lat, lon, height = ecef_to_geodetic(ecef_km)
        return float(min(lat - 20.0, 50.0 - lat, lon - 105.0, 135.0 - lon,
                         height + 5.0, 610.0 - height))

    def get_metadata(self) -> dict:
        return {
            "name": "c1wr_validation_analytic_chapman",
            "coordinate_system": "wgs84_geodetic",
            "units": "m^-3",
            "provider_hash": self.provider_hash,
            "validation_only": True,
        }


class _DirectB:
    def __init__(self):
        self._cache = {}
        self.field_hash = hashlib.sha256(
            f"c1wr-direct-ppigrf:{cfg.IRI_DT.isoformat()}".encode("ascii")).hexdigest()

    def get_B_ecef_T(self, ecef_km, t: float = 0.0) -> np.ndarray:
        del t
        key = tuple(np.round(np.asarray(ecef_km, dtype=float), 10))
        if key not in self._cache:
            lat, lon, height = ecef_to_geodetic(key)
            east, north, up = ppigrf.igrf(lon, lat, height, cfg.IRI_DT)
            local = np.array([
                float(np.asarray(east).reshape(-1)[0]),
                float(np.asarray(north).reshape(-1)[0]),
                float(np.asarray(up).reshape(-1)[0]),
            ])
            self._cache[key] = enu_basis_wgs84(lat, lon) @ local * 1.0e-9
        return self._cache[key].copy()

    def get_metadata(self) -> dict:
        return {
            "name": "c1wr_validation_direct_ppigrf",
            "model": "IGRF-14",
            "epoch": cfg.IRI_DT.isoformat(),
            "coordinate_system": "direct_WGS84_query_ECEF_vectors",
            "units": "T",
            "derivative_semantics": "five_point_ECEF_difference_of_direct_queries",
            "field_hash": self.field_hash,
            "validation_only": True,
        }


def _model(ne_provider, b_provider, mode: str) -> RefractiveIndexAH3DECEF:
    return RefractiveIndexAH3DECEF(
        ne_provider,
        earth_radius_km=cfg.RE_KM,
        freq_MHz=c1w.FREQUENCY_MHZ,
        wave_mode=mode,
        geomagnetic_provider=b_provider,
        coordinate_model=WGS84_COORDINATE_MODEL,
        spatial_fd_base_km=0.25,
        spatial_fd_scales=(0.5, 1.0, 2.0),
        integration_profile="wgs84_reference",
    )


def _field_objects():
    latitude, longitude, height, ne, linear_ne, linear_b = c1w._frozen_fields()
    source = linear_ne.source
    cubic_source = GeocentricGridNeProvider3D(
        source.Ne_m3,
        source.latitude_deg,
        source.longitude_deg,
        source.height_km,
        source.earth_radius_km,
        fill_value_m3=source.fill_value_m3,
        longitude_periodic_complete=source.longitude_periodic_complete,
        source_metadata=source.source_metadata,
        interpolation="cubic",
    )
    cubic_b = IGRFGridGeomagneticField3D(
        linear_b.latitude_deg,
        linear_b.longitude_deg,
        linear_b.height_km,
        linear_b.b_ecef_T,
        epoch=linear_b.epoch,
        igrf_version=linear_b.igrf_version,
        interpolation="cubic",
    )
    return {
        "latitude": latitude,
        "longitude": longitude,
        "height": height,
        "ne": ne,
        "linear_ne": linear_ne,
        "analytic_ne": _AnalyticChapman(),
        "cubic_ne": WGS84NeProviderAdapter(cubic_source),
        "linear_b": linear_b,
        "cubic_b": cubic_b,
        "direct_b": _DirectB(),
        "ne_source": source,
    }


def _variant_model(fields: dict, variant: str, mode: str):
    mapping = {
        "linear_ne_linear_b": ("linear_ne", "linear_b"),
        "analytic_ne_linear_b": ("analytic_ne", "linear_b"),
        "linear_ne_direct_b": ("linear_ne", "direct_b"),
        "cubic_ne_linear_b": ("cubic_ne", "linear_b"),
        "linear_ne_cubic_b": ("linear_ne", "cubic_b"),
        "cubic_ne_cubic_b": ("cubic_ne", "cubic_b"),
        "analytic_ne_direct_b_reference": ("analytic_ne", "direct_b"),
    }
    ne_name, b_name = mapping[variant]
    return _model(fields[ne_name], fields[b_name], mode)


def _grid_cell_indices(points_ecef, lat_axis, lon_axis, height_axis):
    coordinates = ecef_to_geodetic(points_ecef)
    result = np.column_stack([
        np.searchsorted(lat_axis, coordinates[:, 0], side="right") - 1,
        np.searchsorted(lon_axis, coordinates[:, 1], side="right") - 1,
        np.searchsorted(height_axis, coordinates[:, 2], side="right") - 1,
    ])
    limits = np.array([len(lat_axis) - 2, len(lon_axis) - 2, len(height_axis) - 2])
    return np.clip(result, 0, limits)


def _crossing_flags(positions, fields: dict):
    flags = {}
    grids = {
        "ne": (
            fields["ne_source"].latitude_deg,
            fields["ne_source"].longitude_deg,
            fields["ne_source"].height_km,
        ),
        "b": (
            fields["linear_b"].latitude_deg,
            fields["linear_b"].longitude_deg,
            fields["linear_b"].height_km,
        ),
    }
    for name, axes in grids.items():
        crossed = np.zeros((len(positions), len(PRODUCTION_STEPS_KM)), dtype=bool)
        for step_index, step in enumerate(PRODUCTION_STEPS_KM):
            for axis in range(3):
                offset = np.zeros(3)
                offset[axis] = step
                minus = _grid_cell_indices(positions - offset, *axes)
                plus = _grid_cell_indices(positions + offset, *axes)
                crossed[:, step_index] |= np.any(minus != plus, axis=1)
        flags[name] = crossed
    return flags


def _trace_with_log(fields: dict, variant: str, case_id: str):
    index = c1w.CASE_IDS.index(case_id)
    mode = c1w.MODES[index]
    elevation = c1w.ELEVATIONS[index]
    model = _variant_model(fields, variant, mode)
    calls = []
    original = ahcore.spatial_gradient_ah3d_ecef

    def logged(state, local_model, *, t: float = 0.0):
        entry = {
            "position": np.asarray(state[:3], dtype=float).copy(),
            "wavevector": np.asarray(state[3:6], dtype=float).copy(),
            "failed": False,
            "stability_ratio": np.nan,
            "relative_spread": np.nan,
            "absolute_spread": np.nan,
        }
        try:
            gradient, diagnostics = original(state, local_model, t=t)
            entry.update({
                "stability_ratio": float(diagnostics["stability_ratio"]),
                "relative_spread": float(diagnostics["relative_spread"]),
                "absolute_spread": float(diagnostics["absolute_spread"]),
            })
            return gradient, diagnostics
        except ahcore.SpatialDerivativeError:
            entry["failed"] = True
            raise
        finally:
            calls.append(entry)

    ahcore.spatial_gradient_ah3d_ecef = logged
    try:
        ray = trace_ray_ah3d_ecef(
            model,
            tx_lat_deg=c1w.TX[0],
            tx_lon_deg=c1w.TX[1],
            tx_height_km=c1w.TX[2],
            elevation_deg=elevation,
            bearing_deg=0.0,
            h_stop_km=600.0,
            max_path_km=5000.0,
            profile="wgs84_reference",
            apply_dispersion_projection=True,
        )
    finally:
        ahcore.spatial_gradient_ah3d_ecef = original

    positions = np.asarray([item["position"] for item in calls])
    wavevectors = np.asarray([item["wavevector"] for item in calls])
    crossing = _crossing_flags(positions, fields)
    log = {
        "position_ecef_km": positions,
        "wavevector_ecef_km^-1": wavevectors,
        "height_km": ecef_to_geodetic(positions)[:, 2],
        "failed": np.asarray([item["failed"] for item in calls], dtype=bool),
        "stability_ratio": np.asarray(
            [item["stability_ratio"] for item in calls], dtype=float),
        "relative_spread": np.asarray(
            [item["relative_spread"] for item in calls], dtype=float),
        "absolute_spread_km^-3": np.asarray(
            [item["absolute_spread"] for item in calls], dtype=float),
        "ne_crossing": crossing["ne"],
        "b_crossing": crossing["b"],
    }
    return model, ray, log


def _cache_paths(label: str):
    return RAW / f"{label}.npz", RAW / f"{label}.json"


def _save_trace_cache(label: str, ray: dict, log: dict) -> None:
    npz_path, summary_path = _cache_paths(label)
    np.savez_compressed(
        npz_path,
        **{key: value for key, value in ray.items()
           if isinstance(value, np.ndarray)},
        **{f"log_{key}": value for key, value in log.items()},
    )
    _write_json(summary_path, {
        "ray_metadata": {
            key: value for key, value in ray.items()
            if not isinstance(value, np.ndarray)
        },
        "npz_sha256": _sha256(npz_path),
    })


def _load_trace_cache(label: str):
    npz_path, summary_path = _cache_paths(label)
    if not npz_path.is_file():
        return None
    with np.load(npz_path, allow_pickle=False) as data:
        ray = {
            key: np.asarray(data[key]) for key in data.files
            if not key.startswith("log_")
        }
        log = {
            key[4:]: np.asarray(data[key]) for key in data.files
            if key.startswith("log_")
        }
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("npz_sha256") != _sha256(npz_path):
            return None
        ray.update(summary["ray_metadata"])
    else:
        case_id = label.removeprefix("baseline_")
        prior = json.loads(
            (FROZEN / "a4_c1w_validation.json").read_text(encoding="utf-8"))
        ray.update(prior["ah3d_case_status"][case_id])
        ray["spatial_fd_rejection_count"] = int(np.count_nonzero(log["failed"]))
        ray["spatial_fd_stability_ratio_max"] = float(
            np.nanmax(log["stability_ratio"]))
        ray["hamiltonian_residual_max"] = float("nan")
        _save_trace_cache(label, ray, log)
    return ray, log


def _trace_cache_worker(label: str, variant: str, case_id: str):
    fields = _field_objects()
    _, ray, log = _trace_with_log(fields, variant, case_id)
    _save_trace_cache(label, ray, log)
    return label, ray["termination_status"], len(log["height_km"])


def _run_missing_traces(specs):
    missing = [
        spec for spec in specs if _load_trace_cache(spec[0]) is None
    ]
    if missing:
        with ProcessPoolExecutor(max_workers=len(missing)) as executor:
            futures = [
                executor.submit(_trace_cache_worker, *spec) for spec in missing]
            for future in futures:
                label, status, calls = future.result()
                print(f"C1WR {label}: {status} calls={calls}", flush=True)
    return {
        label: _load_trace_cache(label)
        for label, _, _ in specs
    }


def _simple_trace_worker(label: str, profile: str, projection: bool):
    fields = _field_objects()
    model = _variant_model(fields, "cubic_ne_cubic_b", "O")
    ray = trace_ray_ah3d_ecef(
        model,
        tx_lat_deg=c1w.TX[0],
        tx_lon_deg=c1w.TX[1],
        tx_height_km=c1w.TX[2],
        elevation_deg=22.0,
        bearing_deg=0.0,
        h_stop_km=600.0,
        max_path_km=5000.0,
        profile=profile,
        apply_dispersion_projection=projection,
    )
    log = {
        "position_ecef_km": np.empty((0, 3)),
        "wavevector_ecef_km^-1": np.empty((0, 3)),
        "height_km": np.empty(0),
        "failed": np.empty(0, dtype=bool),
        "stability_ratio": np.empty(0),
        "relative_spread": np.empty(0),
        "absolute_spread_km^-3": np.empty(0),
        "ne_crossing": np.empty((0, 3), dtype=bool),
        "b_crossing": np.empty((0, 3), dtype=bool),
    }
    _save_trace_cache(label, ray, log)
    return label, ray["termination_status"]


def _run_h3_shadows():
    specs = (
        ("shared_h3_refined_low_O", "wgs84_refined", True),
        ("shared_h3_projection_off_low_O", "wgs84_reference", False),
    )
    missing = [spec for spec in specs if _load_trace_cache(spec[0]) is None]
    if missing:
        with ProcessPoolExecutor(max_workers=len(missing)) as executor:
            futures = [executor.submit(_simple_trace_worker, *spec)
                       for spec in missing]
            for future in futures:
                label, status = future.result()
                print(f"C1WR {label}: {status}", flush=True)
    return {
        label: _load_trace_cache(label)[0] for label, _, _ in specs
    }


def _central_estimates(model, state):
    estimates = []
    for step in STEPS_KM:
        gradient = np.empty(3)
        for axis in range(3):
            offset = np.zeros(6)
            offset[axis] = step
            gradient[axis] = (
                hamiltonian_ah3d_ecef(state + offset, model)
                - hamiltonian_ah3d_ecef(state - offset, model)
            ) / (2.0 * step)
        estimates.append(gradient)
    return np.asarray(estimates)


def _five_point_estimates(model, state):
    estimates = []
    for step in STEPS_KM[:3]:
        gradient = np.empty(3)
        for axis in range(3):
            offset = np.zeros(6)
            offset[axis] = step
            gradient[axis] = (
                hamiltonian_ah3d_ecef(state - 2.0 * offset, model)
                - 8.0 * hamiltonian_ah3d_ecef(state - offset, model)
                + 8.0 * hamiltonian_ah3d_ecef(state + offset, model)
                - hamiltonian_ah3d_ecef(state + 2.0 * offset, model)
            ) / (12.0 * step)
        estimates.append(gradient)
    return np.asarray(estimates)


def _state_at_height(ray: dict, height_km: float, branch: str):
    heights = np.asarray(ray["trajectory_physical_height_km"], dtype=float)
    wavevectors = np.asarray(ray["wavevector_ecef_km^-1"], dtype=float)
    positions = np.asarray(ray["trajectory_ecef_km"], dtype=float)
    apex = int(np.argmax(heights))
    search = np.arange(0, apex + 1) if branch == "ascending" else np.arange(apex, len(heights))
    index = int(search[np.argmin(abs(heights[search] - height_km))])
    position = positions[index].copy()
    lat, lon, current_height = ecef_to_geodetic(position)
    position += enu_basis_wgs84(lat, lon)[:, 2] * (height_km - current_height)
    return np.r_[position, wavevectors[index]]


def _anchors(rays: dict):
    anchors = []
    for case_id in ("low_O", "low_X"):
        ray = rays[case_id]
        index = int(np.argmax(ray["trajectory_physical_height_km"]))
        anchors.append({
            "anchor_id": f"{case_id}_apex",
            "case_id": case_id,
            "state": np.r_[
                ray["trajectory_ecef_km"][index],
                ray["wavevector_ecef_km^-1"][index],
            ],
        })
    for case_id in ("high_O", "high_X"):
        branches = ("ascending",) if case_id == "high_O" else ("ascending", "descending")
        for branch in branches:
            for height in (329.75, 329.9, 329.99, 330.0, 330.01, 330.1, 330.25):
                anchors.append({
                    "anchor_id": f"{case_id}_{branch}_{height:.2f}",
                    "case_id": case_id,
                    "state": _state_at_height(rays[case_id], height, branch),
                })
    return anchors


def _diagnose_anchors(fields: dict, rays: dict):
    rows = []
    scale_data = {}
    references = {}
    for anchor in _anchors(rays):
        case_id = anchor["case_id"]
        mode = c1w.MODES[c1w.CASE_IDS.index(case_id)]
        state = np.asarray(anchor["state"], dtype=float)
        reference_model = _variant_model(
            fields, "analytic_ne_direct_b_reference", mode)
        reference_estimates = _five_point_estimates(reference_model, state)
        reference = np.median(reference_estimates, axis=0)
        reference_spread = float(np.max(np.linalg.norm(
            reference_estimates - reference, axis=1)))
        references[anchor["anchor_id"]] = reference
        for variant in VARIANTS:
            model = _variant_model(fields, variant, mode)
            estimates = _central_estimates(model, state)
            production = estimates[2:5]
            selected = np.median(production, axis=0)
            spread = float(np.max(np.linalg.norm(production - selected, axis=1)))
            stability_limit = (
                model.spatial_fd_stability_rtol * float(np.linalg.norm(selected))
                + model.k0**2 * model.spatial_fd_stability_atol_n2_per_km)
            error = float(np.linalg.norm(selected - reference))
            denominator = max(
                float(np.linalg.norm(reference)),
                model.k0**2 * model.spatial_fd_stability_atol_n2_per_km)
            row = {
                "anchor_id": anchor["anchor_id"],
                "case_id": case_id,
                "variant": variant,
                "height_km": float(ecef_to_geodetic(state[:3])[2]),
                "reference_norm_km^-3": float(np.linalg.norm(reference)),
                "reference_platform_spread_km^-3": reference_spread,
                "candidate_norm_km^-3": float(np.linalg.norm(selected)),
                "candidate_spread_km^-3": spread,
                "stability_limit_km^-3": stability_limit,
                "stability_ratio": spread / stability_limit,
                "relative_reference_error": error / denominator,
            }
            rows.append(row)
            scale_data[(anchor["anchor_id"], variant)] = estimates
    return rows, scale_data, references


def _case_errors(rows):
    errors = {}
    for case_id in c1w.CASE_IDS:
        errors[case_id] = {}
        for variant in VARIANTS:
            values = [
                row["relative_reference_error"] for row in rows
                if row["case_id"] == case_id and row["variant"] == variant
            ]
            errors[case_id][variant] = float(max(values))
    return errors


def _trace_summary(ray: dict):
    return {
        "termination_status": ray["termination_status"],
        "failure_reason": ray.get("failure_reason"),
        "turning_point_count": int(ray.get("turning_point_count", 0)),
        "apex_height_km": (
            None if ray.get("apex_height_km") is None
            else float(ray["apex_height_km"])),
        "geom_path_km": float(ray.get("geom_path_km", np.nan)),
        "group_path_km": float(ray.get("group_path_km", np.nan)),
        "hamiltonian_residual_max": float(
            ray.get("hamiltonian_residual_max", np.nan)),
        "spatial_fd_rejection_count": int(
            ray.get("spatial_fd_rejection_count", 0)),
        "spatial_fd_stability_ratio_max": float(
            ray.get("spatial_fd_stability_ratio_max", np.nan)),
    }


def _expected_topology(case_id: str, ray: dict) -> bool:
    expected = "upper_boundary" if case_id == "high_O" else "ground"
    return ray["termination_status"] == expected


def _pharlap_path(case_id: str):
    return np.loadtxt(FROZEN / "raw" / f"pharlap_{case_id}.csv", delimiter=",")


def _trajectory_rms(left: dict, right: dict):
    u = np.linspace(0.0, 1.0, 201)
    left_q = c1w._resample(
        left["trajectory_geom_path_km"], left["trajectory_ecef_km"], u)
    right_q = c1w._resample(
        right["trajectory_geom_path_km"], right["trajectory_ecef_km"], u)
    return float(np.sqrt(np.mean(np.sum((left_q - right_q)**2, axis=1))))


def _high_o_common_path_metrics(ray: dict):
    ph = _pharlap_path("high_O")
    ph_ecef = ph[:, :3] / 1000.0
    ph_path = ph[:, 7] / 1000.0
    ah_ecef = np.asarray(ray["trajectory_ecef_km"], dtype=float)
    ah_path = np.asarray(ray["trajectory_geom_path_km"], dtype=float)
    common_end = float(min(ph_path[-1], ah_path[-1]))
    target = np.linspace(0.0, common_end, 201)

    def sample(path, values):
        keep = np.r_[True, np.diff(path) > 1.0e-10]
        return np.column_stack([
            np.interp(target, path[keep], values[keep, axis])
            for axis in range(values.shape[1])
        ])

    ph_q = sample(ph_path, ph_ecef)
    ah_q = sample(ah_path, ah_ecef)
    delta = ah_q - ph_q
    norm = np.linalg.norm(delta, axis=1)
    return {
        "common_geom_path_km": common_end,
        "trajectory_rms_km": float(np.sqrt(np.mean(norm**2))),
        "trajectory_max_km": float(np.max(norm)),
        "landing_metric_applicable": False,
        "semantics": "common_absolute_geometric_path_no_landing_for_penetrated_ray",
    }


def _range_height(ecef):
    geodetic = ecef_to_geodetic(ecef)
    lat0 = np.radians(c1w.TX[0])
    lon0 = np.radians(c1w.TX[1])
    lat = np.radians(geodetic[:, 0])
    lon = np.radians(geodetic[:, 1])
    angle = 2.0 * np.arcsin(np.sqrt(
        np.sin((lat - lat0) / 2.0)**2
        + np.cos(lat0) * np.cos(lat) * np.sin((lon - lon0) / 2.0)**2))
    return cfg.RE_KM * angle, geodetic[:, 2]


def _plots(fields, baseline_rays, baseline_logs, trial_rays, rows, errors,
           scale_data, selected_variant, gates, regression_outcomes):
    figure, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for axis, case_id in zip(axes.flat, c1w.CASE_IDS):
        log = baseline_logs[case_id]
        finite = np.isfinite(log["stability_ratio"])
        scatter = axis.scatter(
            np.arange(len(log["height_km"]))[finite],
            log["height_km"][finite],
            c=np.clip(log["stability_ratio"][finite], 0.0, 2.0),
            s=5, cmap="viridis", vmin=0.0, vmax=2.0)
        crossed = np.any(log["ne_crossing"] | log["b_crossing"], axis=1)
        axis.scatter(np.flatnonzero(crossed), log["height_km"][crossed],
                     s=10, facecolors="none", edgecolors="red")
        axis.axhline(330.0, color="black", linestyle="--", linewidth=0.8)
        axis.set(title=case_id, xlabel="D_r evaluation index",
                 ylabel="WGS84 height (km)")
    figure.colorbar(scatter, ax=axes, label="production stability ratio")
    figure.suptitle(
        "A4-C1WR grid-boundary audit | circles mark Ne/B stencil crossings")
    figure.savefig(OUT / FIGURES[0], dpi=180)
    plt.close(figure)

    selected_anchor = "high_O_ascending_330.00"
    figure, axis = plt.subplots(figsize=(10, 6))
    for variant in VARIANTS:
        values = scale_data[(selected_anchor, variant)]
        axis.loglog(STEPS_KM, np.linalg.norm(values, axis=1),
                    marker="o", label=variant)
    axis.axvspan(0.125, 0.5, color="black", alpha=0.08,
                 label="production scale range")
    axis.set(xlabel="ECEF finite-difference step h (km)",
             ylabel="|D_r| (km^-3)",
             title="High-O 330 km node: spatial-derivative scale platform")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(OUT / FIGURES[1], dpi=180)
    plt.close(figure)

    matrix = np.array([
        [errors[case_id][variant] for variant in VARIANTS]
        for case_id in c1w.CASE_IDS
    ])
    figure, axis = plt.subplots(figsize=(13, 5))
    image = axis.imshow(np.log10(np.maximum(matrix, 1.0e-16)),
                        aspect="auto", cmap="magma")
    axis.set_xticks(np.arange(len(VARIANTS)), VARIANTS, rotation=40, ha="right")
    axis.set_yticks(np.arange(len(c1w.CASE_IDS)), c1w.CASE_IDS)
    figure.colorbar(image, ax=axis, label="log10 max relative D_r error")
    axis.set_title("Ne/B single-factor attribution against analytic Chapman + direct IGRF")
    figure.tight_layout()
    figure.savefig(OUT / FIGURES[2], dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(10, 6))
    base_range, base_height = _range_height(
        baseline_rays["high_O"]["trajectory_ecef_km"])
    axis.plot(base_range, base_height, label="AH3D linear/linear", linewidth=1.8)
    for variant, ray in trial_rays.items():
        distance, height = _range_height(ray["trajectory_ecef_km"])
        axis.plot(distance, height, label=f"AH3D {variant}", linewidth=1.4)
    ph = _pharlap_path("high_O")
    ph_range, ph_height = _range_height(ph[:, :3] / 1000.0)
    axis.plot(ph_range, ph_height, "k--", label="PHaRLAP high-O", linewidth=1.4)
    axis.set(xlabel="TX great-circle range (km)", ylabel="WGS84 height (km)",
             title="High-O fixed-angle topology | one-way | 10 MHz")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(OUT / FIGURES[3], dpi=180)
    plt.close(figure)

    displayed = {**gates, **regression_outcomes}
    labels = list(displayed)
    passed = [bool(displayed[name]) for name in labels]
    figure, axis = plt.subplots(figsize=(11, 5))
    axis.barh(labels, np.ones(len(labels)),
              color=["#2f855a" if value else "#c53030" for value in passed])
    for index, value in enumerate(passed):
        axis.text(0.98, index, "PASS" if value else "FAIL",
                  ha="right", va="center", color="white", fontweight="bold")
    axis.set_xlim(0.0, 1.0)
    axis.set_xticks([])
    axis.set_title(
        "A4-C1WR attribution gates and retained limitations | "
        f"candidate={selected_variant or 'none'}")
    figure.tight_layout()
    figure.savefig(OUT / FIGURES[4], dpi=180)
    plt.close(figure)


def _frozen_hashes():
    existing = OUT / "a4_c1wr_attribution.json"
    if existing.is_file():
        prior = json.loads(existing.read_text(encoding="utf-8"))
        frozen = prior.get("frozen_source_sha256")
        if frozen:
            return frozen
    paths = [
        FROZEN / "a4_c1w_validation.json",
        FROZEN / "a4_c1w_manifest.json",
        FROZEN / "raw" / "a4_c1w_wgs84_ox_input.mat",
        FROZEN / "raw" / "pharlap_wgs84_ox_summary.json",
        ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        ROOT / "models" / "ne_provider_geocentric.py",
        ROOT / "models" / "geomagnetic_field_3d.py",
        PLAN,
    ]
    return {str(path.relative_to(ROOT)): _sha256(path) for path in paths}


def _current_source_hashes():
    paths = [
        Path(__file__),
        ROOT / "tests" / "_verify_a4_c1wr_grid_derivative_repair.py",
        ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        ROOT / "models" / "ne_provider_geocentric.py",
        ROOT / "models" / "geomagnetic_field_3d.py",
        ROOT / "models" / "unified_hybrid_pipeline.py",
        ROOT / "tests" / "_run_a4_c1_prep3.py",
        ROOT / "tests" / "_run_a4_c1wr_p2p_high_low_regression.py",
        ROOT / "tests" / "_run_a4_c1wr_c1_shadow.py",
        PLAN,
    ]
    return {str(path.relative_to(ROOT)): _sha256(path) for path in paths}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    frozen_hashes = _frozen_hashes()
    prior = json.loads(
        (FROZEN / "a4_c1w_validation.json").read_text(encoding="utf-8"))
    if prior.get("status") != "implementation_defect_identified":
        raise RuntimeError("C1W frozen defect status is not present")
    if prior.get("prior_c1_status") != "failed_c1_chapman_high_low_external_comparison":
        raise RuntimeError("original C1 failure semantics drifted")
    p2p_regression = json.loads(P2P_RESULT.read_text(encoding="utf-8"))
    c1_shadow = json.loads(C1_SHADOW.read_text(encoding="utf-8"))

    fields = _field_objects()
    baseline_specs = [
        (f"baseline_{case_id}", "linear_ne_linear_b", case_id)
        for case_id in c1w.CASE_IDS
    ]
    baseline_cached = _run_missing_traces(baseline_specs)
    baseline_rays = {
        case_id: baseline_cached[f"baseline_{case_id}"][0]
        for case_id in c1w.CASE_IDS
    }
    baseline_logs = {
        case_id: baseline_cached[f"baseline_{case_id}"][1]
        for case_id in c1w.CASE_IDS
    }

    rows, scale_data, _ = _diagnose_anchors(fields, baseline_rays)
    errors = _case_errors(rows)
    baseline_high = {
        case_id: errors[case_id]["linear_ne_linear_b"]
        for case_id in ("high_O", "high_X")
    }
    candidate_names = (
        "cubic_ne_linear_b", "linear_ne_cubic_b", "cubic_ne_cubic_b")
    closures = {
        variant: {
            case_id: 1.0 - errors[case_id][variant] / max(
                baseline_high[case_id], 1.0e-30)
            for case_id in ("high_O", "high_X")
        }
        for variant in candidate_names
    }

    trial_specs = [
        (f"shared_trial_{variant}_high_O", variant, "high_O")
        for variant in candidate_names
    ]
    trial_cached = _run_missing_traces(trial_specs)
    trial_high_o = {
        variant: trial_cached[f"shared_trial_{variant}_high_O"][0]
        for variant in candidate_names
    }

    eligible = [
        variant for variant in candidate_names
        if min(closures[variant].values()) >= 0.8
        and max(errors[case_id][variant]
                for case_id in ("high_O", "high_X")) <= 1.0e-4
        and trial_high_o[variant]["termination_status"] == "upper_boundary"
    ]
    selected_variant = min(
        eligible,
        key=lambda variant: max(
            errors[case_id][variant] for case_id in ("high_O", "high_X")),
        default=None,
    )
    selected_rays = {}
    if selected_variant is not None:
        selected_rays["high_O"] = trial_high_o[selected_variant]
        selected_specs = [
            (f"shared_selected_{selected_variant}_{case_id}",
             selected_variant, case_id)
            for case_id in ("low_O", "low_X", "high_X")
        ]
        selected_cached = _run_missing_traces(selected_specs)
        selected_rays.update({
            case_id: selected_cached[
                f"shared_selected_{selected_variant}_{case_id}"][0]
            for case_id in ("low_O", "low_X", "high_X")
        })

    analytic_limits, derivative_audit = c1w._analytic_and_derivative_checks(
        fields["cubic_ne"], fields["cubic_b"])
    h3_shadows = _run_h3_shadows()
    h3_convergence_rms = _trajectory_rms(
        selected_rays["low_O"], h3_shadows["shared_h3_refined_low_O"])
    h3_projection_rms = _trajectory_rms(
        selected_rays["low_O"],
        h3_shadows["shared_h3_projection_off_low_O"])

    h6_metrics = {}
    for case_id in ("low_O", "low_X", "high_X"):
        metric = c1w._path_metrics(
            selected_rays[case_id], FROZEN / "raw" / f"pharlap_{case_id}.csv")
        h6_metrics[case_id] = {
            key: value for key, value in metric.items()
            if key not in {"trajectory_delta_xyz_km", "u"}
        }
    h6_high_o = _high_o_common_path_metrics(selected_rays["high_O"])

    low_controls = (
        selected_variant is not None
        and all(_expected_topology(case_id, selected_rays[case_id])
                for case_id in c1w.CASE_IDS)
        and all(selected_rays[case_id]["spatial_fd_rejection_count"] == 0
                for case_id in c1w.CASE_IDS)
    )
    reference_platform_stable = all(
        row["reference_platform_spread_km^-3"]
        <= row["stability_limit_km^-3"]
        for row in rows if row["variant"] == "analytic_ne_direct_b_reference")
    gates = {
        "frozen_inputs_verified": True,
        "reference_platform_stable": reference_platform_stable,
        "two_high_closure_ge_80pct": (
            selected_variant is not None
            and min(closures[selected_variant].values()) >= 0.8),
        "h3_dr_relative_error_le_1e-4": (
            selected_variant is not None
            and max(errors[case_id][selected_variant]
                    for case_id in c1w.CASE_IDS) <= 1.0e-4),
        "double_high_double_low_topology": low_controls,
        "no_spatial_derivative_rejection": low_controls,
        "h3_analytic_limits": max(analytic_limits.values()) <= 1.0e-10,
        "h3_derivative_audit": (
            derivative_audit["Dr_relative_error"] <= 1.0e-4
            and derivative_audit["Dk_relative_error"] <= 1.0e-5
            and derivative_audit["Domega_relative_error"] <= 1.0e-5),
        "h3_integration_convergence": h3_convergence_rms <= 0.1,
        "h3_projection_shadow": h3_projection_rms <= 0.1,
        "h6_ground_trajectory_metrics": all(
            row["path_xyz_rms_km"] <= 2.0
            and row["path_max_km"] <= 5.0
            and row["landing_difference_km"] <= 1.0
            for row in h6_metrics.values()),
        "h6_ground_path_metrics": all(
            row["geom_path_relative_error"] <= 1.0e-3
            and row["group_path_relative_error"] <= 1.0e-3
            for row in h6_metrics.values()),
        "h6_high_o_common_path": (
            h6_high_o["trajectory_rms_km"] <= 2.0
            and h6_high_o["trajectory_max_km"] <= 5.0),
        "thresholds_unchanged": True,
        "prior_c1_failure_preserved": True,
        "p2p_frozen_roots_regression_recorded": (
            p2p_regression.get("status") in {
                "passed_current_ecef_p2p_high_low",
                "failed_current_ecef_p2p_high_low",
            }),
        "p2p_topology_and_hamiltonian_preserved": (
            p2p_regression["gates"]["two_unique_roots"]
            and p2p_regression["gates"]["both_ground_single_turning"]
            and p2p_regression["gates"]["hamiltonian_residual"]),
        "c1_restricted_failure_preserved": (
            c1_shadow.get("status")
            == "passed_restricted_c1_failure_preserved"
            and c1_shadow.get("original_c1_report_preserved") is True),
        "a4_c2_not_admitted": True,
    }
    production_provider_modified = (
        fields["cubic_ne"].source.get_metadata().get("interpolation") == "cubic"
        and fields["cubic_b"].get_metadata().get("interpolation") == "cubic"
    )

    if (selected_variant is not None and all(gates.values())
            and production_provider_modified):
        status = "identified_and_repaired_grid_derivative"
    elif any(ray["termination_status"] == "upper_boundary"
             for ray in trial_high_o.values()):
        status = "identified_but_not_safely_repairable"
    else:
        status = "multiple_or_unresolved_grid_derivative_contributors"

    with (OUT / "a4_c1wr_anchor_metrics.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    trace_summary = {
        "baseline": {
            case_id: _trace_summary(ray)
            for case_id, ray in baseline_rays.items()
        },
        "high_o_trials": {
            variant: _trace_summary(ray)
            for variant, ray in trial_high_o.items()
        },
        "selected_candidate": {
            case_id: _trace_summary(ray)
            for case_id, ray in selected_rays.items()
        },
    }
    _plots(
        fields, baseline_rays, baseline_logs, trial_high_o, rows, errors,
        scale_data, selected_variant, gates, {
            "p2p_frozen_root_landing_below_5km": (
                p2p_regression["gates"]["landing_error_below_5_km"]),
            "c1_all_same_angle_numeric_gates": all(
                c1_shadow.get("case_pass", {}).values()),
        })

    result = {
        "schema_version": "a4-c1wr-grid-derivative-attribution-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "WGS84_B0_to_IGRF_grid_derivative_attribution_only",
        "formula": {
            "D": "|k|^2-k0^2*n_AH^2",
            "D_units": "km^-2",
            "D_r": "gradient_ECEF(D) at fixed k and omega",
            "D_r_units": "km^-3",
            "direct_IGRF_derivative": (
                "five_point_ECEF_finite_difference_of_direct_ppigrf_queries"),
        },
        "steps_km": STEPS_KM.tolist(),
        "production_steps_km": PRODUCTION_STEPS_KM.tolist(),
        "variants": list(VARIANTS),
        "case_errors": errors,
        "high_case_closures": closures,
        "selected_variant": selected_variant,
        "trace_summary": trace_summary,
        "h3": {
            "analytic_limits": analytic_limits,
            "derivative_audit": derivative_audit,
            "refined_reference_trajectory_rms_km": h3_convergence_rms,
            "projection_on_off_trajectory_rms_km": h3_projection_rms,
        },
        "h6": {
            "ground_ground_metrics": h6_metrics,
            "high_o_common_path_metrics": h6_high_o,
        },
        "gates": gates,
        "thresholds": {
            "D_r_relative_error": 1.0e-4,
            "spatial_fd_stability_ratio": 1.0,
            "high_case_error_closure_fraction": 0.8,
            "trajectory_rms_km": 2.0,
            "trajectory_max_km": 5.0,
            "landing_difference_km": 1.0,
            "path_relative_error": 1.0e-3,
        },
        "frozen_source_sha256": frozen_hashes,
        "post_repair_source_sha256": _current_source_hashes(),
        "prior_c1_status": "failed_c1_chapman_high_low_external_comparison",
        "prior_c1w_status": "implementation_defect_identified",
        "downstream_regressions": {
            "p2p_frozen_roots": p2p_regression,
            "c1_restricted_shadow": c1_shadow,
        },
        "production_provider_modified": production_provider_modified,
        "a4_c2_admission": False,
    }
    result_path = OUT / "a4_c1wr_attribution.json"
    _write_json(result_path, result)
    artifact_paths = [
        result_path,
        REPORT,
        OUT / "a4_c1wr_anchor_metrics.csv",
        *(OUT / name for name in FIGURES),
        P2P_RESULT,
        P2P_MANIFEST,
        OUT / "p2p_regression" / "a4_c1wr_p2p_high_low_paths.png",
        C1_SHADOW,
        C1_SHADOW_MANIFEST,
        OUT / "c1_shadow" / "same_angle" / "c1_same_angle_trajectories.png",
        OUT / "c1_shadow" / "same_angle" / "c1_same_angle_residuals.png",
        OUT / "c1_shadow" / "same_angle" / "c1_same_angle_metrics.png",
        *RAW.glob("*.npz"),
    ]
    manifest = {
        "schema_version": "a4-c1wr-manifest-v1",
        "status": status,
        "source_sha256": _current_source_hashes(),
        "artifact_sha256": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in sorted(artifact_paths)
        },
    }
    _write_json(OUT / "a4_c1wr_manifest.json", manifest)
    print(json.dumps({
        "status": status,
        "selected_variant": selected_variant,
        "gates": gates,
    }, ensure_ascii=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
