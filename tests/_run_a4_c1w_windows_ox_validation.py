"""Windows AH3D WGS84/OX validation and PHaRLAP behavior traceability."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg  # noqa: E402
from models.absorption_collision import complex_absorption_ecef_ray_dB  # noqa: E402
from models.dispersion_ah import (  # noqa: E402
    C_KM_S,
    appleton_hartree_roots,
    dispersion_D,
)
from models.geomagnetic_field_3d import (  # noqa: E402
    GYRO_MHZ_PER_T,
    IGRFGridGeomagneticField3D,
)
from models.ne_provider_geocentric import (  # noqa: E402
    GeocentricGridNeProvider3D,
    WGS84NeProviderAdapter,
)
from models.ray_tracer_ah_3d_ecef import (  # noqa: E402
    WGS84_COORDINATE_MODEL,
    RefractiveIndexAH3DECEF,
    dispersion_gradients_ah3d_ecef,
    hamiltonian_ah3d_ecef,
    launch_state_ah3d_ecef,
    trace_ray_ah3d_ecef,
)
from utils import (  # noqa: E402
    K_FP,
    ecef_to_geodetic,
    enu_basis_wgs84,
    geodetic_to_ecef,
)

OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1w"
RAW = OUT / "raw"
MATLAB = Path(r"D:\Matlab2024a\bin\matlab.exe")
PHARLAP = Path(r"D:\pharlap\pharlap_4.7.4")
MATLAB_SCRIPT = ROOT / "tests" / "matlab" / "a4_c1w_pharlap_wgs84_ox.m"
NM = Path(r"D:\Simply Fortran 3\mingw-w64\bin\nm.exe")

TX = (30.0, 120.0, 0.0)
RX = (32.5, 120.0, 0.0)
FREQUENCY_MHZ = 10.0
ELEVATIONS = (22.0, 22.0, 39.5, 39.5)
MODES = ("O", "X", "O", "X")
CASE_IDS = ("low_O", "low_X", "high_O", "high_X")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _chapman(height_km):
    peak = (7.0e6 / 8.98) ** 2
    u = (np.asarray(height_km, dtype=float) - 350.0) / 100.0
    exponent = 0.5 * (1.0 - u - np.exp(np.clip(-u, -700.0, 700.0)))
    return peak * np.exp(np.clip(exponent, -745.0, 700.0))


def _frozen_fields():
    latitude = np.arange(20.0, 51.0, 1.0)
    longitude = np.arange(105.0, 136.0, 1.0)
    height = np.arange(-5.0, 610.0, 5.0)
    profile = _chapman(height)
    ne = np.broadcast_to(
        profile[None, None, :],
        (latitude.size, longitude.size, height.size)).copy()
    geomagnetic = IGRFGridGeomagneticField3D.from_ppigrf(
        np.arange(20.0, 51.0, 2.0),
        np.arange(105.0, 136.0, 2.0),
        np.arange(-5.0, 610.0, 10.0),
        epoch=cfg.IRI_DT)
    source = GeocentricGridNeProvider3D(
        ne, latitude, longitude, height, cfg.RE_KM,
        source_metadata={
            "model": "Chapman_F2",
            "coordinate_use": "WGS84_geodetic_frozen_validation_grid",
        })
    return latitude, longitude, height, ne, WGS84NeProviderAdapter(source), geomagnetic


def _model(ne_provider, geomagnetic, mode, profile="wgs84_reference"):
    return RefractiveIndexAH3DECEF(
        ne_provider, earth_radius_km=cfg.RE_KM, freq_MHz=FREQUENCY_MHZ,
        wave_mode=mode, geomagnetic_provider=geomagnetic,
        coordinate_model=WGS84_COORDINATE_MODEL,
        spatial_fd_base_km=0.25, spatial_fd_scales=(0.5, 1.0, 2.0),
        integration_profile=profile)


def _trace_worker(case_id, mode, elevation, profile, projection=True):
    *_, ne_provider, geomagnetic = _frozen_fields()
    model = _model(ne_provider, geomagnetic, mode, profile)
    ray = trace_ray_ah3d_ecef(
        model, tx_lat_deg=TX[0], tx_lon_deg=TX[1],
        tx_height_km=TX[2], elevation_deg=elevation, bearing_deg=0.0,
        h_stop_km=600.0, max_path_km=5000.0, profile=profile,
        apply_dispersion_projection=projection)
    return case_id, ray


def _resample(path, values, u):
    path = np.asarray(path, dtype=float)
    values = np.asarray(values, dtype=float)
    keep = np.r_[True, np.diff(path) > 1.0e-10]
    q = path[keep]
    v = values[keep]
    target = q[0] + u * (q[-1] - q[0])
    return np.column_stack([
        np.interp(target, q, v[:, axis]) for axis in range(v.shape[1])])


def _path_metrics(ah_ray, ph_csv):
    ph = np.loadtxt(ph_csv, delimiter=",")
    ph = np.atleast_2d(ph)
    ph_ecef = ph[:, :3] / 1000.0
    ph_path = ph[:, 7] / 1000.0
    ah_ecef = np.asarray(ah_ray["trajectory_ecef_km"], dtype=float)
    ah_path = np.asarray(ah_ray["trajectory_geom_path_km"], dtype=float)
    u = np.linspace(0.0, 1.0, 201)
    ph_q = _resample(ph_path, ph_ecef, u)
    ah_q = _resample(ah_path, ah_ecef, u)
    delta = ah_q - ph_q
    xyz = np.empty_like(delta)
    for index, (point, difference) in enumerate(zip(ph_q, delta)):
        lat, lon, _ = ecef_to_geodetic(point)
        up = enu_basis_wgs84(float(lat), float(lon))[:, 2]
        if index == 0:
            tangent = ph_q[1] - ph_q[0]
        elif index == len(ph_q) - 1:
            tangent = ph_q[-1] - ph_q[-2]
        else:
            tangent = ph_q[index + 1] - ph_q[index - 1]
        tangent -= np.dot(tangent, up) * up
        tangent /= max(float(np.linalg.norm(tangent)), 1.0e-30)
        cross = np.cross(up, tangent)
        xyz[index] = [
            np.dot(difference, tangent),
            np.dot(difference, cross),
            np.dot(difference, up),
        ]
    norm = np.linalg.norm(xyz, axis=1)
    landing = xyz[-1]
    ph_group = float(ph[-1, 6] / 1000.0)
    ph_geom = float(ph[-1, 7] / 1000.0)
    return {
        "path_x_rms_km": float(np.sqrt(np.mean(xyz[:, 0] ** 2))),
        "path_y_rms_km": float(np.sqrt(np.mean(xyz[:, 1] ** 2))),
        "path_z_rms_km": float(np.sqrt(np.mean(xyz[:, 2] ** 2))),
        "path_xyz_rms_km": float(np.sqrt(np.mean(norm**2))),
        "path_max_km": float(np.max(norm)),
        "landing_x_km": float(landing[0]),
        "landing_y_km": float(landing[1]),
        "landing_z_km": float(landing[2]),
        "landing_difference_km": float(np.linalg.norm(landing)),
        "ah_geom_path_km": float(ah_ray["geom_path_km"]),
        "pharlap_geom_path_km": ph_geom,
        "geom_path_relative_error": float(
            abs(ah_ray["geom_path_km"] - ph_geom) / max(ph_geom, 1.0e-30)),
        "ah_group_path_km": float(ah_ray["group_path_km"]),
        "pharlap_group_path_km": ph_group,
        "group_path_relative_error": float(
            abs(ah_ray["group_path_km"] - ph_group) / max(ph_group, 1.0e-30)),
        "trajectory_delta_xyz_km": xyz,
        "u": u,
    }


def _write_pharlap_input(latitude, longitude, height, ne, geomagnetic):
    input_path = RAW / "a4_c1w_wgs84_ox_input.mat"
    grid = np.array([
        latitude[0], latitude[1] - latitude[0], latitude.size,
        longitude[0], longitude[1] - longitude[0], longitude.size,
        height[0], height[1] - height[0], height.size,
    ], dtype=float)
    geomag_grid = np.array([
        geomagnetic.latitude_deg[0],
        geomagnetic.latitude_deg[1] - geomagnetic.latitude_deg[0],
        geomagnetic.latitude_deg.size,
        geomagnetic.longitude_deg[0],
        geomagnetic.longitude_deg[1] - geomagnetic.longitude_deg[0],
        geomagnetic.longitude_deg.size,
        geomagnetic.height_km[0],
        geomagnetic.height_km[1] - geomagnetic.height_km[0],
        geomagnetic.height_km.size,
    ], dtype=float)
    savemat(input_path, {
        "origin_lat_deg": TX[0],
        "origin_lon_deg": TX[1],
        "origin_height_km": TX[2],
        "bearing_deg": 0.0,
        "frequency_MHz": FREQUENCY_MHZ,
        "Ne_cm3": ne / 1.0e6,
        "collision_Hz": np.zeros_like(ne),
        "iono_grid_parms": grid,
        "Bx_T": geomagnetic.b_ecef_T[..., 0],
        "By_T": geomagnetic.b_ecef_T[..., 1],
        "Bz_T": geomagnetic.b_ecef_T[..., 2],
        "geomag_grid_parms": geomag_grid,
    }, do_compression=True)
    return input_path


def _run_pharlap(input_path):
    summary = RAW / "pharlap_wgs84_ox_summary.json"
    quote = lambda value: str(value).replace("'", "''")
    expression = (
        f"addpath('{quote(MATLAB_SCRIPT.parent)}'); "
        f"a4_c1w_pharlap_wgs84_ox('{quote(input_path)}',"
        f"'{quote(RAW)}','{quote(summary)}')")
    env = {**os.environ, "NRT_NUM_THREADS": "1",
           "DIR_MODELS_REF_DAT": str(PHARLAP / "dat")}
    result = subprocess.run(
        [str(MATLAB), "-batch", expression], env=env,
        text=True, capture_output=True, errors="replace", timeout=3600)
    (RAW / "pharlap_wgs84_ox.log").write_text(
        (result.stdout or "") + "\n" + (result.stderr or ""),
        encoding="utf-8")
    if result.returncode != 0 or not summary.is_file():
        raise RuntimeError("PHaRLAP WGS84 O/X export failed")
    return json.loads(summary.read_text(encoding="utf-8"))


def _traceability():
    files = {
        "readme": PHARLAP / "00README.txt",
        "release_limit": PHARLAP / "RELEASE_LIMITATION.txt",
        "help": PHARLAP / "src" / "matlab" / "raytrace_3d.m",
        "wrapper": PHARLAP / "src" / "C" / "raytrace-3d_matlab_wrapper.c",
        "mex": PHARLAP / "mex" / "raytrace_3d.mexw64",
        "library": PHARLAP / "lib" / "win" / "libpropagation.lib",
    }
    if not all(path.is_file() for path in files.values()):
        raise FileNotFoundError("PHaRLAP traceability inputs are incomplete")
    symbols = subprocess.run(
        [str(NM), "-g", str(files["library"])], text=True,
        capture_output=True, errors="replace", timeout=120, check=True).stdout
    required = {
        "haselgrove_rhs": "hasel_rhs_3d_" in symbols,
        "electron_magnetic_derivatives": "elec_mag_pderivs_3d_" in symbols,
        "bogacki_shampine_54": "ode_bosh54_" in symbols,
    }
    return {
        "classification": "behavior_traceable_core_source_unavailable",
        "fortran_source_files_present": bool(list(PHARLAP.rglob("*.f90"))
                                             or list(PHARLAP.rglob("*.for"))),
        "binary_symbols": required,
        "sha256": {name: _sha256(path) for name, path in files.items()},
        "wrapper_calls_precompiled_raytrace_3d": True,
        "matlab_file_is_help_only": True,
        "redistribution_restricted": True,
    }


def _analytic_and_derivative_checks(ne_provider, geomagnetic):
    omega = 2.0 * np.pi * FREQUENCY_MHZ * 1.0e6
    k = np.array([omega / C_KM_S, 0.0, 0.0])
    vacuum = appleton_hartree_roots(
        omega, 0.0, k, np.array([1.0, 0.0, 0.0]), 1.2)
    b0 = appleton_hartree_roots(
        omega, 5.0e11, k, np.array([1.0, 0.0, 0.0]), 0.0)
    X = K_FP**2 * 5.0e11 / (FREQUENCY_MHZ * 1.0e6) ** 2
    Y = 1.2 / FREQUENCY_MHZ
    parallel = appleton_hartree_roots(
        omega, 5.0e11, k, np.array([1.0, 0.0, 0.0]), 1.2)
    perpendicular = appleton_hartree_roots(
        omega, 5.0e11, k, np.array([0.0, 0.0, 1.0]), 1.2)
    analytic = {
        "vacuum_O_error": float(abs(vacuum["O"] - 1.0)),
        "vacuum_X_error": float(abs(vacuum["X"] - 1.0)),
        "B0_OX_error": float(abs(b0["O"] - b0["X"])),
        "parallel_O_error": float(abs(
            parallel["O"] - (1.0 - X / (1.0 + Y)))),
        "parallel_X_error": float(abs(
            parallel["X"] - (1.0 - X / (1.0 - Y)))),
        "perpendicular_O_error": float(abs(
            perpendicular["O"] - (1.0 - X))),
    }

    model = _model(ne_provider, geomagnetic, "O")
    state = launch_state_ah3d_ecef(
        model, tx_lat_deg=30.0, tx_lon_deg=120.0, tx_height_km=300.0,
        elevation_deg=35.0, bearing_deg=5.0)
    dr, dk, dw, _ = dispersion_gradients_ah3d_ecef(state, model)

    def five_point(component, step):
        values = []
        for multiplier in (-2.0, -1.0, 1.0, 2.0):
            shifted = state.copy()
            shifted[component] += multiplier * step
            values.append(hamiltonian_ah3d_ecef(shifted, model))
        return (values[0] - 8.0 * values[1] + 8.0 * values[2]
                - values[3]) / (12.0 * step)

    dr_ref = np.array([five_point(axis, 0.25) for axis in range(3)])
    dk_ref = np.array([five_point(axis + 3, 1.0e-7) for axis in range(3)])

    def D_omega(value):
        shifted = state.copy()
        n2 = model.n2_raw(
            shifted[:3], shifted[3:6], omega=value)
        return float(dispersion_D(shifted[3:6], value, n2).real)

    hw = 1.0e3
    dw_ref = (
        D_omega(omega - 2.0 * hw) - 8.0 * D_omega(omega - hw)
        + 8.0 * D_omega(omega + hw) - D_omega(omega + 2.0 * hw)
    ) / (12.0 * hw)
    derivatives = {
        "Dr_relative_error": float(
            np.linalg.norm(dr - dr_ref) / max(np.linalg.norm(dr_ref), 1.0e-30)),
        "Dk_relative_error": float(
            np.linalg.norm(dk - dk_ref) / max(np.linalg.norm(dk_ref), 1.0e-30)),
        "Domega_relative_error": float(
            abs(dw - dw_ref) / max(abs(dw_ref), 1.0e-30)),
    }
    return analytic, derivatives


def _plot_capability(capability):
    labels = list(capability)
    values = [capability[label]["code"] for label in labels]
    colors = {2: "#2f855a", 1: "#d69e2e", 0: "#c53030"}
    figure, axis = plt.subplots(figsize=(11, 6))
    axis.barh(np.arange(len(labels)), values,
              color=[colors[value] for value in values])
    axis.set_yticks(np.arange(len(labels)), labels)
    axis.set_xticks([0, 1, 2], ["unsupported", "controlled boundary", "complete"])
    axis.set_xlim(-0.1, 2.1)
    axis.set_title("AH3D physical capability and model-boundary audit")
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_ah3d_capability_matrix.png", dpi=180)
    plt.close(figure)


def _plots(latitude, longitude, height, geomagnetic, geometry, analytic,
           derivatives, rays, metrics, traceability):
    del latitude, longitude, height, rays
    mid_lon = int(np.argmin(abs(geomagnetic.longitude_deg - 120.0)))
    field_uT = np.linalg.norm(
        geomagnetic.b_ecef_T[:, mid_lon, :, :], axis=-1) * 1.0e6
    figure, axis = plt.subplots(figsize=(9, 5))
    image = axis.pcolormesh(
        geomagnetic.latitude_deg, geomagnetic.height_km, field_uT.T,
                            shading="auto", cmap="viridis")
    figure.colorbar(image, ax=axis, label="|B| (microtesla)")
    axis.set(xlabel="WGS84 latitude (deg)", ylabel="Ellipsoid height (km)",
             title="Frozen IGRF-14 field at 120 E")
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_igrf3d_field.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].bar(list(geometry), list(geometry.values()))
    axes[0].set_yscale("log")
    axes[0].set_title("WGS84 coordinate closure")
    axes[0].tick_params(axis="x", rotation=25)
    checks = {**analytic, **derivatives}
    axes[1].bar(list(checks), list(checks.values()))
    axes[1].set_yscale("log")
    axes[1].set_title("O/X analytic limits and derivative audit")
    axes[1].tick_params(axis="x", rotation=55)
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_wgs84_ox_equation_audit.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for axis, case_id in zip(axes.flat, CASE_IDS):
        row = metrics[case_id]
        if not row["comparison_valid"]:
            axis.text(
                0.5, 0.5,
                f"Not ground-ground comparable\nPHaRLAP label {row['pharlap_ray_label']}",
                ha="center", va="center", transform=axis.transAxes)
            axis.set_title(case_id)
            axis.set_ylabel("AH3D - PHaRLAP (km)")
            axis.grid(alpha=0.25)
            continue
        xyz = row["trajectory_delta_xyz_km"]
        axis.plot(row["u"], xyz[:, 0], label="x along-track")
        axis.plot(row["u"], xyz[:, 1], label="y cross-track")
        axis.plot(row["u"], xyz[:, 2], label="z vertical")
        axis.axhline(0.0, color="black", linewidth=0.6)
        axis.set_title(case_id)
        axis.set_ylabel("AH3D - PHaRLAP (km)")
        axis.grid(alpha=0.25)
    axes[1, 0].set_xlabel("Normalized geometric path")
    axes[1, 1].set_xlabel("Normalized geometric path")
    axes[0, 0].legend()
    figure.suptitle("WGS84 same-angle O/X behavior comparison")
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_ah3d_pharlap_ox_xyz.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    labels = list(CASE_IDS)
    axes[0].bar(labels, [
        metrics[key]["path_xyz_rms_km"]
        if metrics[key]["comparison_valid"] else np.nan for key in labels],
                label="trajectory RMS")
    axes[0].axhline(2.0, color="black", linestyle="--", label="2 km gate")
    axes[0].set_ylabel("km")
    axes[0].set_title("Trajectory agreement")
    axes[0].legend()
    axes[1].bar(labels, [
        100.0 * metrics[key]["group_path_relative_error"]
        if metrics[key]["comparison_valid"] else np.nan for key in labels],
        label="group path")
    axes[1].bar(labels, [
        100.0 * metrics[key]["geom_path_relative_error"]
        if metrics[key]["comparison_valid"] else np.nan for key in labels],
        alpha=0.6, label="geometric path")
    axes[1].axhline(0.1, color="black", linestyle="--", label="0.1% gate")
    axes[1].set_ylabel("relative error (%)")
    axes[1].set_title("Path observables")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_ox_acceptance_summary.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(11, 4))
    entries = [
        ("MATLAB help", 2), ("C wrapper", 2), ("MEX behavior", 2),
        ("static symbols", 1), ("Fortran core source", 0),
        ("formula-level audit", 0),
    ]
    axis.barh([item[0] for item in entries], [item[1] for item in entries],
              color=["#2f855a" if item[1] == 2 else
                     "#d69e2e" if item[1] == 1 else "#c53030"
                     for item in entries])
    axis.set_xticks([0, 1, 2], ["unavailable", "partial", "available"])
    axis.set_title(traceability["classification"])
    figure.tight_layout()
    figure.savefig(OUT / "a4_c1w_pharlap_traceability.png", dpi=180)
    plt.close(figure)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    latitude, longitude, height, ne, ne_provider, geomagnetic = _frozen_fields()

    points = np.array([
        [0.0, 0.0, 0.0], [30.0, 120.0, 0.0],
        [-45.0, -75.0, 300.0], [89.0, 40.0, 10.0],
    ])
    reconstructed = np.array([
        ecef_to_geodetic(geodetic_to_ecef(*point)) for point in points])
    geometry = {
        "latitude_max_error_deg": float(np.max(abs(
            reconstructed[:, 0] - points[:, 0]))),
        "longitude_max_error_deg": float(np.max(abs(
            (reconstructed[:, 1] - points[:, 1] + 180.0) % 360.0 - 180.0))),
        "height_max_error_km": float(np.max(abs(
            reconstructed[:, 2] - points[:, 2]))),
        "enu_orthogonality_error": float(np.max(abs(
            enu_basis_wgs84(30.0, 120.0).T
            @ enu_basis_wgs84(30.0, 120.0) - np.eye(3)))),
    }
    analytic, derivatives = _analytic_and_derivative_checks(
        ne_provider, geomagnetic)

    rays = {}
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [
            executor.submit(
                _trace_worker, case_id, mode, elevation, "wgs84_reference")
            for case_id, mode, elevation in zip(CASE_IDS, MODES, ELEVATIONS)
        ]
        for future in futures:
            case_id, ray = future.result()
            rays[case_id] = ray
    for case_id, ray in rays.items():
        np.savez_compressed(
            RAW / f"ah3d_{case_id}.npz",
            trajectory_ecef_km=ray["trajectory_ecef_km"],
            trajectory_geom_path_km=ray["trajectory_geom_path_km"],
            trajectory_B_ecef_T=ray["trajectory_B_ecef_T"])

    with ProcessPoolExecutor(max_workers=2) as executor:
        refined_future = executor.submit(
            _trace_worker, "refined", "O", 22.0, "wgs84_refined")
        shadow_future = executor.submit(
            _trace_worker, "shadow", "O", 22.0,
            "wgs84_reference", False)
        _, research = refined_future.result()
        _, shadow = shadow_future.result()
    reference = rays["low_O"]
    u = np.linspace(0.0, 1.0, 101)
    convergence_rms = float(np.sqrt(np.mean(np.sum((
        _resample(research["trajectory_geom_path_km"],
                  research["trajectory_ecef_km"], u)
        - _resample(reference["trajectory_geom_path_km"],
                    reference["trajectory_ecef_km"], u)) ** 2, axis=1))))
    projection_shadow_rms = float(np.sqrt(np.mean(np.sum((
        _resample(shadow["trajectory_geom_path_km"],
                  shadow["trajectory_ecef_km"], u)
        - _resample(reference["trajectory_geom_path_km"],
                    reference["trajectory_ecef_km"], u)) ** 2, axis=1))))

    absorption = complex_absorption_ecef_ray_dB(
        rays["low_O"], ne_provider,
        {**cfg.UPDATE2.get("ABSORPTION_COMPLEX", {}),
         "collision_model": "d_region_exponential"},
        geomagnetic_provider=geomagnetic)

    input_path = _write_pharlap_input(
        latitude, longitude, height, ne, geomagnetic)
    pharlap_summary = _run_pharlap(input_path)
    pharlap_records = {
        item["case_id"]: item for item in pharlap_summary["records"]}
    pharlap_terminations_valid = all(
        int(item["ray_label"]) in {1, -2, -3}
        for item in pharlap_records.values())
    metrics = {}
    serializable_metrics = {}
    for case_id in CASE_IDS:
        metrics[case_id] = _path_metrics(
            rays[case_id], Path(pharlap_records[case_id]["csv"]))
        metrics[case_id]["comparison_valid"] = (
            rays[case_id]["termination_status"] == "ground"
            and int(pharlap_records[case_id]["ray_label"]) == 1)
        metrics[case_id]["pharlap_ray_label"] = int(
            pharlap_records[case_id]["ray_label"])
        serializable_metrics[case_id] = {
            key: value for key, value in metrics[case_id].items()
            if key not in {"trajectory_delta_xyz_km", "u"}}
    comparable_metrics = [
        row for row in serializable_metrics.values()
        if row["comparison_valid"]]

    traceability = _traceability()
    capability = {
        "ECEF canonical Hamiltonian state": {"code": 2},
        "DOP853 adaptive integration": {"code": 2},
        "dispersion-surface projection": {"code": 2},
        "WGS84 ground and launch geometry": {"code": 2},
        "position-dependent IGRF-14 B(r)": {"code": 2},
        "collisionless O/X real-ray dynamics": {"code": 2},
        "collisional absorption path integral": {"code": 2},
        "near-cutoff classification/rejection": {"code": 1},
        "evanescent tunneling and mode conversion": {"code": 0},
        "full-wave diffraction/scattering": {"code": 0},
        "absolute received-power calibration": {"code": 0},
    }
    _plot_capability(capability)
    _plots(latitude, longitude, height, geomagnetic, geometry, analytic,
           derivatives, rays, metrics, traceability)

    gates = {
        "wgs84_coordinate_closure": max(geometry.values()) <= 1.0e-8,
        "analytic_ox_limits": max(analytic.values()) <= 1.0e-10,
        "spatial_derivative": derivatives["Dr_relative_error"] <= 1.0e-4,
        "wavevector_derivative": derivatives["Dk_relative_error"] <= 1.0e-5,
        "frequency_derivative": derivatives["Domega_relative_error"] <= 1.0e-5,
        "integration_convergence": convergence_rms <= 0.1,
        "projection_shadow_stable": projection_shadow_rms <= 0.1,
        "all_ah3d_ground": all(
            ray["termination_status"] == "ground" for ray in rays.values()),
        "pharlap_terminations_valid": pharlap_terminations_valid,
        "local_B_absorption": (
            absorption["geomagnetic_field_semantics"]
            == "position_dependent_3d_provider"),
        "pharlap_traceable": (
            traceability["classification"]
            == "behavior_traceable_core_source_unavailable"),
        "pharlap_core_source_unavailable": (
            not traceability["fortran_source_files_present"]),
        "ox_trajectory_gate": all(
            row["path_xyz_rms_km"] <= 2.0
            and row["path_max_km"] <= 5.0
            and row["landing_difference_km"] <= 1.0
            for row in comparable_metrics) and len(comparable_metrics) >= 3,
        "ox_path_observable_gate": all(
            row["geom_path_relative_error"] <= 1.0e-3
            and row["group_path_relative_error"] <= 1.0e-3
            for row in comparable_metrics) and len(comparable_metrics) >= 3,
    }
    core_gates = (
        "wgs84_coordinate_closure", "analytic_ox_limits",
        "spatial_derivative", "wavevector_derivative",
        "frequency_derivative", "integration_convergence",
        "projection_shadow_stable", "local_B_absorption",
        "pharlap_traceable", "pharlap_core_source_unavailable",
        "pharlap_terminations_valid",
    )
    numerical_domain_rejections = [
        case_id for case_id, ray in rays.items()
        if ray["termination_status"] == "physics_domain_rejected"
        and not str(ray["failure_reason"]).startswith(
            "non-propagating AH branch:")
    ]
    physical_model_boundaries = [
        case_id for case_id, ray in rays.items()
        if ray["termination_status"] == "physics_domain_rejected"
        and str(ray["failure_reason"]).startswith(
            "non-propagating AH branch:")
    ]
    if all(gates.values()):
        final_status = "complete_collisionless_ox_hamiltonian_validated"
    elif numerical_domain_rejections or not all(
            gates[key] for key in core_gates):
        final_status = "implementation_defect_identified"
    elif physical_model_boundaries:
        final_status = "blocked_by_model_validity_boundary"
    elif not gates["ox_trajectory_gate"] or not gates["ox_path_observable_gate"]:
        final_status = "pharlap_behavior_mismatch_unresolved"
    else:
        final_status = "blocked_by_model_validity_boundary"

    result = {
        "schema_version": "a4-c1w-windows-ox-validation-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": final_status,
        "scope": (
            "cold_collisionless_geometric_optics_WGS84_position_dependent_IGRF_OX"),
        "geometry": geometry,
        "analytic_limits": analytic,
        "derivative_audit": derivatives,
        "integration": {
            "wgs84_refined_reference_trajectory_rms_km": convergence_rms,
            "projection_on_off_trajectory_rms_km": projection_shadow_rms,
        },
        "absorption": absorption,
        "pharlap_traceability": traceability,
        "pharlap_behavior_metrics": serializable_metrics,
        "ah3d_case_status": {
            case_id: {
                "termination_status": ray["termination_status"],
                "failure_reason": ray["failure_reason"],
                "turning_point_count": int(ray["turning_point_count"]),
                "geom_path_km": float(ray["geom_path_km"]),
                "group_path_km": float(ray["group_path_km"]),
                "apex_height_km": ray["apex_height_km"],
            }
            for case_id, ray in rays.items()
        },
        "capability_matrix": capability,
        "gates": gates,
        "prior_c1_status": "failed_c1_chapman_high_low_external_comparison",
        "a4_c2_admission": False,
        "model_validity_boundary_encountered": any(
            physical_model_boundaries),
        "numerical_domain_rejections": numerical_domain_rejections,
    }
    json_path = OUT / "a4_c1w_validation.json"
    json_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8")
    with (OUT / "a4_c1w_ox_metrics.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        rows = [{"case_id": case_id, **serializable_metrics[case_id]}
                for case_id in CASE_IDS]
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    sources = [
        Path(__file__), MATLAB_SCRIPT,
        ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
        ROOT / "models" / "geomagnetic_field_3d.py",
        ROOT / "models" / "dispersion_ah.py",
        input_path,
        json_path,
    ]
    manifest = {
        "schema_version": "a4-c1w-manifest-v1",
        "status": final_status,
        "sha256": {str(path.relative_to(ROOT)): _sha256(path)
                   for path in sources},
        "artifacts": sorted(
            str(path.relative_to(ROOT)) for path in OUT.glob("*")
            if path.is_file()),
    }
    (OUT / "a4_c1w_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8")
    print(json.dumps({"status": final_status, "gates": gates},
                     ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
