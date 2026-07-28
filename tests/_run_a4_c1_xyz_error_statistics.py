"""Decompose frozen C1 same-angle differences in link-local x/y/z."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests._run_a4_c1_validation import (  # noqa: E402
    CASES, PH_SUMMARY, RAW, _load_pharlap_case, _refined_input, _resample,
)


OUT = ROOT / "output2" / "a4_pharlap" / "c1" / "c1p" / "xyz_error"
CSV_PATH = OUT / "c1p_xyz_error_statistics.csv"
JSON_PATH = OUT / "c1p_xyz_error_statistics.json"
CURVES_PATH = OUT / "c1p_xyz_path_residuals.png"
SUMMARY_PATH = OUT / "c1p_xyz_error_summary.png"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _link_coordinates(ecef: np.ndarray, tx: np.ndarray, rx: np.ndarray,
                      radius_km: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tx_unit = tx / np.linalg.norm(tx)
    rx_unit = rx / np.linalg.norm(rx)
    normal = np.cross(tx_unit, rx_unit)
    normal /= np.linalg.norm(normal)
    tangent = np.cross(normal, tx_unit)
    unit = ecef / np.linalg.norm(ecef, axis=1)[:, None]
    x = radius_km * np.arctan2(unit @ tangent, unit @ tx_unit)
    y = radius_km * np.arcsin(np.clip(unit @ normal, -1.0, 1.0))
    z = np.linalg.norm(ecef, axis=1) - radius_km
    return x, y, z


def _load_paths(case_id: str, radius_km: float, ph_cases: dict) -> tuple[dict, dict]:
    ph = _load_pharlap_case(Path(ph_cases[case_id]["raw_mat"]), radius_km)[0]
    with np.load(RAW / f"ah3d_{case_id}_raw.npz", allow_pickle=False) as source:
        ah = {
            "trajectory_ecef_km": np.asarray(
                source["repeat1_trajectory_ecef_km"], dtype=float),
            "trajectory_geom_path_km": np.asarray(
                source["repeat1_trajectory_geom_path_km"], dtype=float),
        }
    return ph, ah


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    input_path, _ = _refined_input()
    data = loadmat(input_path, squeeze_me=True, struct_as_record=False)
    radius = float(np.asarray(data["earth_radius_m"]).item()) / 1000.0

    def ecef(lat_key: str, lon_key: str, height_key: str) -> np.ndarray:
        lat = np.radians(float(np.asarray(data[lat_key]).item()))
        lon = np.radians(float(np.asarray(data[lon_key]).item()))
        r = radius + float(np.asarray(data[height_key]).item())
        return r * np.array([np.cos(lat) * np.cos(lon),
                             np.cos(lat) * np.sin(lon), np.sin(lat)])

    tx = ecef("origin_lat_deg", "origin_lon_deg", "origin_height_km")
    rx = ecef("receiver_lat_deg", "receiver_lon_deg", "receiver_height_km")
    rx_x, rx_y, rx_z = _link_coordinates(rx[None, :], tx, rx, radius)
    target = float(np.asarray(data["target_range_km"]).item())
    ph_summary = json.loads(PH_SUMMARY.read_text(encoding="utf-8"))
    ph_cases = {item["case_id"]: item for item in ph_summary["cases"]}
    u = np.linspace(0.0, 1.0, 201)
    rows, curves = [], {}

    for case_id, branch, anchor, angle in CASES:
        ph, ah = _load_paths(case_id, radius, ph_cases)
        ph_path = np.asarray(ph["trajectory_geom_path_km"], dtype=float)
        ah_path = np.asarray(ah["trajectory_geom_path_km"], dtype=float)
        ph_ecef = _resample(ph_path, ph["trajectory_ecef_km"],
                            u * (ph_path[-1] - ph_path[0]))
        ah_ecef = _resample(ah_path, ah["trajectory_ecef_km"],
                            u * (ah_path[-1] - ah_path[0]))
        ph_x, ph_y, ph_z = _link_coordinates(ph_ecef, tx, rx, radius)
        ah_x, ah_y, ah_z = _link_coordinates(ah_ecef, tx, rx, radius)
        dx, dy, dz = ah_x - ph_x, ah_y - ph_y, ah_z - ph_z
        row = {
            "case_id": case_id,
            "branch_id": branch,
            "anchor_solver": anchor,
            "elevation_deg": angle,
            "landing_dx_km": float(dx[-1]),
            "landing_dy_km": float(dy[-1]),
            "landing_dz_km": float(dz[-1]),
            "pharlap_receiver_x_residual_km": float(ph_x[-1] - rx_x[0]),
            "pharlap_receiver_y_residual_km": float(ph_y[-1] - rx_y[0]),
            "pharlap_receiver_z_residual_km": float(ph_z[-1] - rx_z[0]),
            "ah3d_receiver_x_residual_km": float(ah_x[-1] - rx_x[0]),
            "ah3d_receiver_y_residual_km": float(ah_y[-1] - rx_y[0]),
            "ah3d_receiver_z_residual_km": float(ah_z[-1] - rx_z[0]),
            "landing_abs_dx_km": float(abs(dx[-1])),
            "landing_abs_dy_km": float(abs(dy[-1])),
            "landing_abs_dz_km": float(abs(dz[-1])),
            "landing_dx_percent_tx_rx": float(100.0 * abs(dx[-1]) / target),
            "path_x_rms_km": float(np.sqrt(np.mean(dx ** 2))),
            "path_y_rms_km": float(np.sqrt(np.mean(dy ** 2))),
            "path_z_rms_km": float(np.sqrt(np.mean(dz ** 2))),
            "path_xz_rms_km": float(np.sqrt(np.mean(dx ** 2 + dz ** 2))),
            "path_xyz_rms_km": float(np.sqrt(np.mean(dx ** 2 + dy ** 2 + dz ** 2))),
        }
        if not np.isfinite(list(row.values())[3:]).all():
            raise RuntimeError(f"non-finite x-z statistic: {case_id}")
        if not np.isclose(row["path_xz_rms_km"] ** 2,
                          row["path_x_rms_km"] ** 2 + row["path_z_rms_km"] ** 2,
                          rtol=1e-12, atol=1e-12):
            raise RuntimeError(f"x-z RMS closure failed: {case_id}")
        if not np.isclose(row["path_xyz_rms_km"] ** 2,
                          row["path_x_rms_km"] ** 2 + row["path_y_rms_km"] ** 2
                          + row["path_z_rms_km"] ** 2, rtol=1e-12, atol=1e-12):
            raise RuntimeError(f"x-y-z RMS closure failed: {case_id}")
        rows.append(row)
        curves[case_id] = (dx, dy, dz)

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    JSON_PATH.write_text(json.dumps({
        "schema_version": "a4-c1p-frozen-same-angle-xyz-statistics-v1",
        "status": "completed_read_only_xyz_error_statistics",
        "semantics": {
            "comparison": "AH3D_minus_PHaRLAP_at_same_frozen_launch_angle",
            "sampling": "201_samples_at_equal_normalized_geometric_path",
            "x": "signed_great_circle_range_from_TX_toward_RX",
            "y": "signed_horizontal_cross_track_distance",
            "z": "geocentric_height_above_frozen_spherical_Earth",
            "landing": "difference_between_each_solver_ground_event_endpoint",
            "receiver_residual": "each_solver_endpoint_minus_frozen_receiver",
        },
        "source_sha256": {
            "refined_input": _sha256(input_path),
            "pharlap_summary": _sha256(PH_SUMMARY),
            **{f"pharlap_{case_id}": _sha256(Path(ph_cases[case_id]["raw_mat"]))
               for case_id, *_ in CASES},
            **{f"ah3d_{case_id}": _sha256(RAW / f"ah3d_{case_id}_raw.npz")
               for case_id, *_ in CASES},
        },
        "metrics": rows,
    }, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    figure, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for axis, (case_id, branch, anchor, _) in zip(axes.flat, CASES):
        dx, dy, dz = curves[case_id]
        axis.plot(u, dx, color="#1261a0", label="x: along-link range")
        axis.plot(u, dy, color="#16836b", label="y: horizontal cross-track")
        axis.plot(u, dz, color="#c23b22", label="z: height")
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.set_title(f"{branch}, {anchor} anchor")
        axis.set_ylabel("AH3D - PHaRLAP (km)")
        axis.grid(alpha=0.25)
    axes[1, 0].set_xlabel("Normalized geometric path u")
    axes[1, 1].set_xlabel("Normalized geometric path u")
    axes[0, 0].legend(loc="best")
    figure.suptitle("Frozen same-angle path differences in TX-RX x-y-z coordinates")
    figure.tight_layout()
    figure.savefig(CURVES_PATH, dpi=180)
    plt.close(figure)

    labels = [f"{row['branch_id'].replace('_ray', '')}\n{row['anchor_solver']}"
              for row in rows]
    index = np.arange(len(rows))
    width = 0.36
    figure, axes = plt.subplots(1, 3, figsize=(17, 5))
    axes[0].bar(index - width / 2,
                [abs(row["pharlap_receiver_x_residual_km"]) for row in rows],
                width, label="PHaRLAP |x residual|")
    axes[0].bar(index + width / 2,
                [abs(row["ah3d_receiver_x_residual_km"]) for row in rows],
                width, label="AH3D |x residual|")
    width3 = 0.24
    axes[1].bar(index - width3, [abs(row["landing_dx_km"]) for row in rows],
                width3, label="|landing dx|")
    axes[1].bar(index, [abs(row["landing_dy_km"]) for row in rows],
                width3, label="|landing dy|")
    axes[1].bar(index + width3, [abs(row["landing_dz_km"]) for row in rows],
                width3, label="|landing dz|")
    axes[2].bar(index - width3, [row["path_x_rms_km"] for row in rows],
                width3, label="path x RMS")
    axes[2].bar(index, [row["path_y_rms_km"] for row in rows],
                width3, label="path y RMS")
    axes[2].bar(index + width3, [row["path_z_rms_km"] for row in rows],
                width3, label="path z RMS")
    for axis, title in zip(axes, ("Receiver landing residual",
                                         "Cross-solver landing difference",
                                         "Cross-solver path RMS")):
        axis.set_xticks(index, labels)
        axis.set_ylabel("km")
        axis.set_title(title)
        axis.grid(axis="y", alpha=0.25)
        axis.legend()
    axes[1].set_yscale("log")
    axes[2].set_yscale("log")
    figure.suptitle("AH3D versus PHaRLAP: link-local x/y/z differences")
    figure.tight_layout()
    figure.savefig(SUMMARY_PATH, dpi=180)
    plt.close(figure)
    print(json.dumps(rows, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
