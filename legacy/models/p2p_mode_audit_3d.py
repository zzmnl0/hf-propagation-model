"""A3Rnew-P2 full-angle landing-map audit for AH3D point-to-point rays."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np


SCHEMA_VERSION = "a3rnew-p2-v1"


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return value


def _hash(value) -> str:
    payload = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _code_hash() -> str:
    digest = hashlib.sha256()
    root = Path(__file__).resolve().parent
    for name in ("p2p_mode_audit_3d.py", "point_to_point_ah_3d.py",
                 "ray_tracer_ah_3d.py", "ray_tracer_ah_3d_v2.py"):
        digest.update((root / name).read_bytes())
    return digest.hexdigest()


def _provider_metadata(n3d) -> dict:
    provider = getattr(n3d, "ne_provider", None)
    if provider is not None and hasattr(provider, "get_metadata"):
        return dict(provider.get_metadata())
    digest = hashlib.sha256()
    for value in (getattr(n3d, "x_km", []), getattr(n3d, "z_km", []),
                  getattr(n3d, "Ne_2d", [])):
        arr = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(str(arr.shape).encode("ascii"))
        digest.update(arr.tobytes())
    return {"provider_hash": digest.hexdigest(),
            "provider_name": type(n3d).__name__}


def _normalize_az(az: float) -> float:
    value = (float(az) + 180.0) % 360.0 - 180.0
    return -180.0 if np.isclose(value, 180.0) else value


def _axis(lower: float, upper: float, step: float, periodic: bool = False):
    if step <= 0.0 or upper <= lower:
        raise ValueError("invalid audit angle bounds or spacing")
    values = list(np.arange(lower, upper, step, dtype=float))
    if not periodic and (not values or values[-1] < upper):
        values.append(float(upper))
    return values


def _winding(residuals) -> int:
    vectors = np.asarray(residuals, dtype=float)
    if vectors.shape != (4, 2) or not np.all(np.isfinite(vectors)):
        return 0
    if np.any(np.linalg.norm(vectors, axis=1) < 1e-12):
        return 1
    angles = np.arctan2(vectors[:, 1], vectors[:, 0])
    closed = np.r_[angles, angles[0]]
    delta = (np.diff(closed) + np.pi) % (2.0 * np.pi) - np.pi
    return int(np.rint(np.sum(delta) / (2.0 * np.pi)))


def _cell_risks(corners, dbeta, daz, target_near_km, fold_condition):
    statuses = [item["branch_status"] for item in corners]
    valid = all(status == "ground" for status in statuses)
    residuals = np.asarray([item["residual_km"] for item in corners], dtype=float)
    triggers = []
    winding = _winding(residuals) if valid else 0
    if winding:
        triggers.append("nonzero_topological_degree")
    norms = np.linalg.norm(residuals, axis=1)
    if np.any(np.isfinite(norms)) and float(np.nanmin(norms)) <= target_near_km:
        triggers.append("target_near")
    if len(set(statuses)) > 1:
        triggers.append("branch_status_change")
    condition = np.inf
    sigma_min = 0.0
    determinants = [np.nan, np.nan]
    if valid and np.all(np.isfinite(residuals)):
        p00, p10, p11, p01 = residuals
        j0 = np.column_stack(((p10 - p00) / dbeta, (p01 - p00) / daz))
        j1 = np.column_stack(((p11 - p01) / dbeta, (p11 - p10) / daz))
        determinants = [float(np.linalg.det(j0)), float(np.linalg.det(j1))]
        jmid = 0.5 * (j0 + j1)
        singular = np.linalg.svd(jmid, compute_uv=False)
        sigma_min = float(singular[-1])
        condition = float(np.linalg.cond(jmid))
        if determinants[0] * determinants[1] <= 0.0:
            triggers.append("landing_map_fold")
        if not np.isfinite(condition) or condition >= fold_condition:
            triggers.append("small_jacobian_singular_value")
    return {
        "triggers": sorted(set(triggers)), "winding": winding,
        "condition": condition, "sigma_min": sigma_min,
        "determinants": determinants,
    }


def _atomic_json(path: Path, value: dict):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_jsonable(value), indent=2, ensure_ascii=True,
                                    allow_nan=False), encoding="utf-8")
    os.replace(temporary, path)


def _atomic_npz(path: Path, records: list[dict]):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(
            stream,
            beta=np.asarray([r["beta_deg"] for r in records]),
            az=np.asarray([r["az_deg"] for r in records]),
            x=np.asarray([r["x_land_km"] for r in records]),
            y=np.asarray([r["y_land_km"] for r in records]),
            h=np.asarray([r["h_reflect_km"] for r in records]),
            tau=np.asarray([r["tau_ms"] for r in records]),
            turning=np.asarray([r["turning_point_count"] for r in records]),
            status=np.asarray([r["branch_status"] for r in records], dtype="U64"),
        )
    os.replace(temporary, path)


def _load_npz(path: Path) -> dict:
    if not path.exists():
        return {}
    data = np.load(path, allow_pickle=False)
    cache = {}
    for index in range(len(data["beta"])):
        item = {
            "beta_deg": float(data["beta"][index]),
            "az_deg": float(data["az"][index]),
            "x_land_km": float(data["x"][index]),
            "y_land_km": float(data["y"][index]),
            "h_reflect_km": float(data["h"][index]),
            "tau_ms": float(data["tau"][index]),
            "turning_point_count": int(data["turning"][index]),
            "branch_status": str(data["status"][index]),
        }
        item["residual_km"] = [item["x_land_km"], item["y_land_km"]]
        cache[(item["beta_deg"], item["az_deg"])] = item
    return cache


def adaptive_landing_audit(tx, rx, n_ah, freq_MHz, p2p_params, t, trace_ray):
    """Return bounded root seeds from a periodic, adaptively refined 2-D map."""
    beta_bounds = tuple(map(float, p2p_params["audit_beta_bounds_deg"]))
    az_bounds = tuple(map(float, p2p_params["audit_az_bounds_deg"]))
    dbeta = float(p2p_params["audit_beta_step_deg"])
    daz = float(p2p_params["audit_az_step_deg"])
    max_level = int(p2p_params["audit_max_refinement"])
    provider = _provider_metadata(n_ah.n3d)
    frozen = {
        "schema": SCHEMA_VERSION, "branch_id": n_ah.wave_mode,
        "freq_MHz": float(freq_MHz), "tx": tx, "rx": rx, "t_s": float(t),
        "geomag": n_ah.geomag, "integration_profile": "smoke",
        "p2p": {key: value for key, value in p2p_params.items()
                if key not in ("audit_resume", "audit_progress")},
        "provider_hash": provider.get("provider_hash", ""),
        "config_hash": provider.get("config_hash", ""),
        "code_hash": _code_hash(),
    }
    reference_id = _hash(frozen)[:20]
    out = Path(p2p_params["audit_checkpoint_dir"]) / reference_id / n_ah.wave_mode
    out.mkdir(parents=True, exist_ok=True)
    metadata_path = out / "audit.json"
    grid_path = out / "landing_nodes.npz"
    if metadata_path.exists():
        old = json.loads(metadata_path.read_text(encoding="utf-8"))
        if old.get("frozen_hash") != _hash(frozen):
            raise RuntimeError("P2 checkpoint hash mismatch; refusing resume")
    cache = _load_npz(grid_path) if p2p_params.get("audit_resume", True) else {}
    for item in cache.values():
        item["residual_km"] = [item["x_land_km"] - float(rx[0]),
                               item["y_land_km"]]
    status_counts = {}
    trace_calls = 0

    def sample(beta, az):
        nonlocal trace_calls
        beta = float(beta)
        az = _normalize_az(az)
        key = (beta, az)
        if key in cache:
            return cache[key]
        ray = trace_ray(beta, az)
        item = {
            "beta_deg": beta, "az_deg": az,
            "x_land_km": float(ray.get("x_land_km", np.nan)),
            "y_land_km": float(ray.get("y_land_km", np.nan)),
            "h_reflect_km": float(ray.get("h_reflect_km", np.nan)),
            "tau_ms": float(ray.get("tau_ms", np.nan)),
            "turning_point_count": int(ray.get("turning_point_count", 0)),
            "branch_status": str(ray.get("branch_status", "unknown")),
        }
        item["residual_km"] = [item["x_land_km"] - float(rx[0]),
                               item["y_land_km"]]
        cache[key] = item
        trace_calls += 1
        status = item["branch_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
        return item

    betas = _axis(*beta_bounds, dbeta)
    azimuths = _axis(*az_bounds, daz, periodic=True)
    cells = []
    for i in range(len(betas) - 1):
        for j, az0 in enumerate(azimuths):
            az1 = azimuths[j + 1] if j + 1 < len(azimuths) else az_bounds[1]
            cells.append((betas[i], betas[i + 1], az0, az1, 0))

    seeds = []
    levels = []
    for level in range(max_level + 1):
        next_cells = []
        risk_count = 0
        for beta0, beta1, az0, az1, cell_level in cells:
            corners = [sample(beta0, az0), sample(beta1, az0),
                       sample(beta1, az1), sample(beta0, az1)]
            risks = _cell_risks(corners, beta1 - beta0, az1 - az0,
                                float(p2p_params["audit_target_near_km"]),
                                float(p2p_params["audit_fold_condition"]))
            if not risks["triggers"]:
                continue
            risk_count += 1
            seed = {
                "beta_deg": 0.5 * (beta0 + beta1),
                "az_deg": _normalize_az(0.5 * (az0 + az1)),
                "provenance": "adaptive_2d",
                "candidate_id": "{}-L{}-{}".format(n_ah.wave_mode, level, len(seeds)),
                "audit_level": level, "audit_triggers": risks["triggers"],
                "audit_winding": risks["winding"],
                "angle_bounds": ((beta0, beta1), (az0, az1)),
                "support_cells": [], "also_provenance": [],
                "covered_by_candidate_ids": [],
            }
            seeds.append(seed)
            if level < max_level:
                bm = 0.5 * (beta0 + beta1)
                am = 0.5 * (az0 + az1)
                next_cells.extend([
                    (beta0, bm, az0, am, level + 1),
                    (bm, beta1, az0, am, level + 1),
                    (bm, beta1, am, az1, level + 1),
                    (beta0, bm, am, az1, level + 1),
                ])
        levels.append({"level": level, "cells": len(cells),
                       "risk_cells": risk_count, "nodes": len(cache)})
        _atomic_npz(grid_path, list(cache.values()))
        _atomic_json(metadata_path, {
            "schema": SCHEMA_VERSION, "reference_id": reference_id,
            "frozen_hash": _hash(frozen), "frozen": frozen,
            "levels": levels, "status_counts_current_run": status_counts,
            "node_count": len(cache), "trace_calls_current_run": trace_calls,
            "completed": level == max_level,
        })
        if p2p_params.get("audit_progress", True):
            print("  p2_audit branch={} level={}/{} cells={} risk={} nodes={} new={}".format(
                n_ah.wave_mode, level, max_level, len(cells), risk_count,
                len(cache), trace_calls), flush=True)
        cells = next_cells
        if not cells:
            break
    diagnostics = {
        "reference_id": reference_id, "checkpoint_dir": str(out),
        "full_domain": True, "periodic_azimuth": True,
        "levels": levels, "node_count": len(cache), "new_trace_calls": trace_calls,
        "status_counts": status_counts, "candidate_count": len(seeds),
        "provider_hash": frozen["provider_hash"], "code_hash": frozen["code_hash"],
    }
    return seeds, diagnostics

